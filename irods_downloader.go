package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"strings"
)

// PIPELINE STEPS
// 0. Assess what CRAM files are being requested
// 1. Download CRAM files
// 2. Download imeta files for each cram file
// 3. Convert the CRAM files to fastq
// 4. Symlink the fastqs to two different folders depending on 'Library_type'
// 5. Align extracted fastqs with STAR or BWA depending on 'Library_type'
// 6. Samtools Quickcheck generated bams
// 7. Index realigned bam files
// 8. Sort realigned bam files
// 9. Generate counts matrix if RNA

// fileExists checks if a file exists and is not a directory before we
// try using it to prevent further errors.
func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}
func stringInSlice(a string, list []string) bool {
	for _, b := range list {
		if b == a {
			return true
		}
	}
	return false
}

func writeCheckpoint(cram_list []cram_file, step int) {
	checkpoint_file := fmt.Sprintf("checkpoint_%d.json", step)
	rankingsJson, _ := json.MarshalIndent(cram_list, "", "  ")
	err := ioutil.WriteFile(checkpoint_file, rankingsJson, 0644)
	if err != nil {
		panic(err)
	}
	log.Println(fmt.Sprintf("Checkpoint saved for step %d", step))
}

type cram_file struct {
	Filename             string
	Runid                string
	Runlane              string
	Irods_path           string
	File_exists_in_irods bool
	Cram_is_phix         bool
}

var cram_list []cram_file

func main() {

	var run string
	var lane string
	var current_step int

	// flags declaration using flag package
	flag.StringVar(&run, "r", "run", "Specify sequencing run")
	flag.StringVar(&lane, "l", "lane", "Specify sequencing lane")

	flag.Parse() // after declaring flags we need to call it
	if (strings.TrimSpace(run) == "run") || (strings.TrimSpace(lane) == "lane") {
		log.Fatalln("No lane or run argument was provided")
	}

	current_step = 0
	// if iRODS has already been polled and list of objects to download has been assesed
	// then load session information and continue
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println(fmt.Sprintf("Starting step %d", current_step))

		log.Println(fmt.Sprintf("Polling iRODS for crams associated with run: %s, and lane: %s", run, lane))
		output, err := exec.Command(
			"imeta", "qu", "-z", "seq", "-d",
			"id_run", "=", run, "and", "lane", "=", lane,
			"and", "type", "=", "cram").CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			fmt.Println("Error when running command.  Output:")
			fmt.Println(string(output))
			fmt.Printf("Got command status: %s\n", err.Error())
			return
		}

		if strings.TrimSpace(string(output)) == "No rows found" {
			log.Fatalln("No iRODS data retrieved with given lane and run")
		}

		var cram_list []cram_file
		var collection, filename string

		split_output := bytes.Split(output, []byte("----"))

		// for each cram file returned by iRODS, parse into its own object and write its run, lane,
		// and iRODS path as object metadata. Add each of these objects to the cram_list array
		log.Println("Parsing iRODS output to generate list of crams")
		for _, line := range split_output {
			collection = ""
			filename = ""
			for _, l := range strings.Split(strings.TrimSuffix(string(line), "\n"), "\n") {
				if strings.HasPrefix(l, "collection:") {
					collection = strings.TrimSpace(strings.ReplaceAll(l, "collection:", ""))
					if !strings.HasPrefix(collection, "/seq/") {
						log.Fatalf("Revieved unexpected collection '%s' for file '%s'", collection, line)
					}
				}

				if strings.HasPrefix(l, "dataObj:") {
					filename = strings.TrimSpace(strings.ReplaceAll(l, "dataObj:", ""))
					if !strings.HasSuffix(filename, ".cram") {
						log.Fatalf("Revieved unexpected filename '%s' when expecting cram", line)
					}
				}
			}

			if filename == "" || collection == "" {
				log.Fatalf("Filename or collection was empty for line: %s", line)
			}

			filename = strings.ReplaceAll(filename, ":", "")
			split_filename := strings.Split(filename, "_")
			phix_status := false
			if stringInSlice("phix.cram", split_filename) {
				phix_status = true
			}
			run_lane := strings.Split(split_filename[1], "#")[0] // this gets between _ and # which is lane
			run_lane = strings.TrimSpace(run_lane)

			cram_list = append(cram_list, cram_file{
				Filename:     filename,
				Runid:        collection,
				Runlane:      run_lane,
				Irods_path:   strings.TrimSpace(collection + "/" + filename),
				Cram_is_phix: phix_status,
			})
		}

		cram_list_len := len(cram_list)
		if cram_list_len < 1 {
			log.Fatalln("There are less than 1 items in cram list")
		}

		// for each cram in iRODS check it exists with the "ils" command and write result to object metadata
		log.Println("Verifying each iRODS cram file exists")
		for i := range cram_list {
			cram := &cram_list[i]
			output, err := exec.Command("ils", cram.Irods_path).CombinedOutput()
			if err == nil {
				cram.File_exists_in_irods = true
			} else if err != nil {
				// Display everything we got if error.
				fmt.Println("Error when running command.  Output:")
				fmt.Println(string(output))
				fmt.Printf("Got command status: %s\n", err.Error())
				return
			}
		}

		// write copy of array of cram objects to JSON file
		writeCheckpoint(cram_list, current_step)
	}
}
