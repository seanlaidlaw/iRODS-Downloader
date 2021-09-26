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
	"reflect"
	"strings"
	"sync"
	"time"

	"github.com/spf13/viper"
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
// 8. Generate counts matrix of RNA bams

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

func bjobsIsCompleted(submitted_jobs_map map[string]string, attribute_name string, cram_list *[]cram_file) {
	// while there are still jobs in submitted_jobs_map, iterate over reading through all job outputs and only leave when all have completed successfully
	for len(submitted_jobs_map) > 0 {
		for i := range *cram_list {
			func_cram := &((*cram_list)[i])
			bjobs_output_filename := submitted_jobs_map[func_cram.Filename]

			dat, err := ioutil.ReadFile(bjobs_output_filename)
			if err == nil {
				// when job has finished (either successfully or with exit code, remove from the waiting list 'submitted_jobs_map'
				if strings.Contains(string(dat), "Terminated at") {
					delete(submitted_jobs_map, func_cram.Filename)

					// if job has finished and successfully completed then set the specified attribute_name to true
					if strings.Contains(string(dat), "Successfully completed.") {
						reflect.ValueOf(func_cram).Elem().FieldByName(attribute_name).SetBool(true)
					} else {
						log.Println(fmt.Sprintf("Error with bsub job: %s", bjobs_output_filename))
					}
				}
			}
		}
		// sleep for 5 seconds after going through every job's output before retrying
		time.Sleep(5 * time.Second)
	}
}

func quickcheck_alignments(cram_list []cram_file, i int, samtools_exec string) {
	cram := &cram_list[i]

	bam_filename := cram.Realigned_bam_path
	output, err := exec.Command(samtools_exec, "quickcheck", bam_filename).CombinedOutput()

	if err != nil {
		// Display everything we got if error.
		log.Println("Error when running command.  Output:")
		log.Println(string(output))
		log.Printf("Got command status: %s\n", err.Error())
		cram.Realigned_quickcheck_success = false
	} else {
		cram.Realigned_quickcheck_success = true
	}
	wg.Done()
}

func indexBam(cram_list []cram_file, i int, samtools_exec string) {
	cram := &cram_list[i]

	bam_filename := cram.Realigned_bam_path
	output, err := exec.Command(samtools_exec, "index", bam_filename).CombinedOutput()

	if err != nil {
		// Display everything we got if error.
		log.Println("Error when running command.  Output:")
		log.Println(string(output))
		log.Printf("Got command status: %s\n", err.Error())
		cram.Realigned_index_success = false
	} else {
		cram.Realigned_index_success = true
	}
	wg.Done()
}

type cram_file struct {
	Filename                     string
	Runid                        string
	Runlane                      string
	Irods_path                   string
	File_exists_in_irods         bool
	Cram_is_phix                 bool
	Cram_dl_path                 string
	Cram_download_success        bool
	Imeta_path                   string
	Library_type                 string
	Sample_name                  string
	Imeta_downloaded             bool
	Imeta_parsed                 bool
	Fastq_1_path                 string
	Fastq_2_path                 string
	Fastq_extracted_success      bool
	Symlinked_fq_1               string
	Symlinked_fq_2               string
	Realigned_bam_path           string
	Realigned_succesful          bool
	Realigned_quickcheck_success bool
	Realigned_index_success      bool
}

var cram_list []cram_file

var wg sync.WaitGroup

func main() {
	// we want to load a config file named "irods_downloader_config.yaml" if it exists in WD or in ~/.config
	viper.SetConfigName("irods_downloader_config")
	viper.SetConfigType("yaml")
	viper.AddConfigPath(".")              // look for config in the working directory first
	viper.AddConfigPath("$HOME/.config/") // if not found then look in .config folder

	viper.SetDefault("star_align_libraries", []string{"GnT scRNA"})
	viper.SetDefault("bwa_align_libraries", []string{"GnT Picoplex"})

	viper.SetDefault("attribute_with_sample_name", "sample_supplier_name")
	viper.SetDefault("samtools_exec", "/software/CASM/modules/installs/samtools/samtools-1.11/bin/samtools")
	viper.SetDefault("star_exec", "/nfs/users/nfs_r/rr11/Tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR")
	viper.SetDefault("star_genome_dir", "/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/75/")
	viper.SetDefault("bwa_exec", "/software/CASM/modules/installs/bwa/bwa-0.7.17/bin/bwa")
	viper.SetDefault("bwa_genome_ref", "/lustre/scratch119/casm/team78pipelines/reference/human/GRCH37d5/genome.fa")
	viper.SetDefault("featurecounts_exec", "/nfs/users/nfs_s/sl31/Tools/subread-2.0.1-Linux-x86_64/bin/featureCounts")
	viper.SetDefault("genome_annot", "/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCH37d5/star/e75/ensembl.gtf")

	// read in config file if found, else use defaults
	if err := viper.ReadInConfig(); err != nil {
		log.Fatalln("Unable to read config file")
	}

	// Config file found and successfully parsed
	//star_align_libraries := viper.GetStringSlice("star_align_libraries")
	//bwa_align_libraries := viper.GetStringSlice("bwa_align_libraries")

	attribute_with_sample_name := viper.GetString("attribute_with_sample_name")
	samtools_exec := viper.GetString("samtools_exec")
	star_exec := viper.GetString("star_exec")
	star_genome_dir := viper.GetString("star_genome_dir")
	bwa_exec := viper.GetString("bwa_exec")
	bwa_genome_ref := viper.GetString("bwa_genome_ref")
	featurecounts_exec := viper.GetString("featurecounts_exec")
	genome_annot := viper.GetString("genome_annot")

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
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
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
				log.Println("Error when running command.  Output:")
				log.Println(string(output))
				log.Printf("Got command status: %s\n", err.Error())
				return
			}
		}

		// write copy of array of cram objects to JSON file
		writeCheckpoint(cram_list, current_step)
	}

	current_step = 1
	// if CRAMS downloaded then load session information and continue
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

		// download each of the cram files
		log.Println("Starting download of iRODS CRAM files")
		cram_dl_dir := "1_iRODS_CRAM_Downloads/"
		err := os.Mkdir(cram_dl_dir, 0755)
		if err != nil {
			log.Fatal(err)
		}

		undownloaded_cram_map := make(map[string]string)
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.File_exists_in_irods {
				if !cram.Cram_is_phix {

					undownloaded_cram_map[cram.Filename] = cram_dl_dir + "/" + cram.Filename + ".o"
					cram.Cram_dl_path = cram_dl_dir + "/" + cram.Filename
					output, err := exec.Command(
						"bsub",
						"-o", cram_dl_dir+"/"+cram.Filename+".o",
						"-e", cram_dl_dir+"/"+cram.Filename+".e",
						"-R'select[mem>2000] rusage[mem=2000]'", "-M2000",
						"iget", "-K", cram.Irods_path, cram.Cram_dl_path).CombinedOutput()

					if err != nil {
						// Display everything we got if error.
						log.Println("Error when running command.  Output:")
						log.Println(string(output))
						log.Printf("Got command status: %s\n", err.Error())
						return
					}
				}
			}
		}

		// verify the cram files downloaded correctly and write download status to object metadata
		log.Println("Verifying success of downloads")
		// this updates in place the specified attribute for objects in cram_list for the
		// jobs that have finished successfully
		bjobsIsCompleted(undownloaded_cram_map, "Cram_download_success", &cram_list)

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 2
	// if imeta already downloaded then load session information and continue
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

		log.Println("Getting imeta for each downloaded cram")
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Cram_download_success {

				cmd := exec.Command("imeta", "ls", "-d", cram.Irods_path)

				cram.Imeta_path = cram.Cram_dl_path + ".imeta"
				imeta_file, err := os.Create(cram.Imeta_path)
				if err != nil {
					log.Fatal(err)
				}
				defer imeta_file.Close()

				// Send stdout to the outfile. cmd.Stdout will take any io.Writer.
				cmd.Stdout = imeta_file
				err = cmd.Run()
				if err != nil {
					// Display everything we got if error.
					log.Println("Error when running command.  Output:")
					log.Printf("Got command status: %s\n", err.Error())
					return
				}
				cram.Imeta_downloaded = true
			}
		}

		log.Println("Parsing imeta file to obtain library_type and sample name")
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Cram_download_success {

				library_type := ""
				sample_name := ""

				imeta, _ := ioutil.ReadFile(cram.Imeta_path)
				split_imeta := bytes.Split(imeta, []byte("----"))
				for _, line := range split_imeta {
					if bytes.Contains(line, []byte("attribute: library_type")) {
						line_split := strings.Split(strings.TrimSuffix(string(line), "\n"), "\n")
						for _, l := range line_split {
							if strings.Contains(l, "value:") {
								library_type = strings.ReplaceAll(l, "value:", "")
								library_type = strings.TrimSpace(library_type)
								cram.Library_type = library_type
								break
							}
						}
						if library_type != "" {
							break
						}
					}
				}
				for _, line := range split_imeta {
					if bytes.Contains(line, []byte("attribute: "+attribute_with_sample_name)) {
						line_split := strings.Split(strings.TrimSuffix(string(line), "\n"), "\n")
						for _, l := range line_split {
							if strings.Contains(l, "value:") {
								sample_name = strings.ReplaceAll(l, "value:", "")
								sample_name = strings.TrimSpace(sample_name)
								cram.Sample_name = sample_name
								cram.Imeta_parsed = true
								break
							}
						}
						// as this loop goes line by line, this break means it stops when the first sample_name is found.
						if sample_name != "" {
							break
						}
					}
				}
			}
		}

		log.Println("Checking there are no duplicate sample names in parsed imeta information")
		sample_names_map := make(map[string][]string)
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Library_type != "" {
				sample_names_map[cram.Library_type] = []string{}
			}
		}

		if len(sample_names_map) == 0 {
			log.Fatalln("There are no library_type information for samples")
		}

		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Library_type != "" {
				if stringInSlice(cram.Sample_name, sample_names_map[cram.Library_type]) {
					log.Printf("Duplicate sample_names found for %s", cram.Sample_name)
					log.Fatalln("There are duplicate values in sample_names, double check your choice of 'attribute_with_sample_name'")
				}
				sample_names_map[cram.Library_type] = append(sample_names_map[cram.Library_type], cram.Sample_name)
			}
		}

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 3
	// if fastq have already been split then load checkpoint instead of rerunning
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

		log.Println("Extracting fastq from downloaded crams")
		err := os.Mkdir("3_Fastq_Extraction", 0755)
		if err != nil {
			log.Fatal(err)
		}

		// extract fastq from downloaded cram files
		fastq_cram_map := make(map[string]string)
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Imeta_parsed {

				fastq_cram_map[cram.Filename] = "3_Fastq_Extraction/3_cram_to_fastq_" + cram.Filename + ".o"
				fq_filename := strings.ReplaceAll(cram.Filename, ".cram", "")
				cram.Fastq_1_path = "3_Fastq_Extraction/" + fq_filename + ".1.fq.gz"
				cram.Fastq_2_path = "3_Fastq_Extraction/" + fq_filename + ".2.fq.gz"

				output, err := exec.Command(
					"bsub",
					"-o", "3_Fastq_Extraction/3_cram_to_fastq_"+cram.Filename+".o",
					"-e", "3_Fastq_Extraction/3_cram_to_fastq_"+cram.Filename+".e",
					"-R'select[mem>2000] rusage[mem=2000]'", "-M2000",
					"-n", "4",
					samtools_exec, "fastq", "-c", "7", "-@", "4",
					"-1", cram.Fastq_1_path,
					"-2", cram.Fastq_2_path,
					"-0", "/dev/null",
					"-s", "/dev/null",
					"-n", cram.Cram_dl_path).CombinedOutput()

				if err != nil {
					// Display everything we got if error.
					log.Println("Error when running command.  Output:")
					log.Println(string(output))
					log.Printf("Got command status: %s\n", err.Error())
					return
				}
			}
		}

		// verify extracting the crams into fastq finished successfully
		log.Println("Verifying success of extracting fastq")

		bjobsIsCompleted(fastq_cram_map, "Fastq_extracted_success", &cram_list)

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 4
	// if symlinks already created then load session information and continue
	// if fastq have already been split then load checkpoint instead of rerunning
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {

		log.Println("Symlinking fastq into different folders based on library_type")
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Fastq_extracted_success && cram.Library_type != "" {
				lib_type_dir := strings.ReplaceAll(cram.Library_type, " ", "_")
				lib_type_dir = "4_Split_by_Library_Type/" + lib_type_dir
				_ = os.MkdirAll(lib_type_dir, 0755)
				os.Symlink("../../"+cram.Fastq_1_path, lib_type_dir+"/"+cram.Sample_name+".1.fq.gz")
				os.Symlink("../../"+cram.Fastq_2_path, lib_type_dir+"/"+cram.Sample_name+".2.fq.gz")
				cram.Symlinked_fq_1 = lib_type_dir + "/" + cram.Sample_name + ".1.fq.gz"
				cram.Symlinked_fq_2 = lib_type_dir + "/" + cram.Sample_name + ".2.fq.gz"
			}
		}

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 5
	// if alignments already performed then load session information and continue
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {

		_ = os.MkdirAll("5_realignments/", 0755)
		log.Println("Running alignments between extracted fastq and specified reference")
		realignment_map := make(map[string]string)
		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Symlinked_fq_1 != "" && cram.Symlinked_fq_2 != "" {

				out_folder := "5_realignments/" + strings.ReplaceAll(cram.Library_type, " ", "_") + "/"
				_ = os.MkdirAll(out_folder, 0755)

				bam_output := out_folder + cram.Sample_name + ".bam"
				job_out := out_folder + "/5_realignement_RNA_" + cram.Sample_name + ".o"
				job_err := out_folder + "/5_realignement_RNA_" + cram.Sample_name + ".e"

				if stringInSlice(cram.Library_type, star_align_libraries) {
					output, err := exec.Command(
						"bsub",
						"-o", job_out,
						"-e", job_err,
						"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
						"-n", "10",
						star_exec, "--runThreadN", "10",
						"--outSAMattributes", "NH", "HI", "NM", "MD",
						"--limitBAMsortRAM", "31532137230",
						"--outSAMtype", "BAM", "SortedByCoordinate",
						"--genomeDir", star_genome_dir,
						"--readFilesCommand", "zcat",
						"--outFileNamePrefix", out_folder+"/"+cram.Filename,
						"--readFilesIn", cram.Symlinked_fq_1, cram.Symlinked_fq_2,
						"--outStd", "BAM_SortedByCoordinate",
						"|", samtools_exec, "sort", "-@3", "-l7", "-o", bam_output).CombinedOutput()

					if err != nil {
						// Display everything we got if error.
						log.Println("Error when running command.  Output:")
						log.Println(string(output))
						log.Printf("Got command status: %s\n", err.Error())
						return
					}

					cram.Realigned_bam_path = bam_output
					realignment_map[cram.Filename] = job_out

				} else if stringInSlice(cram.Library_type, bwa_align_libraries) {
					output, err := exec.Command(
						"bsub",
						"-o", job_out,
						"-e", job_err,
						"-R'select[mem>50000] rusage[mem=50000]'", "-M50000",
						"-n", "10",
						bwa_exec, "mem", "-t", "10",
						bwa_genome_ref,
						cram.Symlinked_fq_1,
						cram.Symlinked_fq_2,
						"|", samtools_exec, "sort", "-@3", "-l7", "-o", bam_output).CombinedOutput()

					if err != nil {
						// Display everything we got if error.
						log.Println("Error when running command.  Output:")
						log.Println(string(output))
						log.Printf("Got command status: %s\n", err.Error())
						return
					}

					cram.Realigned_bam_path = bam_output
					realignment_map[cram.Filename] = job_out
				}
			}
		}

		// check on alignment jobs until they have finished
		bjobsIsCompleted(realignment_map, "Realigned_succesful", &cram_list)

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 6
	// if quickcheck has already been performed then load session information and continue
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println("Running samtools quickcheck on completed bams")

		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Realigned_succesful {
				wg.Add(1)
				go quickcheck_alignments(cram_list, i, samtools_exec)
			}
		}
		wg.Wait() // wait until all quickcheck processes have finished

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 7
	// if quickcheck has already been performed then load session information and continue
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println("Indexing quickchecked bams")

		for i := range cram_list {
			cram := &cram_list[i]
			if cram.Realigned_quickcheck_success {
				wg.Add(1)
				indexBam(cram_list, i, samtools_exec)
			}
		}
		wg.Wait() // wait until all quickcheck processes have finished

		writeCheckpoint(cram_list, current_step)
	}

	current_step = 8
	// if counts matrix has already been built then load session info
	if fileExists(fmt.Sprintf("checkpoint_%d.json", current_step)) {
		jsonFile, err := os.Open(fmt.Sprintf("checkpoint_%d.json", current_step))
		byteValue, _ := ioutil.ReadAll(jsonFile)
		err = json.Unmarshal([]byte(byteValue), &cram_list)
		if err != nil {
			panic(err)
		}

		log.Println(fmt.Sprintf("Checkpoint exists for step %d, loading progress", current_step))

	} else {
		log.Println("Running featurecounts on completed RNA bams")
		var rna_bams_featurecounts_input []string

		for i := range cram_list {
			cram := &cram_list[i]
			// if quickcheck worked then add its realigned and sorted bam path to list of bams to include in counts matrix
			if cram.Realigned_quickcheck_success {
				if stringInSlice(cram.Library_type, star_align_libraries) {
					rna_bams_featurecounts_input = append(rna_bams_featurecounts_input, cram.Realigned_bam_path)

				}
			}
		}

		if len(rna_bams_featurecounts_input) < 1 {
			log.Fatalln("Less than  1 bams in RNA category, not enough for featurecounts, aborting.")
		}

		err := os.Mkdir("6_Counts_matrix_RNA", 0755)
		if err != nil {
			log.Fatal(err)
		}

		matrix_out := "6_Counts_matrix_RNA/featurecounts_matrix.tsv"
		job_out := "6_Counts_matrix_RNA/featurecounts_run.o"
		job_err := "6_Counts_matrix_RNA/featurecounts_run.e"

		featureCountsCmd := []string{
			"-o", job_out,
			"-e", job_err,
			"-R'select[mem>20000] rusage[mem=20000]'", "-M20000",
			"-n", "14",
			featurecounts_exec,
			"-Q", "30",
			"-p",
			"-t", "exon",
			"-g", "gene_name",
			"-F", "GTF",
			"-a", genome_annot,
			"-o", matrix_out}

		// append bam paths to end of command options, as this is what featureCounts expects
		featureCountsCmd = append(featureCountsCmd, rna_bams_featurecounts_input...)

		output, err := exec.Command("bsub", featureCountsCmd...).CombinedOutput()

		if err != nil {
			// Display everything we got if error.
			log.Println("Error when running command.  Output:")
			log.Println(string(output))
			log.Printf("Got command status: %s\n", err.Error())
			return
		}

		// wait for featurecounts job to finish
		for {
			dat, err := ioutil.ReadFile(job_out)
			if err == nil {
				if strings.Contains(string(dat), "Terminated at") {
					break
				}
			}
			time.Sleep(5 * time.Second)
		}

		// if featurecounts exited successfuly write new checkpoint file
		// this doesn't have any new information but its presence will indicate not to repeat the featurecounts step
		dat, err := ioutil.ReadFile(job_out)
		if err == nil {
			if strings.Contains(string(dat), "Successfully completed.") {
				writeCheckpoint(cram_list, current_step)
			}
		}
	}
}
