[![Go Report Card](https://goreportcard.com/badge/github.com/seanlaidlaw/iRODS-Downloader)](https://goreportcard.com/report/github.com/seanlaidlaw/iRODS-Downloader)
[![Maintainability](https://api.codeclimate.com/v1/badges/9d887eaebd349e28260e/maintainability)](https://codeclimate.com/github/seanlaidlaw/iRODS-Downloader/maintainability)

# iRODS-Downloader

## Installation

To download the latest binary release the following one liner can be used

```{bash}
curl -s https://api.github.com/repos/seanlaidlaw/iRODS-Downloader/releases/latest | grep browser_download_url | cut -d '"' -f 4 | wget -qi -
```

## Dependencies

### Server side Dependencies

As this is essentially a pipeline that calls other functions, quite a few
dependencies are required to be installed on the server it's run on.

- [IRODS](https://irods.org) - The data management system used to store the raw
  sequencing data
- [LSF](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=overview-lsf-introduction)
  \- the Job Scheduler used on the Sanger cluster that the wrapper uses to submit
  jobs and check job completion status

### Optional Dependencies

- Go - the language the pipeline is written in and required for compilation from
  source but _not_ if using the binaries from
  [the release](https://github.com/seanlaidlaw/iRODS-Downloader/releases)

- [viper library](https://github.com/spf13/viper) for Golang config management

## Usage

### Authentication

Before the script can be used, iRODS authentication is required. Failiure to do
so will result in the error:
`failed with error -993000 PAM_AUTH_PASSWORD_FAILED`. iRODS authentication can
be done by running `iinit` just before running the script, e.g.

```{bash}
$ iinit
Enter your current PAM password:
```

### Command line arguments

There are two required arguments, `-r` for specifying the run and `-l` for
specifying the lane. Both must be provided for the script to run properly.

```{bash}
$ ./irods_downloader -r 1234 -l 1
```

each run will download and process the samples in the working directory, so
make sure to create a specific directory before running. Additionally, if
multiple lanes need to be downloaded, the command will have to be run multiple
times, each time in a separate directories. If errors occur, rerunning the
command in the same directory will attempt to pick up where the downloader left
off thanks to the checkpoint json files that irods_downloader produces as it
goes.

### Configuration

irods_downloader will look for a configuration file named
`irods_downloader_config.yaml` to know where to look for the program
dependencies as well as what library_types it should class as RNA vs DNA. Config
file matching this filename are looked for first in the working directory (thus
allowing for project specific configs), then `$HOME/.config/` is searched, and
finally if neither location contains a config the default versions are used.

Example YAML configuration file. This is a valid config file with the default
values:

```{yaml}
bwa_align_libraries: ["GnT Picoplex"]
attribute_with_sample_name: "sample_supplier_name"
samtools_exec: "/software/CASM/modules/installs/samtools/samtools-1.11/bin/samtools"
star_exec: "/nfs/users/nfs_r/rr11/Tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR"
star_genome_dir: "/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/75/"
bwa_exec: "/software/CASM/modules/installs/bwa/bwa-0.7.17/bin/bwa"
bwa_genome_ref: "/lustre/scratch119/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"
featurecounts_exec: "/nfs/users/nfs_s/sl31/Tools/subread-2.0.1-Linux-x86_64/bin/featureCounts"
genome_annot: "/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh37d5_ERCC92/cgpRna/e75/ensembl.gtf"
```

Additionally if the default memory usage is not appropriate the config can take
the optional additional settings for each tool's ram usage:

```{yaml}
bwa_ram: "50000"
star_ram: "50000"
featurecounts_ram: "20000"
```

### Outputs

- A_iRODS_CRAM_Downloads

the downloaded CRAM and imeta files are stored here

- B_Fastq_Extraction

this is the location the gz compressed fastq files, extracted from the crams in
`A_iRODS_CRAM_Downloads`

- C_Split_by_Library_Type

here is where symlinks to the split fastqs are stored, in separate folders for
each library_type. Additionally they are named no longer by iRODS filename but
by the sample name obtained from imeta

- D_realignments

here is where the realigned bam files are output, following the library_type
separated folder structure like before. The realigned bams are sorted before
writing to disk, and are indexed in step 7 of analysis.

- E_Counts_matrix_RNA

if there are bams that have a library_type specified as RNA, the produced counts
matrix for those bams is computed and stored here.

