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

### Outputs

- 1_iRODS_CRAM_Downloads

the downloaded CRAM and imeta files are stored here

- 3_Fastq_Extraction

this is the location the gz compressed fastq files, extracted from the crams in
`1_iRODS_CRAM_Downloads`

- 4_Split_by_Library_Type

here is where symlinks to the split fastqs are stored, in separate folders for
each library_type. Additionally they are named no longer by iRODS filename but
by the sample name obtained from imeta

- 5_realignments

here is where the realigned bam files are output, following the library_type
separated folder structure like before
