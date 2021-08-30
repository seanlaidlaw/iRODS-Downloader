# iRODS-Downloader

## Usage

### Authentication

Before the script can be used, iRODS authentication is required. Failiure to do
so will result in the error:
`failed with error -993000 PAM_AUTH_PASSWORD_FAILED`. iRODS authentication can
be done by running `iinit` just before running the script, e.g.

```{bash}
$ iinit
Enter your current PAM password:
$ go run irods_downloader.go -r 1234 -l 1
```

### Command line arguments

There are two required arguments, `-r` for specifying the run and `-l` for
specifying the lane. Both must be provided for the script to run properly.

```{bash}
$ go run irods_downloader.go -r 1234 -l 1
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
