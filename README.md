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
