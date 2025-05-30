# Job: split-sdf

This describes how to run the `split-sdf` job from the `file utils` category in the `file-utils` collection.

## What the job does

This job reads a SD-file and splits it into chunks.

## Implementation details

This job is implemented as a bash script, using the sdsplit utility from rdock.

* Job definition: `jobs.split-sdf` in [file-utils.yaml](/data-manager/file-utils.yaml)

## How to run the job

### Inputs

* **SDFile to split**: the file to split into chunks

### Options
* **Output file name**: the base name of the output files. e.g. specify `foo` and you get files
  line `foo_00001.sdf`, `foo_00002.sdf` etc.
* **Chunk size**: the number of molecules in each output file. 
