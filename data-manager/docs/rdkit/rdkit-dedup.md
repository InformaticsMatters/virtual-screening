# Job: rdkit-dedup

This describes how to run the `rdkit-dedup` job from the `comp chem` category in the `rdkit` collection.

## What the job does

This job de-duplicates molecules by comparing their canonical SMILES.
Molecule are also optionally "desalted", keeping the biggest fragment.

The jobs can handle SD files or delimited text files (e.g. tab separated) as input and output.
When using delimited text files tab separation is recommended, though other separators are supported.
When using delimited text files the molecules are read and written as SMILES.

## Implementation details

* Job implementation: [rdk_dedup.py](/rdk_dedup.py)
* Job definition: `jobs.rdkit-dedup` in [rdkit.yaml](../rdkit.yaml)

## How to run the job

### Inputs

* **Molecules**: The molecules to process, in SDF or delimited text files. If the file has a `.sdf` or `.sd` extension it is handled as a SD file, otherwise as delimited text.

### Options

* **Output file name**: The name of the output file. The format to output is determined using the file extension, `.sdf` for a SD file, otherwise delimited text.
* **Biggest fragment strategy**: The strategy for picking the biggest fragment. Option are 'hac' or 'mw' for using the 
  Heavy Atom Count or Molecular weight to identify the biggest fragment, or 'none' if you want to keep all the fragments.
* **Output has header line**: when writing delimited text files writer the first line as a header line containing the field names (1, 2).
* **Separator for text formats**:  The delimiter used for delimited text files. Tab is recommended (1).
* **Index or name to use for ID field**:  Optional field specifying the index (delimited text files, zero based) or name (SD file) for the field that will be used as the ID of the record (e.g. written as the first line of the SD file record).
* **Number of SDF records to read field names**: When reading a SD file read this number of records to determine the fields that are present (3, 4).

## Outputs

The file specified by the *Output file name* option is created. 
The type of file is determined from the file extension, `.sdf` for SD file, otherwise delimited text.

## Notes
(1). These options apply only to delimited text files and are ignored for SD files.
(2). If the input does not have a header line the field names are not know, so ones are created using the pattern *field1*, *field2* ...
(3). These options apply only to SD files and are ignored for delimited text files.
(4). There is no way of knowing all the fields that are present in a SD file without reading the data, and not all fields have to be present in every record. This means that new fields can be encountered as you process the file. Normally you can get the full list of fields by reading a small number of records, but this is not guaranteed. By default 100 records are read which is usually enough, but in rare cases you might need to read more. A new field could in theory appear in the last record!
