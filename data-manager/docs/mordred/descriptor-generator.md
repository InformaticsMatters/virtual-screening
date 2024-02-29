# Job: descriptor-generator

This describes how to run the `descriptor-generator` job from the `comp chem` category in the `mordred` collection.

## What the job does

This job calculates the Mordred descriptors for molecules.
Currently this entails 1613 descriptors for 2D molecules and an additional 213 for 3D molecules.

The jobs can handle SD files or delimited text files (e.g. tab separated) as input and output.
For 3D descriptors molecules with 3D coordinates mys tbe provided in SD-file format.
When using delimited text files tab separation is recommended, though other separators are supported.
When using delimited text files the molecules are read and written as SMILES.

## Implementation details

* Job implementation: [descriptor_calc.py](/im_mordred/descriptor_calc.py)
* Job definition: `jobs.descriptor-generator` in [mordred.yaml](../mordred.yaml)
* Original Modred repo (now unsupported): https://github.com/mordred-descriptor/mordred
* Community supported repo: https://github.com/JacksonBurns/mordred-community

* This job uses the code from the community supported repo.

## How to run the job

### Inputs

* **Molecules**: The molecules to calculate, in SDF or delimited text files. If the file has a `.sdf`
  extension it is handled as a SD file, otherwise as delimited text.

### Options

* **Output file name**: The name of the output file. The format to output is determined using the file extension, `.sdf` for a SD file, otherwise delimited text.
* **Fragment method**: How to determine the largest fragment. One of Heavy Atom Count (hac), Molecular Weight (mw) or none. (1).
* **Calculate 3D descriptors**: Whether to calculate 3D descriptors (2)
* **Input has header line**: when reading delimited text files read the first line as a header line containing the field names.
* **Output file name**: The name of the output file. The format to output is determined using the file extension, `.sdf` for a SD file, otherwise delimited text.
* **Output has header line**: when writing delimited text files writer the first line as a header line containing the field names (3, 4).
* **Separator for text formats**:  The delimiter used for delimited text files. Tab is recommended (3).
* **Index or name to use for ID field**:  Optional field specifying the index (delimited text files, zero based) or name (SD file) for the field that will be used as the ID of the record (e.g. written as the first line of the SD file record).
* **Number of records to read field names**: When reading an input file read this number of records to determine the fields that are present (5, 6).

## Outputs

The file specified by the *Output file name* option is created containing all the original fields plus additional ones
for the calculated properties.
The type of file is determined from the file extension, `.sdf` for SD file, otherwise delimited text.
It is perfectly OK to use `.smi` when generating 3D descriptors, but the 3D coordinates will be lost as the molecules are
written as SMILES.

## Notes
(1). If the molecules have multiple fragments many of the descriptors will not be generated correctly, so only use the 
'none' option if you are certain that all molecules only contain a single fragment.
(2). 3D descriptors can only be generated from molecules with 3D coordinates in a SD-file 
(3). These options apply only to delimited text files and are ignored for SD files.
(4). If the input does not have a header line the field names are not know, so ones are created using the pattern *field1*, *field2* ...
(5). These options apply only to SD files and are ignored for delimited text files.
(6). There is no way of knowing all the fields that are present in a SD file without reading the data, and not all
fields have to be present in every record. This means that new fields can be encountered as you process the file. Normally you can get the full list of fields by reading a small number of records, but this is not guaranteed. By default 100 records are read which is usually enough, but in rare cases you might need to read more. A new field could in theory appear in the last record!

## Related topics

- [generate-low-energy-conformers job](../rdkit/generate-low-energy-conformers.md) can be used to generate conformers.