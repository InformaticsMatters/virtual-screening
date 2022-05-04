This describes how to run the `max-min-picker` job from the `virtual screening` category in the `rdkit` collection.

## What the job does

This job allows to select a diverse subset of molecules using the RDKit MaxMinPicker algorithm which picks molecules sequentially, with each new pick being as dissimilar as possible to the already picked molecules. Typically this job uses the output of the [filter](https://discourse.squonk.it/t/job-filter/68) job as its input.

For more on the MaxMinPicker see this [blog post](http://rdkit.blogspot.com/2017/11/revisting-maxminpicker.html).

## Implementation details

* Job implementation: [/max_min_picker.py]()
* Job definition: `jobs.max-min-picker` in [/data-manager/rdkit.yaml]()

## How to run the job

### Inputs

* **Molecules to pick from**:  candidate molecules in SMILES format.
* **Molecules that are already picked**: Optional set of SMILES files with candidate molecules that are already picked (the seeds)

### Options

* **Output file name**: name for the output file e.g. `diverse.smi`
* **Number of molecules to pick**: number of molecules to pick e.g. `1000`
* **Similarity threshold**: Optional similarity threshold to stop picking once reached.

## Limitations

1. Currently, only Morgan2 fingerprints with Tanimoto distance is supported.

## Related topics

* [cluster-butina job](cluster-butina.md)
* [filter job](../im-virtual-screening/filter.md)