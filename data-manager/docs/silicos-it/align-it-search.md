# Job: align-it-search

This topic describes how to run the `align-it-search` job from the `virtual screening` category in the `silicos-it` collection.

This performs pharamacophore based 3D alignment of a set of molecules against a target molecule.
For more info [look here](http://www.silicos-it.be/software.html#from-the-original-silicos-it-toolbox)


## How to run the job

### Inputs

* **Molecule pr pharmacophore to align to**:  The molecule (.mol or .sdf) to align the candidate to or the pharmacophore
* defintion of the molecule (.phar).
* **Molecules to align**: The candidates to align (.sdf)

### Options

* **Output file base name**: Base name for output files (without extension).
* **Score cutoff**: Similarity score cuttoff (value between 0 and 1). Only output molecule more similar than this value.
* **Rank scores by**: Sort the outputs by this metric. Options are TANIMOTO, TVERSKY_DB and TVERSKY_REF.
* **Keep only this many best scores**: output only this many molecules.
* **Pharmacophore types**: The types of pharmacophore to consider (default AROM,HDON,HACC,LIPO,CHARGE).
* **Tolerance (epsilon)**: The tolerance to use (epsilon parameter)
* **Merge points**: Merge pharmacophre points that are close together.
* **Do not use normal information**: Do not use normal information.
* **Do not calculate hybrid points**: Do not calculate hybrid points
* **Inputs don't need aligning**: Input molecules are already aligned and just need scoring.
* **Use exclusion spheres during optimisation**: Use exclusion spheres during optimisation.

### Outputs

Two files are created based on the *Output file base name* option. A SD-file (.sdf extension) contains the molecules
and their properties and a tab separated values file (.tab extension) contains just the scores.

## Related topics

* [shape-it job](shape-it-search.md)
* [open3dAlign job](../rdkit/open3dalign.md)