# Job: conformers-for-mol

This describes how to run the *conformers-for-mol* job from the *virtual screening* category in the *rdkit* collection.

## What the job does

This job generates a multiple energy minimised 3D conformers of molecules provided as SMILES input. It is typically used
to generate input for shape similarity screening using USR and it's derivatives.

The conformers are generated using RDKit's ETKDG conformer generator and energy minimised using the MMFF94 force field.
The basics are described [here](http://rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules)

The number of conformers generated is based on the rules defined in Ejeber *et al.* J. Chem. Inf. Model. 2012, 52, 1146-1158
[doi:abs/10.1021/ci2004658](https://pubs.acs.org/doi/abs/10.1021/ci2004658).
Following conformer generation the structures are energy minimised using the MMFF94 force field. Finally, conformers with
very low RMS differences are removed (keeping the lowest energy one).  Note that we use a different RMS threshold as the
default, 1.0 as opposed to 0.35 that is used in Ejeber *et al.*.

The output is a SD-file containing the conformers.

## Implementation details

* Python module: [le_conformers.py](/le_conformers.py). This module is also used by the [generate-low-energy-conformers](generate-low-energy-conformers.md) job.
* Job definition: `jobs.conformers-for-mol` in [rdkit.yaml](/data-manager/rdkit.yaml)

## How to run the job

### Inputs

* **Input molfile**: molfile containing the molecule (or use the **SMILES string for input** option).

### Options

* **SMILES string for input**: input molecule (or use the **Input molfile** input)
* **Filename for output**: Output filename
* **Number of conformers**: Number of conformers to generate. If undefined the detail number of conformers is used (see above)
* **Number of minimization cycles**: Number of energy minimization cycles (default 500)
* **RMS threshold**: Remove molecules with RMSD less than this (lower energy conformer is retained)
* **Remove Hydrogens**: Remove explict hydrogens from the output

### Outputs

Outputs are written to a SD-file.

## Related topics

* [generate-low-energy-conformers](generate-low-energy-conformers.md)