# Job: enumerate-candidates

This describes how to run the `enumerate-candidates` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does
This job enumerates microstates, tautomers and undefined chiral centres of a molecules to generate a set of variants of
the molecules suitable for virtual screening.

The microstates are enumerated using [Dimorphite-DL](https://durrantlab.pitt.edu/dimorphite-dl/).
Tautomers are enumerated using [RDKit's tautomer generator](http://rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html?highlight=tautomerenumerator#rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator).
Chiral centres are enumerated using Python code using RDKit.

The input is a file containing the molecules that need enumerating.
Typically about 10 molecules are generated for each input, but this number varies considerably.

Typically the [generate-low-energy-conformers](../rdkit/generate-low-energy-conformers.md) job is run immediately after this job.

## Implementation details

* Python module: [enumerate.py](/enumerate.py)
* Job definition: `jobs.enumerate-candidates` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml)

## How to run the job

### Inputs

* **Molecules to enumerate**:  typically the  *Molecules needing enumeration* output of the 
[prep-enum-conf-lists](prep-enum-conf-lists.md) job.

### Options

* **Output file name**: the name fo the SD-file with the enumerated structures.
* **enumerateCharges**: whether to enumerate charge forms (microstates)
* **enumerateChirals**: whether to enumerate undefined chiral centres
* **enumerateTautomers**: whether to enumerate tautomers
* **combinatorial**: do combinatorial enumeration of charges and tautomers
* **minHac**, **maxHax**: min and max heavy atom counts for molecules to include
* **minPh**, **maxPh**: min and max pH range for charge enumeration
* **maxTautomers**: maximum numer of tautomers to generate
* **minCharge**, **maxCharge**: min and max absolute charge on the molecule
* **numCharges**: maximum number of atoms with charges
* **addHydrogens**: add explicit hydrogens to the generated molecules
* **tryEmbedding**: try generating a 3D molecule to check that the stereo-chemistry is sane
* **readHeader**: input file (if delineated text) contains a header line with field names
* **separator**: separator for input file (if delineated text)
* **idColumn**: column number (zero indexed, for delineated text) or field name (for SDF) with the molecule ID

### Outputs

A SD-file with the enumerated 3D structures.
Each enumerated molecule is written to the output SD-file which contains fields for the SMILES of the original molecule,
the SMILES for the enumerated molecule and a single letter code identifying which type of action has led to the creation
of the enumerated form (B = base molecule, M = microstate, T = tautomer, C = stereoisomer).
That single letter code is pulled through into the 3D conformers that will be generated next, and hence into the
docking results.

## Related topics

* [generate-low-energy-conformers job](../rdkit/generate-low-energy-conformers.md)