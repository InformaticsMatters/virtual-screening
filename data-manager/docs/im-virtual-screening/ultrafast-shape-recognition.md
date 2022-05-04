This topic describes how to run the `ultrafast-shape-recognition` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job runs one of the variants of Ultrafast Shape Recognition (USR)
- **USR**: Ballester PJ, Richards WG (2007). DOI: [10.1002/jcc.20681](http://dx.doi.org/10.1002/jcc.20681)
- **Electroshape**: Armstrong, M. S. et al. (2010). DOI: [10.1007/s10822-010-9374-0](http://dx.doi.org/doi:10.1007/s10822-010-9374-0)
- **USRCAT**: Schreyer A, Blundell T (2012). DOI: [10.1186/1758-2946-4-27](http://dx.doi.org/10.1186/1758-2946-4-27)

These allow molecular shape similarity to be determined **without** needing to align the molecules.
The original USR algorithm considers just the shape of the molecules. Electroshape also considers electrostatics. USRCAT also considers a range of pharmacophore-style atom types.

All provide a fast way of comparing the molecular shape of a query molecule with a database of 3D conformers. Typical usage is to find molecules similar in shape to a known active molecule.

This implementation is **not** optimised for performance, though the slowest step is the conformer generation that is needed beforehand.

## Implementation details

This job uses the ODDT implementation of these three tools. Details can be found [here](https://oddt.readthedocs.io/en/latest/#molecular-shape-comparison).

* Python module: [/usr.py]()
* Job definition: `jobs.ultrafast-shape-recognition` in [/data-manager/virtual-screening.yaml]()

## How to run the job

Prior to running the job you must generate a list of molecules you want to screen, probably by running the [filter](filter.md) job 
(and possibly the [max-min-picker](max-min-picker.md) job to select a diverse subset of those).
Then you must run the [prep-enum-conf-lists](prep-enum-conf-lists.md) job, then the 
[enumerate-candidates](enumerate-candidates].md) job, then the [generate-low-energy-conformers](generate-low-energy-conformers.md)
job that will generate low energy conformers of the molecules that you want to screen. Generating the conformers can take a long time.

### Inputs

* **Query molecule**: A SD file or Molfile containing the molecule you want to use as a query. If SDF then the first record is used.
* **Molecules to screen**: Filename for conformers to be screened, typically the output of the  [filter](filter.md) job.

### Options

* **Output file name**: Name of the SD file containing the similar molecules.
* **Similarity threshold**: Similarity threshold to use (between 0 an 1). You may want to screen a small sample to determine the appropriate threshold.
* **Group by field**: Optional value for the field to group the input records by. Only the most similar molecule from each group will be reported. Typically use `std_smi` to get the most similar of all enumerated forms of each input molecule, or  `enum_smi` to group by the enumerated form (each tautomer, microstate etc.). If you do not specify a value then all forms above the similarity threshold are reported.

### Outputs

The output is a SD file containing molecules that had similarity greater than the specified threshold.
If a *Group by field* was specified only the most similar conformer for each group is output

## Related topics

* [Description of the sharded molecule system](https://discourse.squonk.it/t/the-sharded-molecule-system/88)
* [filter job](filter.md)
* [prep-enum-conf-lists job](prep-enum-conf-lists.md)
* [enumerate-candidates job](enumerate-candidates.md)
* [generate-low-energy-conformers job](generate-low-energy-conformers.md)
* [assemble-conformers job](assemble-conformers.md)