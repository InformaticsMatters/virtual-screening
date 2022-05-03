This describes how to run the *cluster-butina* job from the *comp chem* category in the *rdkit* collection.

## What the job does

This job takes a  set of molecules and clusters them using the Butina algorithm. A number of RDKit fingerprints and metrics can be used to determine the distance between the molecules. Optionally a diverse subset of those molecules can be selected.

Note that the Butina algorithm does not scale to very large datasets (more than a few thousand). Consider using the [max-min-picker job](https://discourse.squonk.it/t/job-max-min-picker/69) if you want to pick a diverse subset from a very large number of candidates.

## Implementation details

Job implementation: [/cluster_butina.py]()
Job definition: `cluster-butina jobs` in [/data-manager/rdkit.yaml]()

## How to run the job

### Inputs

**Molecules**: molecules to analyse (.smi or .sdf)

### Options

**Output file** : name for the output file (.smi or .sdf)
**Threshold**: clustering threshold (0 - 1, 1 being identical)
** Fingerprint**: fingerprint to use (maccs,  morgan2, morgan3, rdkit)
**Metric**: The similarity measure to use (braunblanquet, cosine, dice, kulczynski, mcconnaughey, rogotgoldberg, russel, sokal, tanimoto)
**Fragment method**: How to determine the biggest fragment (hac, mw)
**Output fragment**: Output the biggest fragment rather than the entire molecule
**Number to pick**: Number of molecules to pick for a diverse subset
**Exclude molecules this similar**: Exclude molecules this similar when picking a diverse subset
**Field for optimising diversity**: specify a field who's values is used to pick values for a diverse subset
**Field values ascending**: The fields values should be treated as being in ascending order (e.g low values are preferred)


#### Fragments
If the molecule is composed of multiple fragments (e.g. salts), the "biggest" fragment is used for clustering. The *Fragment method* option specifies whether to do this by heavy atom count (hac) or by molecular weight (mw). Whether the entire molecule or the biggest fragment is written to the output file is determined by the *Output fragment* option.

#### Diverse subset selection
As well as clustering the molecules options can be specified to pick a diverse subset.
In the simplest form you just specify the number of molecules you want to be output with the *Number to pick* option and the most similar molecules will be excluded giving you the number you specified.
You can also use the *Exclude molecules this similar* to specify a similarity threshold for molecules to exclude.
Also, rather than using the similarity score you can specify a field's value to use to optimise the picked molecules (e.g. pick molecules with a high activity). Do this by specifying the *Field for optimising diversity* and *Field values ascending* options.

## Outputs

The output format is specified using the file extension. If you specify *.sdf* a SD-file is output, if you specify *.smi* the molecules are output as SMILES in tab delineated text format. All input data will be written out as well as a field named *Cluster* which contains a number that specifies which cluster the molecule belongs to.

## Related topics
* [max-min-picker job](max-min-picker.md)  for picking diverse subsets from a large number of molecules