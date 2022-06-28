# Jobs: similarity-screen-smiles and similarity-screen-file

This describes how to run the `similarity-screen-smiles` or `similarity-screen-file` jobs from the `virtual screening`
category in the `rdkit` collection.

## What the job does

This job uses similarity search to filter a set of input molecules.
It does this using RDKit descriptors and metrics, a summary of which can be found 
[here](http://rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity).

Molecules with a similarity greater than a user defined threshold are written out.

Multiple query molecules can be specified, allowing molecules to be found that are similar to one or all of the queries in a user defined manner.

The  job comes in two flavours:
1. `similarity-screen-smiles`  where the queries are specified  as SMILES strings as a job parameter.
2. `similarity-screen-file` where the queries are specified in a file in SMILES format.

The second version is better suited when wanting to use a large number of queries. Although both jobs share the same implementation, they are provided as separate jobs to avoid the input options becoming overly complex.

The definition of which molecules that pass the threshold can be defined by the user along with a number of other options.

## Implementation details

* Job implementation: [screen.py](/screen.py)
* Job definition: `jobs.similarity-screen-smiles` and `jobs.similarity-screen-file`  in [rdkit.yaml](../rdkit.yaml)

## How to run the job

### Inputs

* **Molecules to screen**: a file of molecules in SMILES format.
* **Query molecules**: a file of query molecules (similarity-screen-file only)

### Options

* **Output file name**: the file name for the output
* **Query SMILES**: one or more query molecules (similarity-screen-smiles only)
* **Input file has header line**: the input file has a header line
* **Queries file has header line**: the queries file has a header line (similarity-screen-file only)
* **Output has header line**: write a header line for the output
* **Separator for input file**: separator for the inputs file
* **Separator for queries file**: separator for the queries file (similarity-screen-file only)
* **Descriptor or fingerprint type**: the type of molecular descriptor
* **Similarity metric**: The similarity metric to use
* **Similarity threshold**:  output molecules more similar to this (0 - 1)
* **Similarity score column index**: the index of the score to use to filter [1]
* **Tversky alpha**: The alpha parameter when using Tversky metric [2]
* **Tversky beta**: The beta parameter when using Tversky metric [2]
* **Number of bits for Morgan bit vector**: The fingerprint length when using Morgan bit vectors [3]

### Outputs

A file of moleciles (in SMILES format) that pass the similarity threshold is created.
This contains all the fields from the input plus one or more similarity scores. See [1] for details.

### Notes

[1] When there is only a single query molecule then there is only a single similarity score. This score is used to filter the molecules and is included in the output (as the last column).
When there are multiple query molecules this situation is more complex. Not only is there a similarity score for each query, but we need to combine those scores in a way that allows filtering. We do this by calculating the following:
- Minimum score (score_min)
- Maximum score (score_max)
- Arithmetic mean (score_amean)
- Geometric mean (score_gmean)
- Product of the scores (score_prod)

Those values are appended to the results after the individual similarity scores. e.g. for the similarity-screen-smiles  job in the case of 2 query molecules you will see the following additional fields:
- score_1
- score_2
- score_min
- score_max
- score_amean
- score_gmean
- score_prod

For the similarity-screen-file job the individual scores and the score_prod are not output as it is expected that there 
are a large number of query molecules.

When using a threshold to filter the outputs you must specify the index of the score you want to use to filter. For instance, 
if you have two query molecules and want to use the *Geometric mean* to filter then the value to use for the 
*Similarity score column index* parameter should be 5 (the first score is index zero).

[2] When using the Tversky metric you also need to provide values for the *alpha* and *beta* parameters which determine
the asymmetry of the search. Default values are alpha=1 and beta=0.  Typically these should add up to 1, but don't have to.

[3] When using Morgan fingerprints by default counts are used, but these can only be used with Tanimoto, Dice or Tversky
metrics. If you want to use the other metrics then you can instead generate a bit vector. You do this by specifying a
value for the *Number of bits for Morgan bit vector* parameter e.g. 1024. If you do not specify a value then count based
descriptors are generated. See the 
[RDKit docs](http://rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints) for more on this.