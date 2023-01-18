# Job: sa-score

This describes how to run the `sa-score` job from the `molecular properties` category in the `rdkit` collection.

## What the job does

This job calculates a synthetic accessibility score for a set of molecules.
This is based on the work of Peter Ertl and Greg Landrum that can be found here:
https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score

That in turn is based on this paper:
Peter Ertl and Ansgar Schuffenhauer
Journal of Cheminformatics 1:8 (2009)
http://www.jcheminf.com/content/1/1/8

Scores are between 1 (easy to make) and 10 (very difficult to make).

## Implementation details

* Job implementation: [sa_scorepy](/sa_score.py)
* Job definition: `jobs.sa-score` in [rdkit.yaml](../rdkit.yaml)

## How to run the job

### Inputs

* **Molecules**: The molecules to calculate, in SDF or delimited text files. If the file has a `.sdf` or `.sd` extension it is handled as a SD file, otherwise as delimited text.

### Options