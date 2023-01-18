# Job: sucos

This describes how to run the `sucos` job from the `comp chem` category in the `xchem` collection.

## What the job does

This job takes a set of molecules and calculates SuCOS overlay scores compared to a reference molecule.
The molecules must already be aligned to the reference molecules.

SuCOS a shape and colour overlay technique and is the work of Susan Leung:
* GitHub: https://github.com/susanhleung/SuCOS
* Publication: https://doi.org/10.26434/chemrxiv.8100203.v1

It uses RDKit FeatureMaps to determine the  feature complimentarity:
* http://rdkit.blogspot.com/2017/11/using-feature-maps.html

## Implementation details

* Job implementation: [sucos.py](/sucos.py)
* Job definition: `sucos job` in [xchem.yaml](../xchem.yaml)

## How to run the job

You need to specify a SD-file of 3D molecules that will be scored, and a reference molecule to compare to, in molfile 
or SD-file format. Note that SuCOS scores the 3D molecules assuming they are already aligned. It does not generate an
alignment.

### Inputs

* **Molecules to process**: molecules to analyse (.sdf)
* **Reference molecule(s)**: the molecules to compare to (.mol or .sdf) (1)

Both sets of molecules must be 3D structures.

### Options

* **Output file name**: name for the SD-file output
* **Use Tanimoto distance**: use symmetric Tamimoto distance instead of symmetric distance (2)
* **Score mode**: FeatureMaps score mode (3)

Notes:

(1) Multiple reference molecules can be specified when using a SD-file. For instance, you might want to score molecules
against 2 fragments. In this case the scores generated are the geometric mean of the scores for each reference molecule.

(2) If the Tanimoto option is specified the [ShapeTanimotoDist](http://rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html?highlight=shapetanimotodist#rdkit.Chem.rdShapeHelpers.ShapeTanimotoDist) function is used otherwise the [ShapeProtrudeDist](http://rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html?highlight=shapeprotrudedist#rdkit.Chem.rdShapeHelpers.ShapeProtrudeDist) function is used.

(3) See [FeatMapScoreMode](http://rdkit.org/docs/source/rdkit.Chem.FeatMaps.FeatMaps.html?highlight=featmapscoremode#rdkit.Chem.FeatMaps.FeatMaps.FeatMapScoreMode)

### Output

A SD-file is output that contains the input molecule sand properties plus the SuCOS scores.
Those scores are:
* SuCOS_Score
* SuCOS_FeatureMap_Score
* SuCOS_Tanimoto_Score (only when *Use Tanimoto distance* is true)
* SuCOS_Protrude_Score (only when *Use Tanimoto distance* is false)

All are numbers between 0 (no overlap) and 1 (perfect overlap).