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

You need to specify a SD-file of 3D molecules that will be scored, and a reference molecule to compare to, in molfile or SD-file format. If SD-file is specified the first molecule is used as the reference.

### Inputs

* **Molecules to process**: molecules to analyse (.sdf)
* **Reference molecule**: the molecule to compare to (.mol or .sdf)

### Options

* **Output file name**: name for the SD-file output
* **Use Tanimoto distance**: use symmetric Tamimoto distance instead of symmetric distance(1)
* **Score mode**: FeatureMaps score mode (2)

Notes:
(1) If the Tanimoto option is specified the [ShapeTanimotoDist](http://rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html?highlight=shapetanimotodist#rdkit.Chem.rdShapeHelpers.ShapeTanimotoDist) function is used otherwise the [ShapeProtrudeDist](http://rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html?highlight=shapeprotrudedist#rdkit.Chem.rdShapeHelpers.ShapeProtrudeDist) function is used.

(2) See [FeatMapScoreMode](http://rdkit.org/docs/source/rdkit.Chem.FeatMaps.FeatMaps.html?highlight=featmapscoremode#rdkit.Chem.FeatMaps.FeatMaps.FeatMapScoreMode)

### Output

A SD-file is output that contains the input molecule sand properties plus the SuCOS scores.
Those scores are:
* SuCOS_Score
* SuCOS_FeatureMap_Score
* SuCOS_Tanimoto_Score (only when *Use Tanimoto distance* is true)
* SuCOS_Protrude_Score (only when *Use Tanimoto distance* is false)

All are numbers between 0 (no overlap) and 1 (perfect overlap).