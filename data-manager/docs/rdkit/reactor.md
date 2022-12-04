# Job: reactor

This describes how to run the `reactor` job from the `molecular libraries` category in the `rdkit` collection.

## What the job does

This job takes one or more sets of reactants and one or more reaction definitions and generates reaction products.
A simple example would be to use a set or carboxylic acids and a set of amines to generate amide products. In this
simplest case you would specify a file with the carboxylic acids in SMILES format as one reactants file, a file with
the amines as a second file and a reaction in reaction SMARTS format. You can optionally specify multiple reaction
SMARTS if you have variations of the reaction e.g. one for primary amines and another for secondary amines.

## Implementation details

* Job implementation: [reactor.py](/reactor.py)
* Job definition: `reactor job` in [rdkit.yaml](/data-manager/rdkit.yaml)

## How to run the job

### Inputs

* **Reactant 1**: file with 1st reactants (.smi) 
* **Reactant 2**: file with 2nd reactants (.smi) 
* **Reactant 3**: file with 3rd reactants (.smi)

Only the first reactant file is manadatory.

### Options

* **Reactions**: One or more reaction SMARTS
* **Output file**: name for the output file (.smi)
* **Reactants have header line**: skip the header line in the reactants files
* **Output has header line**: Write a header line in the output file


## Outputs

The molecules are output as SMILES in tab delineated text format. 
The format is as follows:

product index reactant1 reactant2 ...

The fields are tab separated.
The product and reactants are in SMILES format.
The index is the numerical index of the combination of reactants. Some reactants can
generate multiple products (e.g. if there are multiple reactive centres).


## Challenges

Two aspects need careful consideration, and are somewhat inter-dependent.

### 1. Choice of reactants

Only compatible reactants will be reacted to generate a product. Hence, you can use
a brute force approach and use a whole building blocks library as all your inputs.
However, this is very inefficient as a lot of time will be spent reacting unreactive reactants.

It is better to generate sets of reactants that correspond to "compatible" reactants. e.g. if your reaction is for 
acid halides then generate a file that just contains acid halides. The "MolDB" functionality that is currently being
rolled out will provide smart ways of doing this, but right now you are best to generate those reactant files and 
upload them (which will always be a good solution anyway). Sharing those files as "Datasets" is also a good idea.

### 2. Reaction SMARTS specification

Understanding SMARTS expressions is hard! Writing reaction SMARTS is even harder. Do you just want to react primary 
amines or secondary as well? Do you want to exclude amines that are part of an amide?
The [RDKit docs](http://rdkit.org/docs/source/rdkit.Chem.AllChem.html#rdkit.Chem.AllChem.EnumerateLibraryFromReaction)
gives an example for amide formation as `[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]`. That certainly enumerates classical
amide formation from acids and amines, but also enumerates other things that might not be expected (e.g. using the N
of an azide as the "amine"), and it doesn't handle acid halides as the "acid" part.

Yes, you can generate a more sophisticated reaction SMARTS that might handle these aspects, but defining a single SMARTS
that covers every subtlety is challenging if not impossible. Certainly not something you would expect a typical chemist
to be able to do.

Hence, we try to provide some examples that can be used (see next section).

And, that's why we allow multiple reaction SMARTS to be specified so that you can provide multiple variations of your 
reaction rather than try to capture them all in one super-complex reaction SMARTS expression. This is modeled on 
earlier work in the Squonk Computational Notebook that was defined by Anthony Bradley. See
[here](https://github.com/InformaticsMatters/pipelines/blob/b0830631bc77745ee5c71df2ea2c624124594802/src/python/pipelines/rdkit/poised_filter.py#L65-L114)
for examples.

## Examples

Here we try to exemplify a number of "standard" reactions. We welcome input into improving these and note that they are mostly untested.
This need to be community driven. 

| Reaction            | Reaction SMARTS                                                                                                                   | Reactants                                      | 
|---------------------|-----------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------|
| Amide formation     | `[C:1](=[O:2])-[OD1,Cl,Br].[N!H0:3]>>[C:1](=[O:2])[N:3]`                                                                          | 1. Carboxylic acids, Acid halides<br>2. Amines |
| SNAr (1)            | `[#7H2:1].[c:2]Br>>[c:2][PH1:1]`                                                                                                  |                                                |
|                     | `[NH1:1].[c:2]Br>>[c:2][PH0:1]++[#15:1]>>[#7:1]`                                                                                  |                                                |
|                     | `[c:2]-[Cl,Br,I].[#7:1]>>[c:2]-[#7:1]`                                                                                            |                                                |
| Urea                | `[#7H2:1].[#7:2]>>[#7H1:1]C([#7:2])=O`                                                                                            |                                                |
|                     | `[NH1:1].[#7:2]>>[#7H0:1]C([#7:2])=O`                                                                                             |                                                |
|                     | `[#7:2]=C=O.[#7H2:1]>>[#7:1]-[#6](-[#7:2])=O`                                                                                     |                                                |
|                     | `[#7:2]=C=O.[#7;AH1:1]>>[#7:1]-[#6](-[#7:2])=O`                                                                                   |                                                |
| Suzuki Coupling     | `[#6;a:2]-[#35,#53].[#6;a:1]-[#5](-[#8])-[#8]>>[#6;a:1]-[#6;a:2]`                                                                 |                                                |
|                     | `[#6;a:1]-[#5](-[#8])-[#8].[#6;a:2]-[#35,#53]>>[#6;a:1]-[#6;a:2]`                                                                 |                                                |
| Sonogashira         | `C:1]#[C:2].[#6;a:3]-[#35,#53]>>[#6;a:3][C:1]#[C:2]`                                                                              |                                                |
|                     | `[#6;a:1]-[#35,#53].[C:2]#[C:3]>>[#6;a:1][C:2]#[C:3]`                                                                             |                                                |
| Sulfonamide         | `[#17,#35,#53][S:4](=[O:3])=[O:2].[#7:1]>>[#7:1][S:4](=[O:3])=[O:2]`                                                              |                                                |
|                     | `[#7H2:1].Cl[S:2](=[O:3])=[O:4]>>[#7:1][S:2](=[O:3])=[O:4]`                                                                       |                                                |
|                     | `[#7H1:1].Cl[S:2](=[O:3])=[O:4]>>[#7:1][S:2](=[O:3])=[O:4]`                                                                       |                                                |
| Reductive Amination | `[#6H1:2]=O.[#7:1]>>[#6:2]-[P:1]`                                                                                                 |                                                |
|                     | `[#6:3]-[#6:2](-[#6:4])=O.[#7:1]>>[#6:3]-[#6:2](-[#6:4])-[P:1]++[P:1]>>[N:1]`                                                     |                                                |
|                     | `[#7;AH2:1].[#6:4]-[#6:2](-[#6:3])=O>>[#6:4]-[P:2](-[#6:3])-[#7:1]`                                                               |                                                |
|                     | `[#7AH1:1].[#6:2]-[#6:3](-[#6:4])=O>>[#6:2]-[P:3](-[#6:4])-[#7:1]`                                                                |                                                |
|                     | `[#7AH2:1].[#6:4]-[#6H1:2]=O>>[#6:4]-[P:2]-[#7:1]`                                                                                |                                                |
|                     | `[#7AH1:1].[#6:2]-[#6H1:3]=O>>[#6:2]-[P:3]-[#7:1]++[P:1]>>[C:1]`                                                                  |                                                |
| N-Alkylation        | `[#7H1:1].[#6;A:2]Br>>[#7:1]-[#15:2]++[P:1]>>[C:1]`                                                                               |                                                |
|                     | `[C:2][Cl,Br,I].[#7H1:1]>>[#7:1]-[#15:2]++[P:1]>>[C:1]`                                                                           |                                                |
| Ether Coupling      | `[#8H1:1].[#6;A:2]Br>>[#6:2]-[#8:1]`                                                                                              |                                                |
|                     | `[OH1:1].[n:3][#6;a:2]Br>>[O:1]-[#6:2][n:3]`                                                                                      |                                                |
|                     | `[#6;A:2][#17,#35,#53].[#8H1:1]>>[#6:2]-[#8:1]`                                                                                   |                                                |
|                     | `[n:3][c:2][#17,#35,#53].[#8H1:1]>>[n:3][c:2][#8:1]`                                                                              |                                                |
| Ester Coupling      | `[Cl,Br,I][C:2]=[O:3].[#8:1]>>[#8:1][C:2]=[O:3]`                                                                                  |                                                |
|                     | `[OH1][C:2]=[O:3].[#8:1]>>[#8:1][C:2]=[O:3]`                                                                                      |                                                |
| Benzimidazole       | `Cl[#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2]>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                     |                                                |
|                     | `[OH][#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2]>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                   |                                                |
|                     | `[N:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[N:2].Cl[#6:9]=O>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                       |                                                |
|                     | `[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2].Cl[#6:9]=O>>[nH:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                     |                                                |
| Triazole            | `[C:3]#[C:4][#6:6].[#7:1][#7:5]=[#7:2]>>[#7:1]1[#6H:3]=[#6:4]([#6:6])[#7:2]=[#7:5]1`                                              |                                                |
|                     | `[#7:1]=[#7+:2]=[#7-:3]>>[#7:1][#7:2]=[#7:3]++[#7:1][#7:5]=[#7:2].[C:3]#[C:4][#6:6]>>[#7:1]1[#6H:3]=[#6:4]([#6:6])[#7:2]=[#7:5]1` |                                                |
| Benzoxazole         | `Cl[#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2]>>[n:1]1[c:9][o:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                      |                                                |
|                     | `[OH][#6:9]=O.[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2]>>[n:1]1[c:9][o:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                    |                                                |
|                     | `[OH:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[NH2:2].Cl[#6:9]=O>>[o:1]1[c:9][n:2][c:4]2[c:5][c:6][c:7][c:8][c:3]12`                     |  

Notes:
1. SNAr: Nucleophilic Aromatic Substitution