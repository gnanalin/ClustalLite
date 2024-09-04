# PROJECT CLUSTAL

This programm aims at reproducing the algorithm of Clustal. Clustal is a multiple alignment algorithm. In this project, we are simplifiying the algorithm.

The  three major steps are :

1. Pairwise alignment of sequences with Needleman and Wunsch algorithm.

2. Hierarchial clustering with UPGMA algorithm of the scores generated from step 1.

3. Multiple sequence alignement

## INSTALLATION

By cloning this projet :

        git clone git@github.com:gnanalin/ClustalLite.git

## CONDA INITIALIZATION

In order to have reproducible results, I used CONDA. In order to run the code, you have to create an environnement thanks to the YAML file (`environement.yml`) in the directory `ClustalLite`. 

To setup and use the environnement, you can do :

        conda env create -n project_clustallite -f environement.yml
        conda activate project_mc

You can now start the programm.