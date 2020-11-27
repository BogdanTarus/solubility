
# Prediction of aqueous solubility of novel drug-like compounds

## 1. Overview

The aqueous solubility of a chemical compound is its property to dissolve in water. It reflects the strength of the compound affinity to water. It is an important physico-chemical characteristic because it dictates the compound behavior in an aqueous environment, with immediate applications in medicine or industry.

It is challenging to design drug-like molecules that can efficiently reach their target in an aqueous environment. The hydrophobic active part of anti-cancer drugs needs to be pharmacokinetically compensated by groups with high aqueous solubility. Degenerative diseases, like Alzheimer’s or Parkinson diseases, are suspected to be triggered by abnormal protein-protein aggregations in the aqueous inter-cellular environment. Contrary, disruption of specific protein-protein interactions can generate genetic diseases. It is thus important to control the aqueous interaction of proteins, by disrupting their aggregation or enhancing their specific interactions, respectively.

One mean of durable development involves transition to green chemistry. A way to reach this goal is to reduce the use of organic solvents and enhance the use of water-soluble molecules. One example would be the paint industry, by employing water-soluble polymers. To increase the time-stability of the paint it is necessary to control the degree of solubility of the polymers in the initial stage, before their irreversible entanglement and/or cross-linking. Small compounds can be used as cross-linkers, increasing in this way the domain specificity of the paints.

## 2. Prepare the working environment

This systems were prepared on a macOS machine. The training of the physical descriptors based model was done on the same macOS machine. However, the deep learning SMILESX run was deployed on Google Colab. Detailed instructions on how to set the remote working environment are given at the beginig of each jupyter notebook that run on Google Colab and were downloaded [here](https://github.com/BogdanTarus/solubility/tree/master/models/2_solubility_SMILESX).

## 3. Dataset preparation

### 3.1. Individual compounds and datasets with experimental solubility values

New chemical compounds, with their corresponding experimenal aqueous solubilities, are added to the reference Delaney's dataset. Detailed description of the search procedure and the deposit place can be find [here](https://github.com/BogdanTarus/solubility/tree/master/00_database/2_new_compounds).

### 3.2. Build the molecular SMILES specification

This is a light introduction of SMILES. For an in-depth presentation of the SMILES, please see [here](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html).

SMILES (Simplified Molecular Input Line Entry System) is a string representation of chemical molecules. It is a language with a simple vocabulary, including atom and bond symbols) and a few grammar rules. 

SMILES notation consists of a series of characters containing no spaces. Hydrogen atoms may be omitted or included. There are encoding rules for atoms, bonds, branches, ring closures, and disconnections.

| Encoding rule | Structure |       Name       | SMILES |
|:-------------:|:---------:|:----------------:|:------:|
|      atom     |     S     | elemental sulfur |   [S]  |
|      atom     |    CH4    |      methane     |    C   |
|      atom     |    H2O    |       water      |    O   |
|     bonds     |   CH3CH3  |      ethane      |   CC   |
|     bonds     |    CO2    |  carbon dioxide  |  O=C=O |


The SMILE representation of a chemical structure is not necessarily unique, especially in the case of linear structures:


|       Encoding rule      |        Structure        |  Valid SMILES |
|:------------------------:|:-----------------------:|:-------------:|
|                          |                         |   C=CCC=CCO   |
| bonds: linear structures | CH2=CH-CH2-CH=CH-CH2-OH | C=C-C-C=C-C-O |
|                          |                         |   OCC=CCC=C   |

Branches are specified by enclosing them in parenthesis:

| Encoding rule |                Structure                |        SMILES        |
|:-------------:|:---------------------------------------:|:--------------------:|
|    branches   | ![](media/smiles/branches_struct-1.png) |      CC(C)C(=O)O     |
|    branches   | ![](media/smiles/branches_struct-2.png) | C=CC(CCC)C(C(C)C)CCC |


Cyclic structures are represented by breaking one bond in each ring:

|               Structure               |      SMILES      |
|:-------------------------------------:|:----------------:|
| ![](media/smiles/cyclic_struct_1.png) |     C1CCCCC1     |
| ![](media/smiles/cyclic_struct_2.png) |    C1=CC=CC=C1   |
| ![](media/smiles/cyclic_struct_3.png) | O1CCCCC1N1CCCCC1 |


These are the structures of two drug candidates and their SMILES mapping:

|              Structure             |                            SMILES                           |
|:----------------------------------:|:-----------------------------------------------------------:|
| ![](media/smiles/set_001_c001.png) |            CC(C)CCOC1=CC2=C(C=C1)C1=CC=NC(C)=C1N2           |
| ![](media/smiles/set_002_c010.png) | OC(=O)C1CCCN1C1=NC=C(C=N1)C1=CC2=NC=CC(NC3=NC=CN=C3)=C2C=C1 |


## 4. Model training

Two models have been trained in this work. The first one uses reference research of Delaney, who proposed four physical descriptors to predict the aqueous solubility of chemical compounds. The second model tries to predict the aqueous solubility directly from the compound chemical structure, without any intermediary descriptor.

### 4.1. Aqueous solubility prediction based on physical descriptors

Delaney's work to calculate the aqueous solubility, *S(mol/l)*, uses four descriptors: *logP*, *MW*, *RB*, and *AP*. The *logP* is the decimal logarithm of the octanol/water partition coefficient. It is a measure of the relative affinity of the compound for hydrophobic/aqueous solvents. It depends on the polarity of the compound. The *MW* is the molecular weight of the compound. *RB* is the number of the rotatable bonds of the compound. It is a measure of the compound entropic propensity to follow the solvent water molecules. The *AP* is the aromatic proportion, the ratio between the number of compound's aromatic atoms and the number of heavy atoms.

*logS* = C0 + C1 *LogP* + C2 *MW* + C3 *RB* + C4 *AP* 

where the *logS* is the decimal logarithm of the compound's aqueous solubility.


### 4.2. Prediction of the aqueous solubility directly from the chemical structure

This work is based on the research of Lambard and Gracheva (Guillaume Lambard and Ekaterina Gracheva 2020 Mach. Learn.: Sci. Technol. 1 025004 
). It was slightly modified to run on Google Colab.

The major inconvenient of a model that uses physical descriptors to predict the solubility is that they are empirically engineered by humans. It is thus desirable to propose methods that avoid the subjectivity of the hand-encoded features.

The SMILESX workflow propose an algorithm that takes as input the SMILES structure of chemical compounds and predicts their solubility. As a SMILES specification is represented by a string, one can use techniques specific to the Natural Language Processing (NLP) to generate models for de novo drug design. 

Another problem in training models for drug design predictions is the size of the training dataset. In general, the number of the samples in the dataset is of the order of 10^3, too small to apply modern neural network methods directly. The scarcity of the dataset is given by the experimental difficulty to get data associated to the descriptors. The SMILESX algorithm also propose a way to augment the size of the dataset by, in a first step, removing the canonicalization of the SMILES specifications. In a second step, the atoms of a given SMILES are renumbered by rotating their index correlated with a reconstruction of the correct SMILES syntax.

 However, using abstract features to construct a neural architecture makes it difficult to interpret their contribution to the observable of interest, aqueous solubility in our case. Adding an attention mechanism to the algorithm make it possible to both read deeper into the SMILES and interpret, at no extra cost, the output of the model.

## 5. Construct a web application

A flask API that predicts the aqueous solubility can be deployed on AWS/EC2. The are two forms of the API: a command-line based and a HTML based. The user should input the SMILES of the chemical structure with the unknown aqueous solubility and will get the predicted value. Detailed instructions related to the deployment process can be found [here](https://github.com/BogdanTarus/solubility/tree/master/models/1_solubility_physical-descriptors/03_flask_deployment_AWS-EC2).

## 6. Perspectives

* Make an analysis of the training dataset by inspecting how the physical descriptors are distributed among the compounds.

* Perform a classification of the compounds within the dataset based on their size, number of rotatable bonds, and aromaticity proportion. This will give as indication for what kind of compounds we can predict so solubility. 

* Add new compounds into the dataset.

* Enhance my knowledge about the SMILESX algorithm.

* Reverse design compounds based on desired MW, logS, polarity, logP and classify them as potential carcinogenic or not

* Deploy the models as serverless APIs.

* Add collaborative members to the project.

* Introduce metrics to the project core evaluation.


