>Assan:
Bogdan pls try to to organize it with 1. 2. etc.
the style section did not work on my end, just a FYI
- pls add a description on the tab on the right corner thx.
- pls add pictures ( to explain SMiles)
- pls add as well resources for ppl wishing to discover what smiles are (links)
>


# Prediction of aqueous solubility of novel drug-like compounds
> Assan:  Why is it relevant to you , to the community etc.

## 1. Overview

The aqueous solubility of a chemical compound is its property to dissolve in water. It reflects the strength of the compound affinity to water. It is an important physico-chemical characteristic because it dictates the compound behavior in an aqueous environment, with immediate applications in medicine or industry.

It is challenging to design drug-like molecules that can efficiently reach their target in an aqueous environment. The hydrophobic active part of anti-cancer drugs needs to be pharmacokinetically compensated by groups with high aqueous solubility. Degenerative diseases, like Alzheimer’s or Parkinson diseases, are suspected to be triggered by abnormal protein-protein aggregations in the aqueous inter-cellular environment. Contrary, disruption of specific protein-protein interactions can generate genetic diseases. It is thus important to control the aqueous interaction of proteins, by disrupting their aggregation or enhancing their specific interactions, respectively.

One mean of durable development involves transition to green chemistry. A way to reach this goal is to reduce the use of organic solvents and enhance the use of water-soluble molecules. One example would be the paint industry, by employing water-soluble polymers. To increase the time-stability of the paint it is necessary to control the degree of solubility of the polymers in the initial stage, before their irreversible entanglement and/or cross-linking. Small compounds can be used as cross-linkers, increasing in this way the domain specificity of the paints.

## 2. Prepare the working environment

## 3. Dataset preparation

### 3.1. Individual compounds and datasets with experimental solubility values

* from literature
* from dedicated data bases

### 3.2. Update the dataset

> Assan: where will you store the dataset (S3 ??)

The resulting data set will have the duplicates removed and the solubility values unified using units of mol/l at the room temperature and neutral pH.


### 3.3. Build the molecular fingerprint as SMILES
>> please explain to ppeople what smiles are

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


### 3.4. Structure-based solubility calculation. Case study examples

## 4. Decide on the deep learning model
> Assan: with the litterature/bibliography u will get insights for the model to choose :)
## 5. Learn and test the model

## 6. Predict the aqueous solubility of novel compounds w/o experimental solubility values

## 7. Construct a web application

> Assan :This application will be a Flask application deployed in the cloud. It will take as input the SMILES of the new compound and will display the solubility.

