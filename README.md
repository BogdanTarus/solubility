<style type="text/css">
  h2 { margin-left: 10px; }
  h3 { margin-left: 20px; }
</style>

>Assan:
Bogdan pls try to to organize it with 1. 2. etc.
the style section did not work on my end, just a FYI
- pls add a description on the tab on the right corner thx.
- pls add pictures ( to explain SMiles)
- pls add as well resources for ppl wishing to discover what smiles are (links)
>


# Prediction of aqueous solubility of novel drug-like compounds
> Assan:  Why is it relevant to you , to the community etc.

## Overview

## Prepare the working environment

## Dataset preparation

### Individual compounds and datasets with experimental solubility values

* from literature
* from dedicated data bases

### Update the dataset

> Assan: where will you store the dataset (S3 ??)


### Build the molecular fingerprint as SMILES
>> please explain to ppeople what smiles are

The resulting data set will have the duplicates removed and the solubility values unified using units of mol/l at the room temperature and neutral pH.

## Decide on the deep learning model
> Assan: with the litterature/bibliography u will get insights for the model to choose :)
## Learn and test the model

## Predict the aqueous solubility of novel compounds w/o experimental solubility values

## Construct a web application

> Assan :This application will be a Flask application deployed in the cloud. It will take as input the SMILES of the new compound and will display the solubility.

