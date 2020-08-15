<style type="text/css">
  h2 { margin-left: 10px; }
  h3 { margin-left: 20px; }
</style>

# Prediction of aqueous solubility of novel drug-like compounds

## Overview

## Prepare the working environment

## Dataset preparation

### Individual compounds and datasets with experimental solubility values

* from literature
* from dedicated data bases

### Update the dataset

### Build the molecular fingerprint as SMILES

The resulting data set will have the duplicates removed and the solubility values unified using units of mol/l at the room temperature and neutral pH.

## Decide on the deep learning model

## Learn and test the model

## Predict the aqueous solubility of novel compounds w/o experimental solubility values

## Construct a web application

This application will be a Flask application deployed in the cloud. It will take as input the SMILES of the new compound and will display the solubility.

