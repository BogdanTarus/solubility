
# Prediction of the aqueous solubility directly from the chemical strcture

This work is based on the research of Lambard and Gracheva (Guillaume Lambard and Ekaterina Gracheva 2020 Mach. Learn.: Sci. Technol. 1 025004 
). It was slightly modified to run on Google Collab.

The major inconvenient of a model that uses physical descriptors to predict the solubility is that they are empirically engineered by humans. It is thus desirable to propose methods that avoid the subjectivity of the hand-encoded features.

The SMILESX workflow propose an algorithm that takes as input the SMILES structure of chemical compounds and predicts their solubility. As a SMILES specification is represented by a string, one can use techniques specific to the Natural Language Processing (NLP) to generate models for de novo drug design. 

Another problem in training models for drug design predictions is the size of the training dataset. In general, the number of the samples in the dataset is of the order of $10^3$, too small to apply modern neural network methods directly. The scarcity of the dataset is given by the experimental difficulty to get data associated to the descriptors. The SMILESX algorithm also propose a way to augment the size of the dataset by, in a first step, removing the canonicalization of the SMILES specifications. In a second step, the atoms of a given SMILES are renumbered by rotating their index correlated with a reconstruction of the correct SMILES syntax.

 However, using abstract features to construct a neural architecture makes it difficult to interpret their contribution to the observable of interest, aqueous solubility in our case. Adding an attention mechanism to the algorithm make it possible to both read deeper into the SMILES and interpret, at no extra cost, the output of the model.

There are two jupyter notebooks to be run, the first one for training the model (SMILESX_Prediction.ipynb) and the second one (SMILESX_Visualization.ipynb) to interpret the output. Both of them where slightly modified from the original ones to run on Google Collab. One should make sure to have in the working directory on Google drive the two notebooks, the SMILESX folder which contains the python scripts that train the model and activate the attention mechanism. Also, in the same working directory we need to have a folder with the datasets.
