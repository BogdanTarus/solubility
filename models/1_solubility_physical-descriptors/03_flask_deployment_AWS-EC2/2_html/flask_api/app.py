# run it with:
# python3 app.py

#import the necessary libraries
import numpy as np
from flask import Flask, request, jsonify, render_template
import pickle

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

# load the model from disk
loaded_model = pickle.load(open('model_desc_total.pkl', 'rb'))

def AromaticAtoms(m):
    """
    This function calculates the number of aromatic atoms in a molecule.
    The input argument, m, is a rdkit molecular object
    """
    
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
    aa_count = []
    for i in aromatic_atoms:
        if i==True:
            aa_count.append(1)
        
    sum_aa_count = sum(aa_count)
    
    return sum_aa_count


def predictSingle(smiles, model):
    """
    This function predicts the four molecular descriptors: the octanol/water partition coefficient (LogP),
    the molecular weight (Mw), the number of rotatable bonds (NRb), and the aromatic proportion (AP) 
    for a single molecule
    
    The input arguments are SMILES molecular structure and the trained model, respectively.
    """
    
    # define the rdkit moleculat object
    mol = Chem.MolFromSmiles(smiles)
    
    # calculate the log octanol/water partition descriptor
    single_MolLogP = Descriptors.MolLogP(mol)
    
    # calculate the molecular weight descriptor
    single_MolWt   = Descriptors.MolWt(mol)
    
    # calculate of the number of rotatable bonds descriptor
    single_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
    
    # calculate the aromatic proportion descriptor
    single_AP = AromaticAtoms(mol)/Descriptors.HeavyAtomCount(mol)
    
    # put the descriptors in a list
    single_list = [single_MolLogP, single_MolWt, single_NumRotatableBonds, single_AP]
    
    # add the list to a pandas dataframe
    single_df = pd.DataFrame(single_list).T
    
    # rename the header columns of the dataframe
    single_df.columns = ['MolLogP', 'MolWt', 'NumRotatableBonds', 'AromaticProportion']
    #return single_df
    
    return model.predict(single_df)[0]

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict',methods=['POST'])
def predict():

    input_features = [str(x) for x in request.form.values()]
    final_features = input_features[0]

    # check the syntax of the SMILES structure
    try:
      check = Chem.MolFromSmiles(final_features)
      if check is None:        
        return render_template('index.html', prediction_text='Invalid SMILES. Try again')
      else:
        predOUT = predictSingle(final_features, loaded_model)    
        output = str(round(predOUT,5))

        return render_template('index.html', prediction_text='LogS should be {} (mol/L)'.format(output))
    except:
      pass

# comment the next two rows to run the flask locally
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8080)

# uncomment the next two lines to run the flask locally
#if __name__ == '__main__':
#    app.run(debug=True)
