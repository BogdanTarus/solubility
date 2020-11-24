
from flask import Flask, request
import pickle

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

# load the model from disk
loaded_model = pickle.load(open('model_desc_total.pkl', 'rb'))

def AromaticAtoms(m):
  aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic() for i in range(m.GetNumAtoms())]
  aa_count = []
  for i in aromatic_atoms:
    if i==True:
      aa_count.append(1)
  sum_aa_count = sum(aa_count)
  return sum_aa_count

def predictSingle(smiles, model):
    mol = Chem.MolFromSmiles(smiles)
    
    single_MolLogP = Descriptors.MolLogP(mol)
    single_MolWt   = Descriptors.MolWt(mol)
    single_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
    
    single_AP = AromaticAtoms(mol)/Descriptors.HeavyAtomCount(mol)
    
    single_list = [single_MolLogP, single_MolWt, single_NumRotatableBonds, single_AP]
    single_df = pd.DataFrame(single_list).T
    single_df.columns = ['MolLogP', 'MolWt', 'NumRotatableBonds', 'AromaticProportion']
    
    return model.predict(single_df)[0]

@app.route('/', methods=['POST'])
def home():
    smiles = request.json.get("smiles")

    try:
      check = Chem.MolFromSmiles(smiles)
      if check is None:        
        return 'Invalid SMILES. Try again' + "\n"
      else:
        predOUT = predictSingle(smiles, loaded_model)
        predOUT = round(predOUT, 5)

        return str('LogS should be: ' + str(predOUT) + ' (mol/L)') + "\n"
    except:
      pass

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8080)
 
#if __name__ == '__main__':
#    app.run(debug=True)


