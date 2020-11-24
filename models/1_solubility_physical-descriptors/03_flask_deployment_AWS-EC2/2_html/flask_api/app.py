import numpy as np
from flask import Flask, request, jsonify, render_template
import pickle

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

#from sklearn.model_selection import train_test_split
#from sklearn import linear_model
#from sklearn.metrics import mean_squared_error, r2_score

app = Flask(__name__)

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

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict',methods=['POST'])
def predict():

    input_features = [str(x) for x in request.form.values()]
    final_features = input_features[0]

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
    
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8080)

#if __name__ == "__main__":
#    app.run(debug=True)

