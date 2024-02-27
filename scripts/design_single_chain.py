#@title Install colabdesign
import os
import sys
import colabdesign

from colabdesign.mpnn import mk_mpnn_model, clear_mem
from colabdesign.shared.protein import pdb_to_string

import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import HTML
import pandas as pd
import tqdm.notebook
TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'




#from google.colab import files
#from google.colab import data_table
#data_table.enable_dataframe_formatter()

def get_pdb(pdb_code=""):
  if pdb_code is None or pdb_code == "":
    upload_dict = files.upload()
    pdb_string = upload_dict[list(upload_dict.keys())[0]]
    with open("tmp.pdb","wb") as out: out.write(pdb_string)
    return "tmp.pdb"
  elif os.path.isfile(pdb_code):
    return pdb_code
  elif len(pdb_code) == 4:
    os.system(f"wget -qnc https://files.rcsb.org/view/{pdb_code}.pdb")
    return f"{pdb_code}.pdb"
  else:
    os.system(f"wget -qnc https://alphafold.ebi.ac.uk/files/AF-{pdb_code}-F1-model_v3.pdb")
    return f"AF-{pdb_code}-F1-model_v3.pdb"

import warnings,  re
warnings.simplefilter(action='ignore', category=FutureWarning)

os.system("mkdir -p output")

# USER OPTIONS
#@markdown #### ProteinMPNN options
model_name = "v_48_020" #@param ["v_48_002", "v_48_010", "v_48_020", "v_48_030"]
#@markdown #### Input Options
pdb='6MRR' #@param {type:"string"}
#@markdown - leave blank to get an upload prompt
chains = "A,B" #@param {type:"string"}
homooligomer = False #@param {type:"boolean"}
#@markdown #### Design constraints
fix_pos = "B296-300" #@param {type:"string"}  #Remember that B100 is B101 real residue
#@markdown - specify which positions to keep fixed in the sequence (example: `1,2-10`)
#@markdown - you can also specify chain specific constraints (example: `A1-10,B1-20`)
#@markdown - you can also specify to fix entire chain(s) (example: `A`)
#inverse = False #@param {type:"boolean"}
inverse = True #@param {type:"boolean"}
#@markdown - inverse the `fix_pos` selection (define position to "free" [or design] instead of "fix")
rm_aa = "" #@param {type:"string"}
#@markdown - specify amino acid(s) to exclude (example: `C,A,T`)

#@markdown #### Design Options
num_seqs = 32 #@param ["32", "64", "128", "256", "512", "1024"] {type:"raw"}
sampling_temp = 0.1 #@param ["0.0001", "0.1", "0.15", "0.2", "0.25", "0.3", "0.5", "1.0"] {type:"raw"}
#@markdown - Sampling temperature for amino acids, T=0.0 means taking argmax, T>>1.0 means sample randomly.

#@markdown Note: designed sequences are saved to `design.fasta`

# cleaning user options
chains = re.sub("[^A-Za-z]+",",", chains)
if fix_pos == "": fix_pos = None
rm_aa = ",".join(list(re.sub("[^A-Z]+","",rm_aa.upper())))
if rm_aa == "": rm_aa = None

directory = os.getcwd()

#pdb_path = get_pdb(pdb)
print(sys.argv[1])
pdb_path = directory+'/'+sys.argv[1]
jj =sys.argv[2]
if "mpnn_model" not in dir():
  mpnn_model = mk_mpnn_model(model_name)

mpnn_model.prep_inputs(pdb_filename=pdb_path,
                       chain=chains, homooligomer=homooligomer,
                       fix_pos=fix_pos, inverse=inverse,
                       rm_aa=rm_aa, verbose=True)
out = mpnn_model.sample(num=num_seqs//32, batch=32,
                        temperature=sampling_temp,
                        rescore=homooligomer)

with open("design.fasta","w") as fasta:
  for n in range(num_seqs):
    line = f'>score:{out["score"][n]:.3f}_seqid:{out["seqid"][n]:.3f}\n{out["seq"][n]}'
    fasta.write(line+"\n")

labels = ["score","seqid","seq"]
data = [[out[k][n] for k in labels] for n in range(num_seqs)]

df = pd.DataFrame(data, columns=labels)
df.to_csv('output/mpnn_results.csv')
#data_table.DataTable(df.round(3))



#@title ### Get amino acid probabilties from ProteinMPNN (optional)
mode = "conditional" #@param ["unconditional", "conditional", "conditional_fix_pos"]
#@markdown - `unconditional` - P(sequence | structure) 
#@markdown - `conditional` - P(sequence | structure, sequence)
#@markdown - `conditional_fix_pos` - P(sequence[not_fixed] | structure, sequence[fix_pos])
show = "all" 
#import plotly.express as px
from scipy.special import softmax
from colabdesign.mpnn.model import residue_constants
L = sum(mpnn_model._lengths)
fix_pos = mpnn_model._inputs.get("fix_pos",[])
free_pos = np.delete(np.arange(L),fix_pos)

if mode == "conditional":
  ar_mask = 1-np.eye(L)
  logits = mpnn_model.score(ar_mask=ar_mask)["logits"]
  pdb_labels = None
if mode == "conditional_fix_pos":
  assert "fix_pos" in mpnn_model._inputs, "no positions fixed"
  ar_mask = 1-np.eye(L)
  p = np.delete(np.arange(L),mpnn_model._inputs["fix_pos"])
  ar_mask[free_pos[:,None],free_pos[None,:]] = 0
  logits = mpnn_model.score(ar_mask=ar_mask)["logits"]
  logits = logits[free_pos]
  pdb_labels = np.array([f"{i}_{c}" for c,i in zip(mpnn_model.pdb["idx"]["chain"], mpnn_model.pdb["idx"]["residue"])])
  pdb_labels = pdb_labels[free_pos]
else:
  ar_mask = np.zeros((L,L))
  logits = mpnn_model.score(ar_mask=ar_mask)["logits"]
  pdb_labels = None

pssm = softmax(logits,-1)
np.savetxt("output/pssm_"+str(jj)+".txt",pssm)

#fig = px.imshow(np.array(pssm).T,
#               labels=dict(x="positions", y="amino acids", color="probability"),
#               y=residue_constants.restypes + ["X"],
#               x=pdb_labels,
#               zmin=0,
#               zmax=1,
#               template="simple_white",
#              )
#fig.update_xaxes(side="top")
#fig.show()
