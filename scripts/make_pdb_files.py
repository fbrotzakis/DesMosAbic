
import mdtraj as md
import os
from glob import glob 
import subprocess


chainA_l=1555
chainB_l=988

os.system('echo "0"|gmx_mpi trjconv -f run.xtc  -s run.tpr -sep -o .pdb')

pdb_files = glob('./*.pdb')
for f in range(0,len(pdb_files)):
   pdb=str(f)+'.pdb'
   topology = md.load(pdb).topology
   table, bonds = topology.to_dataframe()
   traj=md.load(pdb, top=pdb)
   print(table)
   #Look at the gro file and check how many residues does the chain have. 

   chain_names=[]
   for i in range(0,chainA_l):
      chain_names.append('A')
   for i in range(0,chainB_l):
      chain_names.append('B')    
   print(chain_names)

   for i in range(0,len(table)):
      table.loc[table.index[i],'chainID']=chain_names[i]

   print(table,'n')
   topology_ren = md.Topology.from_dataframe(table, bonds)

   print(topology_ren)
   traj = md.load_pdb(pdb,top=topology_ren)
   traj.save('chain'+str(f)+'.pdb')
   subprocess.call(['python', 'design.py', 'chain'+str(f)+'.pdb',str(f)])
