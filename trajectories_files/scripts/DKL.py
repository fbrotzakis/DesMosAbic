import sys
import os
from glob import glob
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from scipy.special import rel_entr
import  argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-wt", "--wt", help="wild type sequence")
argParser.add_argument("-mt", "--mt", help="variant sequence")
argParser.add_argument("-path", "--path", help="path of trajectory  master folder")

args = argParser.parse_args()
#print("args=%s" % args)

#print("args.wt=%s" % args.wt)
#print("args.mt=%s" % args.mt)



weights=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,0]
weights_rew=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,1]
rmsd=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,2]

weights = np.asarray(weights, dtype=float)
weights_rew = np.asarray(weights_rew, dtype=float)
rmsd = np.asarray(rmsd, dtype=float)

weights /= weights.sum()
weights_rew /= weights_rew.sum()
DKL=sum(rel_entr(weights_rew,weights))
dkl_f=open(str(args.path)+'/DKL_'+str(args.mt)+'.dat','w')
dkl_f.write(str(DKL))
dkl_f.close()

weights=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,0]
weights_rew=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,1]
rmsd=np.loadtxt(str(args.path)+'/'+'weights_'+str(args.mt)+'.dat',dtype='str')[:,2]

weights = np.asarray(weights, dtype=float)
weights_rew = np.asarray(weights_rew, dtype=float)
rmsd = np.asarray(rmsd, dtype=float)


value1, bins1=np.histogram(rmsd, bins=100,  weights=weights, density=False)
value2, bins2=np.histogram(rmsd, bins=100,  weights=weights_rew, density=False)

#print(value2[-1],value1[-1])

psBA=[]
fsB=[]
value1_shift=[]
value2_shift=[]
for i in range(0,len(value1)):
   value1_shift.append(value1[i]-value1[-1]+1)
   value2_shift.append(value2[i]-value2[-1]+1)



plt.plot(bins1[:-1],-np.log(value1_shift),'o-' ,linewidth=2)
plt.plot(bins2[:-1],-np.log(value2_shift),'+-',linewidth=2)

plt.legend( ['wt', 'mt'],loc='upper right')
#for i in range(0,len(value1)):
#   print(bins1[i],value1[i],value2[i])
#print(value1[-2],value2[-2])
plt.xlabel('RMSD' )
plt.xlim(0,2.5 )
plt.ylabel('F (kT)' )
plt.savefig(str(args.path)+'/RMSD_FES_'+str(args.mt)+'.pdf',dpi=400,transparent=True, bbox_inches='tight')

plt.close()
plt.plot(bins1[:-1],value1,'o-' ,linewidth=2)
plt.plot(bins2[:-1],value2,'+-',linewidth=2)

plt.legend( ['wt', 'mt'],loc='upper right')
#for i in range(0,len(value1)):
#   print(bins1[i],value1[i],value2[i])
#print(value1[-2],value2[-2])
plt.xlabel('RMSD' )
plt.xlim(0,2.5 )
plt.ylabel('counts' )
plt.savefig('RMSD_prob_'+str(args.mt)+'.pdf',dpi=400,transparent=True, bbox_inches='tight')

rmsd_minimum=open(str(args.path)+'/rmsd_minimum_'+str(args.mt)+'.dat','w')
rmsd_minimum.write(str((bins2[np.where(value2 == value2.max())][0])))
rmsd_minimum.close()
print(bins2[np.where(value2 == value2.max())],value2.max())

