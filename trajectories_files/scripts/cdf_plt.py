import os
from glob import glob
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import  argparse


argParser = argparse.ArgumentParser()
argParser.add_argument("-wt", "--wt", help="wild type sequence")
argParser.add_argument("-mt", "--mt", help="variant sequence")

args = argParser.parse_args()
print("args=%s" % args)

print("args.wt=%s" % args.wt)
print("args.mt=%s" % args.mt)

cdf_wt=np.loadtxt('CDF_WT.dat')
cdf_mt=np.loadtxt('CDF_'+str(args.mt)+'.dat')
ind_wt=[]
ind_mt=[]
for i in range(0,len(cdf_wt)):
   ind_wt.append(i)
for i in range(0,len(cdf_mt)):
   ind_mt.append(i)
plt.plot(ind_wt,cdf_wt,'o-' ,linewidth=2,label='wt')
plt.plot(ind_mt,cdf_mt,'+-',linewidth=2,label='mt')
plt.xlabel('n')
plt.ylabel('time (ps)')
plt.savefig('cdf_'+str(args.mt)+'.png')
