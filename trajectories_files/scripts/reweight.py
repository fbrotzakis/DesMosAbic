import sys
import os
from glob import glob
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import  argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-wt", "--wt", help="wild type sequence")
argParser.add_argument("-mt", "--mt", help="variant sequence")
argParser.add_argument("-rfolder", "--rfolder", help="specific trajectory folder")

args = argParser.parse_args()
#print("args=%s" % args)

#print("args.wt=%s" % args.wt)
#print("args.mt=%s" % args.mt)


def first_true(iterable, default=False, pred=None):
    """Returns the first true value in the iterable.

    If no true value is found, returns *default*

    If *pred* is not None, returns the first item
    for which pred(item) is true.

    """
    # first_true([a,b,c], x) --> a or b or c or x
    # first_true([a,b], x, f) --> a if f(a) else b if f(b) else x
    return next(filter(pred, iterable), default)


restypes=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
#H11
wt_seq=[]
mt_seq=[]
wt_seq.append(args.wt)
mt_seq.append(args.mt)
#H11
#wt_seq=['RVTRS']
#print(wt_seq,mt_seq)
#H11-H4
#mt_seq=['HYVSY']
#H11-D4 hybrid
#mt_seq=['EYVSY']
def reweight_function(mt_sequence,wt_sequence,pssm):
   arr=np.loadtxt(pssm) 
   mt_seq_p=1
   wt_seq_p=1
   for i in range(0,len(mt_sequence)):
      #print(i, restypes.index(mt_sequence[i]),mt_sequence[i],arr[i,restypes.index(mt_sequence[i])])
      mt_seq_p*=arr[i,restypes.index(mt_sequence[i])]
   summ=sum(arr[1,:])
   #print(summ)
   for i in range(0,len(wt_sequence)):
      #print(i, restypes.index(wt_sequence[i]),wt_sequence[i],arr[i,restypes.index(wt_sequence[i])])
      wt_seq_p*=arr[i,restypes.index(wt_sequence[i])]
   return mt_seq_p,wt_seq_p,mt_seq_p/wt_seq_p

#COLVAR FILE
colvar_arr=np.loadtxt('../COLVAR')
colvar_arr_analysis=np.loadtxt('COLVAR_ANALYSIS')
CV=colvar_arr[:,3].tolist()
distance=colvar_arr_analysis[:,1].tolist()
torsion=colvar_arr_analysis[:,2].tolist()
stateB=1.8
if any(num > stateB for num in CV):
   bias=colvar_arr[:,8].tolist()
   time=colvar_arr[:,0].tolist()
   a=0
   a_rew=0
   kT=2.47
   w_ratio_list=[]
   w_list=[]
   w_rew_list=[]
   #print(indexx)
   indexx=CV.index(first_true(CV, None, lambda x: x>stateB)) 
   #print('index',indexx)
   #Calc boos factors
   for j in range(0,indexx):
      #print(reweight_function(mt_seq[0],wt_seq[0],'pssm_'+str(j)+'.txt'))
      w_ratio=reweight_function(mt_seq[0],wt_seq[0],'pssm_'+str(j)+'.txt')[2]
      w_ratio_list.append(w_ratio)
      a+=np.exp(bias[j]/kT)
      a_rew+=np.exp(bias[j]/kT)*w_ratio
      w_list.append(np.exp(bias[j]/kT))
      w_rew_list.append(np.exp(bias[j]/kT)*w_ratio)

   print('boost',a,'boost rew',a_rew,"\n")
   print('t_esc',time[indexx]*a/indexx,'t_esc rew',time[indexx]*a_rew/indexx,"\n")
   print('t_rew/t_esc',(time[indexx]*a_rew/indexx)/(time[indexx]*a/indexx),"\n")
   print('last t',time[indexx])
   jumptimes_wt = open('../../JUMPTIMES_WT',mode='a')
   jumptimes_wt.write(str(time[indexx]*a/indexx)+"\n")
   jumptimes_wt.close()
   jumptimes_mt = open('../../JUMPTIMES_'+str(mt_seq[0]),mode='a')
   jumptimes_mt.write(str(time[indexx]*a_rew/indexx)+"\n")
   jumptimes_mt.close()

   weights_f = open('../../weights_'+str(mt_seq[0])+'.dat',mode='a')
   rmsd=np.loadtxt('../COLVAR')[:,3]

   for i in range(0,len(w_list)):
      weights_f.write(str(w_list[i])+" "+str(w_rew_list[i])+" "+str(rmsd[i])+" "+str(distance[i])+" "+str(torsion[i])+" "+str(i)+" "+str(args.rfolder)+"\n")
      #weights_f.write(str(w_list[i])+" "+str(w_rew_list[i])+" "+str(rmsd[i])+" "+str(distance[i])+" "+str(torsion[i])+"\n")
else:
   sys.exit("No state B")


