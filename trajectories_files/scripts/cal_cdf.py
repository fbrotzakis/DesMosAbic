#iMetaD KS test
#use python3 <n> 
#t = vector with rescaled times for each event (iMetaD times)
#ko = initial guess, here the average rescaled time
#n = total number of simulations
#input_t = deposition time or PACE of the gaussians [ps]

import scipy
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import ks_2samp
from sys import argv
import matplotlib.pyplot as plt
import  argparse


argParser = argparse.ArgumentParser()
argParser.add_argument("-wt", "--wt", help="wild type sequence")
argParser.add_argument("-mt", "--mt", help="variant sequence")

args = argParser.parse_args()
print("args=%s" % args)

print("args.wt=%s" % args.wt)
print("args.mt=%s" % args.mt)



#n = argv[1]
t=np.loadtxt('JUMPTIMES_WT') #filename of the file with rescaled times
n=len(t)
print(n)
mintime=min(t)
maxtime=max(t)

ko=1/np.average(t)
l=int(len(t))

#Use logscale to facilitate the fit
logbins = np.logspace(np.log10(mintime),np.log10(maxtime),l)
#bins = [mintime,maxtime,l]
#print(bins)

#Calculate ECDF
counts, bin_edges = np.histogram(t, bins=logbins)
ecdf = np.cumsum(counts)

#print(logbins)
#print(ecdf)

#Fit ECDF to TCDF of Poisson process
k,e_k = scipy.optimize.curve_fit(lambda t,b: 1-np.exp(-b*t), logbins[0:-1],  ecdf/float(n), ko)

#Estimate MFPT from the fit
tau = 1/k

theoretical_times=np.random.gamma(1,tau,(1,10**(6)))
fil_theo_t = theoretical_times[theoretical_times < maxtime*0.95 ]

#KS test
s,p_val = scipy.stats.ks_2samp(t,fil_theo_t)

statistics = open('statistics_WT',mode='w')
statistics.write("tau (ps)"+str(tau)+" koff "+str(1/tau)+" pval "+str(p_val)+"\n")
statistics.write("tau (s) "+str(tau*pow(10,-12))+" koff "+str(1/(tau*pow(10,-12)))+" pval "+str(p_val))
statistics.close()
############################### MT #################################
t=np.loadtxt('JUMPTIMES_'+str(args.mt)) #filename of the file with rescaled times
n=len(t)
print(n)
mintime=min(t)
maxtime=max(t)

ko=1/np.average(t)
l=int(len(t))

#Use logscale to facilitate the fit
logbins = np.logspace(np.log10(mintime),np.log10(maxtime),l)
#bins = [mintime,maxtime,l]
#print(bins)

#Calculate ECDF
counts, bin_edges = np.histogram(t, bins=logbins)
ecdf = np.cumsum(counts)

#print(logbins)
#print(ecdf)

#Fit ECDF to TCDF of Poisson process
k,e_k = scipy.optimize.curve_fit(lambda t,b: 1-np.exp(-b*t), logbins[0:-1],  ecdf/float(n), ko)

#Estimate MFPT from the fit
tau = 1/k

theoretical_times=np.random.gamma(1,tau,(1,10**(6)))
fil_theo_t = theoretical_times[theoretical_times < maxtime*0.95 ]

#KS test
s,p_val = scipy.stats.ks_2samp(t,fil_theo_t)
statistics = open('statistics_'+str(args.mt),mode='w')

statistics.write("tau (ps)  "+str(tau)+" koff "+str(1/tau)+" pval "+str(p_val)+"\n")
statistics.write("tau (s) "+str(tau*pow(10,-12))+" koff "+str(1/(tau*pow(10,-12)))+" pval "+str(p_val))
t_esc_file=open('t_esc_file',mode='w')
t_esc_file.write(str(tau[0]*pow(10,-12)))
pval_file=open('pval_file',mode='w')
pval_file.write(str(p_val))

t_esc_file.close()
pval_file.close()
statistics.close()


