#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 15:19:59 2021

@author: gabriela
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os.path
import copy
import statsmodels.stats.multitest as multi

def ttest(cind1,cind_multi):  # computes t and p values from a t-test
	tval, pval = stats.ttest_rel(cind_multi,cind1)
	return tval, pval

# heatmap over a numpy matrix cavm1
def heatm(cavm1,cancers,*args):  # title & save-as name & p-vals optional args
	# optional arguments:
	# arg1: title to the figure
	# arg2: save-as name, if not provided, figure doesn't get saved
	# arg3: p-vals associated to cavm1 values, to know where to put asterisk
	
	np2 = cavm1.shape[1]
	fig, ax = plt.subplots(figsize=(25,20))
	plt.rcParams.update({'font.size': 16})
	plt.imshow(cavm1)
	
	if len(args)>=1: # there is a title provided
		plt.title(args[0],fontsize = 40)
	if len(args)>=3:  # optionally p-vals provided
		pvals = copy.deepcopy(args[2])
	plt.rcParams.update({'font.size': 11})
# Loop over data dimensions and create text annotations.
	for i in range(np2):
		if len(args)>=3:  # multiple test correction
			pvals[i,np.isfinite(pvals[i,:])] = multi.fdrcorrection(pvals[i,np.isfinite(pvals[i,:])])[1]
		for j in range(np2):
			if (len(args)>=3 and (np.abs(cavm1[i,j])<0.01) and (pvals[i,j]<0.05)):  #((pvals[i,j]*(len(cancers)))<0.05) ):
				text = ax.text(j, i,'<0.01*', ha="center", va="center", color="w")
#				text = ax.text(j, i,'<*', ha="center", va="center", color="w")				
			elif (len(args)>=3 and (np.abs(cavm1[i,j])<0.01)):
				text = ax.text(j, i, '<0.01', ha="center", va="center", color="w")
#				text = ax.text(j, i, '<', ha="center", va="center", color="w")				
			elif (len(args)>=3 and (pvals[i,j])<0.05): #(pvals[i,j]*(len(cancers)))<0.05):  # tval significant -> add asterisk

				text = ax.text(j, i, str(np.around(cavm1[i, j],2))+'*', ha="center", va="center", color="w")
			else:

				text = ax.text(j, i, np.around(cavm1[i, j],2), ha="center", va="center", color="w")

# We want to show all ticks...
	ax.set_xticks(np.arange(np2))
	ax.set_yticks(np.arange(np2))
# ... and label them with the respective list entries
	canup = [str.upper(cancers[ii]) for ii in range(len(cancers))]  # upper case cancer names
#	ax.set_xticklabels(cancers, fontsize = 16)
#	ax.set_yticklabels(cancers, fontsize = 16)
	ax.set_xticklabels(canup, fontsize = 30)
	ax.set_yticklabels(canup, fontsize = 30)

	plt.setp(ax.get_xticklabels(), rotation=90)
	
	plt.rcParams.update({'font.size': 20})
	plt.colorbar()

	if len(args) >= 2:  # if save-as name is provided
		name = args[1]
		plt.savefig('heatmap'+str(name)+'.pdf',bbox_inches='tight')
		
	#plt.show()
#####################################################################################
def congr(a,b):  # Tucker congruence coefficient
	rc = np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
	return rc
###############################################################################
def cong_av(bigvec):  # Tucker congruence coefficient: averaged
	nc,nrv,npars = bigvec.shape
	rcvec = np.zeros(nc)
	cong_pairs = [[] for kk in range(nc)]
	for cc in range(nc):
		rc = 0  # initialize
		for ii in range(nrv):
			for jj in range(ii+1,nrv): # all combinations only once
				temp = congr(bigvec[cc,ii,:],bigvec[cc,jj,:])
				rc += temp
				if ((ii%2==0) and (jj==ii+1)):
					cong_pairs[cc] = np.append(cong_pairs[cc],temp)
		rcvec[cc] = 2/nrv/(nrv-1)*rc
	return rcvec, cong_pairs
###############################################################################
## start here #################################################################
###############################################################################
# all unique cancers
cancers = ['acc','blca','brca','cesc','chol','coadread','dlbc','esca','gbm','hnsc',
			  'kich','kirc','kirp','laml','lgg','lihc','luad','lusc','meso',
			  'ov','paad','prad','sarc','skcm','stad','thca','thym',
		   'ucec','ucs','uvm']
NC = len(cancers) # how many unique elements	
reg = '22'   # which coupling term
tval_mat, pval_mat = np.zeros((NC,NC)), np.zeros((NC,NC))
cdiff_av = np.zeros((NC,NC)) # average difference in cvals
cong_single_mat, cong_multi_mat = np.zeros((NC,NC)), np.zeros((NC,NC))
pval_cong = np.zeros((NC,NC))
csing_av , cmulti_av = np.zeros((NC,NC)), np.zeros((NC,NC))
for jj in range(NC):  # all pairs
	for kk in range(jj+1,NC):
		#cancers_type = ['ov','paad']  # cancer tuples
		cancers_type = [cancers[jj],cancers[kk]]
#		tit = maketitle(cancers_type)
		nc = len(cancers_type)

		labels_t = np.zeros(nc)  # t-vals for each bar
		labels_p = np.zeros(nc)  # p-vals
		lab = [[] for cc in range(nc)]

		aas, aas2 = '', ''
		for cc in range(nc):
			aas = aas + cancers_type[cc];
			aas = aas + '_';
		for cc in range(nc-1,-1,-1):  # the data file could have a reversed name
			aas2 = aas2 + cancers_type[cc];
			aas2 = aas2 + '_';

		if (os.path.isfile('./results/cval_single_'+aas+'reg'+str(reg)+'.npy') and os.path.isfile('./results/beta_single_'+aas+'reg'+str(reg)+'.npy')): # is there a file of this name?
			cind1 = np.load('./results/cval_single_'+aas+'reg'+str(reg)+'.npy') 
			cind2 = np.load('./results/cval_multi_'+aas+'reg'+str(reg)+'.npy') 
		elif (os.path.isfile('./results/cval_single_'+aas2+'reg'+str(reg)+'.npy') and (os.path.isfile('./results/beta_single_'+aas2+'reg'+str(reg)+'.npy'))): # arguments could be reversed
			cind1 = np.load('./results/cval_single_'+aas2+'reg'+str(reg)+'.npy') 
			cind2 = np.load('./results/cval_multi_'+aas2+'reg'+str(reg)+'.npy') 
			# reverse cvals
			temp1 = copy.deepcopy(cind1[0,:])
			cind1[0,:] = copy.deepcopy(cind1[1,:])
			cind1[1,:] = copy.deepcopy(temp1)
			temp2 = copy.deepcopy(cind2[0,:])
			cind2[0,:] = copy.deepcopy(cind2[1,:])
			cind2[1,:] = copy.deepcopy(temp2)
		else: # otherwise warning
#			raise ValueError('Cancer combination '+ aas + ' does not exist.')
			tval_mat[jj,kk], tval_mat[kk,jj] = float('NaN'), float('NaN')
			pval_mat[jj,kk], pval_mat[kk,jj] = float('NaN'), float('NaN')
#			print('No c-value data for ' + aas + ' available.')
			raise ValueError('No c-value data for ' + aas + ' available.')
#			break
			continue
		
		# check if there are any 0 and remove them
		indd = ((cind2[0,:]>1e-6) & (cind1[0,:]>1e-6) & (cind2[1,:]>1e-6) & (cind1[1,:]>1e-6))
		if (np.sum(indd)!=len(indd)): # not all cvals are above 0
			cind2 = cind2[:,indd]
			cind1 = cind1[:,indd]
#			can_error.append(can)
			print(str(np.sum(~indd))+' cvals removed from '+ aas)  # print warning

		tt, pp = np.zeros(nc), np.zeros(nc)  # tempoary holder for tvals and pvals
		cdav = np.zeros(nc)  # average cval difference
		csing = np.zeros(nc)
		cmulti = np.zeros(nc)
		for cc in range(nc):
			cind_single = cind1[cc,:]
			cind_multi = cind2[cc,:]

			tt[cc], pp[cc] = ttest(cind_single,cind_multi)
			cdav[cc] = np.mean(cind_multi-cind_single) # avergae c_multi- c_single
			csing[cc] = np.mean(cind_single)
			cmulti[cc] = np.mean(cind_multi)
		# save t and p values
		tval_mat[jj,kk], tval_mat[kk,jj] = tt[0], tt[1]
		pval_mat[jj,kk], pval_mat[kk,jj] = pp[0], pp[1]
		cdiff_av[jj,kk], cdiff_av[kk,jj] = cdav[0], cdav[1]
		csing_av[jj,kk], csing_av[kk,jj] = csing[0],csing[1]
		cmulti_av[jj, kk], cmulti_av[kk, jj] = cmulti[0], cmulti[1]
		# compute average congruence
		if (os.path.isfile('./results/cval_single_'+aas+'reg'+str(reg)+'.npy') and os.path.isfile('./results/beta_single_'+aas+'reg'+str(reg)+'.npy')): # is there a file of this name?
			beta1 = np.load('./results/beta_single_'+aas+'reg'+str(reg)+'.npy') 
			beta2 = np.load('./results/beta_multi_'+aas+'reg'+str(reg)+'.npy') 
			# remove failed computations
			if (np.sum(indd)!=len(indd)): # not all cvals are above 0
				beta1 = beta1[:,indd,:]
				beta2 = beta2[:,indd,:]
#			cong_single_mat[jj,kk], cong_single_mat[kk,jj] = cong_av(beta1)
#			cong_multi_mat[jj,kk], cong_multi_mat[kk,jj] = cong_av(beta2)
			r11, r12 = cong_av(beta1)
			cong_single_mat[jj,kk], cong_single_mat[kk,jj] = r11[0], r11[1]
			r21, r22 = cong_av(beta2)
			cong_multi_mat[jj,kk], cong_multi_mat[kk,jj] = r21[0], r21[1]
						
			pval_cong[jj,kk] = ttest(r12[0],r22[0])[1]
			pval_cong[kk,jj] = ttest(r12[1],r22[1])[1]
		elif (os.path.isfile('./results/cval_single_'+aas2+'reg'+str(reg)+'.npy') and (os.path.isfile('./results/beta_single_'+aas2+'reg'+str(reg)+'.npy'))): # arguments could be reversed
#			print('reversed')
			beta1 = np.load('./results/beta_single_'+aas2+'reg'+str(reg)+'.npy') 
			temp = copy.copy(beta1[0,:,:])
			beta1[0,:,:] = beta1[1,:,:] # reverse variables
			beta1[1,:,:] = temp
			beta2 = np.load('./results/beta_multi_'+aas2+'reg'+str(reg)+'.npy') 
			temp2 = copy.copy(beta2[0,:,:])
			beta2[0,:,:] = beta2[1,:,:] # reverse variables
			beta2[1,:,:] = temp2
			# remove failed computations
			if (np.sum(indd)!=len(indd)): # not all cvals are above 0
				beta1 = beta1[:,indd,:]
				beta2 = beta2[:,indd,:]


			r11, r12 = cong_av(beta1)
			cong_single_mat[jj,kk], cong_single_mat[kk,jj] = r11[0], r11[1]
			r21, r22 = cong_av(beta2)
			cong_multi_mat[jj,kk], cong_multi_mat[kk,jj] = r21[0], r21[1]
						
			pval_cong[jj,kk] = ttest(r12[0],r22[0])[1]
			pval_cong[kk,jj] = ttest(r12[1],r22[1])[1]

		else: # otherwise warning
			cong_single_mat[jj,kk], cong_single_mat[kk,jj] = float('NaN'), float('NaN')
			cong_multi_mat[jj,kk], cong_multi_mat[kk,jj] = float('NaN'), float('NaN')
			print('No beta data for ' + aas + ' available.')
			
		# remove the diagonal elements
	tval_mat[jj,jj] = float('NaN')
	pval_mat[jj,jj] = float('NaN')
	cong_single_mat[jj,jj] = float('NaN')
	cong_multi_mat[jj,jj] = float('NaN')	
	cdiff_av[jj,jj] = float('NaN')	

# plot heatmaps
heatm(csing_av, cancers, 'Single cancer c-value', 'single_'+reg)
heatm(cmulti_av, cancers, 'Multi cancer c-value', 'multi_'+reg)
heatm(tval_mat, cancers, 't-value','tvals_'+reg, pval_mat)
heatm(cdiff_av, cancers, 'Average c-value difference', 'cdiff_av_'+reg, pval_mat)
heatm(pval_mat, cancers, 'p-value','pvals_'+reg,pval_mat)
heatm(cong_single_mat, cancers, 'Average singletask congruence','cong_single_'+reg)
heatm(cong_multi_mat, cancers, 'Average multitask congruence','cong_multi_'+reg)
heatm(cong_multi_mat-cong_single_mat, cancers, 'Average congruence difference','cong_diff_'+reg, pval_cong)
