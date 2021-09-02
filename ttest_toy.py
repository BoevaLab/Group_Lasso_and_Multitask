#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 13:00:15 2021

@author: gabriela
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib as mtplt
import pandas as pd
import copy
import statsmodels.stats.multitest as multi;

def ttest_old(cind1,cind_multi):  # computes t and p values from a t-test
#	corr = np.corrcoef(cind_multi.flatten(),cind1.flatten())
	tval, pval = stats.ttest_rel(cind_multi,cind1)
#	print(corr[0,1], pval)
	return tval, pval

def ttest(cind1,cind_multi,name,col1,col2,cc,*args):
	## cind1
	ll = np.arange(len(cind1))
	corr = np.corrcoef(cind_multi.flatten(),cind1.flatten())
	tval, pval = stats.ttest_rel(cind_multi,cind1)
	print(corr[0,1], pval)
	
	# plot the c-values
	f = plt.figure(1)
#	if cc == 0:  # sort only once
#		indd = np.argsort(cind1)
	if len(args)==1: # if sorting provided
		indd = args[0]
		plt.plot(ll,cind1[indd],'-o',color = col1)
		plt.plot(ll,cind_multi[indd],'*-',color = col2)
	else:
		plt.plot(ll,cind1,'-o',color = col1)
		plt.plot(ll,cind_multi,'*-',color = col2)
#	plt.title('Toy example: c-values')
	plt.title('c-values')
	plt.legend(tit, loc='lower right',fontsize=16)
	plt.xlabel('random seed')
	plt.ylabel('c-value')
	plt.savefig('cor_'+name+'_reg'+str(reg)+'.pdf',bbox_inches='tight')
#	plt.show()
	f.show()
	
# BOXPLOT				
	g = plt.figure(2)
	bplot = plt.boxplot(cind_multi-cind1, positions = [cc], labels = [cancers_type[cc]], patch_artist=True,\
					boxprops=dict(facecolor=col1, color=col2), widths = 0.5)#,	color = ((ng-jjj)/ng,(jjj+1)/2/ng, 0)) # color grade)
#	plt.title('Differences in c-values')
	plt.title('Paired c-value differences')
#	plt.ylim((-0.02,0.027))
	plt.savefig('boxplot_'+name+'_reg'+str(reg)+'.pdf',bbox_inches='tight')
	g.show()
	
	return tval, pval

def maketitle(cancer_type):
	tit = []
	nc = len(cancer_type)
	for ii in range(nc):
#		tit.append(cancer_type[ii]+'_S')
#		tit.append(cancer_type[ii]+'_M')
		tit.append(cancer_type[ii]+'_single')
		tit.append(cancer_type[ii]+'_multi')
	return tit

#####################################################################################
def congr(a,b):  # Tucker congruence coefficient
	rc = np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
	return rc
#####################################################################################
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
def group_list(group):	# map groups into lists
	unq = np.unique(group)	# unique group names
	neunq = len(unq)		# n of unique elements

	D = dict() #  mapping the group names to their index
	for (i,x) in enumerate(unq):
		D[x] = i;

	x = [[] for i in range(neunq)]	# empty list of lists
	for jj in range(len(group)):
		x[D[group[jj]]].append(jj)	# append indices of gloup elements into respective lists
	return x, unq
###################################################################################
############### MAIN ##############################################################
###################################################################################
	
col_multi = ['deepskyblue','orange','darkolivegreen','darkviolet','firebrick','teal']
col_single = ['lightblue','gold','olivedrab','mediumorchid','indianred','cadetblue']
#
# larger font everywhere
plt.rcParams.update({'font.size': 20})
	
cancers_type = ['T1','T2']
tit = maketitle(cancers_type)
reg = '22'  # coupling term number
nc = len(cancers_type)
tval, pval = np.zeros(nc), np.zeros(nc)  # t and p values initialize

aas = ''
for cc in range(nc):
	aas = aas + cancers_type[cc];
	aas = aas + '_';

cind1 = np.load('./cval_single_'+aas+'reg'+str(reg)+'.npy')
cind2 = np.load('./cval_multi_'+aas+'reg'+str(reg)+'.npy') 
#indd = np.argsort(cind1[1,:]) # sort according to first cancer single

for cc in range(nc):
	name = aas
	cind_single = cind1[cc,:]
	cind_multi = cind2[cc,:]
#	tt, pp = ttest(cind_single,cind_multi,name,col_single[cc],col_multi[cc],cc)	
	tt, pp = ttest(cind_single,cind_multi,name,col_single[cc],col_multi[cc],cc,np.argsort(cind_single))
	
	print(cancers_type[cc]+': t-value = '+str(np.around(tt,2))+' , p-value = '+str(np.around(pp,2)))
	tval[cc] = tt
	pval[cc] = pp

# load betas	
beta1 = np.load('./beta_single_'+aas+'reg'+str(reg)+'.npy')
beta2 = np.load('./beta_multi_'+aas+'reg'+str(reg)+'.npy') 

# compute congruence
cong_single, cong_pairs_single = cong_av(beta1)
cong_multi, cong_pairs_multi = cong_av(beta2)

tval1, pval1 = ttest_old(cong_pairs_single[0],cong_pairs_multi[0])
tval2, pval2 = ttest_old(cong_pairs_single[1],cong_pairs_multi[1])

gr = pd.read_csv(cancers_type[0]+'_groups.csv').values[:,1:].flatten()
ng = np.max(gr)

# move the three important pathways 4, 54, 62 -> 1, 2, 3
gr_new = copy.deepcopy(gr)
gr_new[gr==1] = 4
gr_new[gr==2] = 62
gr_new[gr==3] = 54
gr_new[gr==4] = 1
gr_new[gr==62] = 2
gr_new[gr==54] = 3

ind = 5 # which index do you want to plot (0-50)
glist_orig, gunq_orig = group_list(gr_new)
cmap = mtplt.cm.get_cmap('PiYG') # viridis/Spectral/plasma
#rgba = cmap(0.5)
plt.figure() # new figure

for cc in range(nc):

	##### BOX PLOT	: VERSION 2
	fig = plt.figure()
	ax = fig.add_subplot(111)
	beg = 0   # plot begin pointer
	for jj in range(len(glist_orig)):
		if (np.linalg.norm(beta1[cc,ind,glist_orig[jj]]) > 1e-5):
			c = cmap(jj/ng)
			bplot = plt.boxplot(beta1[cc,ind,glist_orig[jj]], positions = [beg], labels = [jj+1], patch_artist=True,\
					boxprops=dict(facecolor=c, color=c), widths = 0.7)#,	color = ((ng-jjj)/ng,(jjj+1)/2/ng, 0)) # color grade)
			beg+=1;
	plt.setp(ax.get_xticklabels(), fontsize = 16)#, rotation=90)
	plt.title('Single tasking: '+ cancers_type[cc])
#	plt.xlabel('pathway')
#	if cc == 0: # print y label only for first
	plt.ylabel(r'$\beta$',fontsize=20)
	plt.xlabel('Pathways with nonzero '+ r'$\beta$',fontsize=20)
	plt.savefig('single'+str(cancers_type[0])+'_'+str(cancers_type[1])+'_'+str(cc)+'_boxplot.pdf',bbox_inches='tight')
	plt.show()

	fig = plt.figure()
	ax = fig.add_subplot(111)
	beg = 0   # plot begin pointer
	for jj in range(len(glist_orig)):
		if (np.linalg.norm(beta2[cc,ind,glist_orig[jj]]) > 1e-5):
			c = cmap(jj/ng)
			bplot = plt.boxplot(beta2[cc,ind,glist_orig[jj]], positions = [beg], labels = [jj+1], patch_artist=True,\
					boxprops=dict(facecolor=c, color=c), widths = 0.7)#,	color = ((ng-jjj)/ng,(jjj+1)/2/ng, 0)) # color grade)
			beg+=1;
	plt.setp(ax.get_xticklabels(), fontsize = 16)#, rotation=90)
	plt.title('Multi tasking: '+ cancers_type[cc])
	plt.xlabel('pathway')
#	if cc == 0: # print y label only for first
#		plt.ylabel(r'$\beta$')
	plt.ylabel(r'$\beta$',fontsize=20)
	plt.xlabel('Pathways with nonzero '+ r'$\beta$',fontsize=20)
	plt.savefig('multi_'+str(cancers_type[0])+'_'+str(cancers_type[1])+'_'+str(cc)+'_boxplot.pdf',bbox_inches='tight')
	plt.show()
