#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 15:19:59 2021

@author: gabriela
"""

from __future__ import division
import numpy as np
import matplotlib as mtplt
import matplotlib.pyplot as plt
from scipy import stats
import os.path
import pandas as pd

def ttest(cind1,cind_multi):  # computes t and p values from a t-test
	tval, pval = stats.ttest_rel(cind_multi,cind1)
	return tval, pval
######################################################################################
def congr(a,b):  # Tucker congruence coefficient
	rc = np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b)
	return rc
#####################################################################################
def cong_av(bigvec):  # Tucker congruence coefficient: averaged
	nc,nrv,npars = bigvec.shape
	if(nc>1):
		print('Something is wrong.')
	rcvec = np.zeros(nc)
	stdvec = np.zeros(nc)
#	rc_all = [[] for ii in range(nc)]  # all congruences
	for cc in range(nc):
		rc = 0  # initialize
		congtemp = []  # placeholder for cong coef to comp the 
		cong_pairs = []
		for ii in range(nrv):
			for jj in range(ii+1,nrv): # all combinations only once
				temp = congr(bigvec[cc,ii,:],bigvec[cc,jj,:])
				rc += temp
				congtemp = np.append(congtemp,temp)
				if ((ii % 2 == 0) and (jj == ii + 1)): # i even and j one larger
					cong_pairs = np.append(cong_pairs,temp)
		rcvec[cc] = np.mean(congtemp)
		stdvec[cc] = np.std(congtemp)
	return rcvec, stdvec, cong_pairs
######################################################################################
def pathways(bigvec,groups):  # signaling pathway selection
	nc,nrv,npars = bigvec.shape
	rcvec = np.zeros(np.max(groups))
#	stdvec = np.zeros(nc)
	for cc in range(nc):
#		rc = 0  # initialize
#		congtemp = []  # placeholder for cong coef to comp the 
		for ii in range(nrv):
#			for jj in range(ii+1,nrv): # all combinations only once
#				temp = congr(bigvec[cc,ii,:],bigvec[cc,jj,:])
#				rc += temp
#				congtemp = np.append(congtemp,temp)
			temp = set(groups[np.abs(bigvec[cc,ii,:]>1e-6)])
			for jj in range(len(temp)):
				rcvec[temp.pop()-1]+=1
#		rcvec[cc] = np.mean(congtemp)
#		stdvec[cc] = np.std(congtemp)
	return rcvec
###############################################################################
## start here #################################################################
###############################################################################

cancers = ['blca','lgg','brca','cesc','coadread','hnsc','kirc','lihc','luad','lusc','ov','prad','skcm','stad','thca','ucec']			
NC = len(cancers) # how many unique elements
plot_sol = False # plot c-vals for each cancer
plot_box = True # plot box-plotted c-vals
can_error = []  # cancers that experienced an error (cval = 0)

tval_vec, pval_vec = np.zeros(NC), np.zeros(NC)
cong, cong_gr = np.zeros(NC), np.zeros(NC)
cong_std, cong_std_gr = np.zeros(NC), np.zeros(NC) # congruence std
av_lasso, av_grlasso = np.zeros(NC), np.zeros(NC) # mean
std_lasso, std_grlasso = np.zeros(NC), np.zeros(NC)  # c-value std
#c_lasso_vec, c_grlasso_vec = [], []  # save all cvals for later
tt_cong, pp_cong = np.zeros(NC), np.zeros(NC)

pnames=pd.read_csv('results/pathway_names.csv').values  # pathway names
pnames = pnames[:,1:].astype('str')
npath = len(pnames)
	
ntests = 0
npathways = np.zeros(len(pnames)) 
for jj in range(NC):  # all cancers
		can = cancers[jj]

		try:
			cind2 = np.load('./results/cval_single_'+can+'.npy').flatten()
			cind1 = np.load('./results/cval_single_nogroup_'+can+'.npy').flatten()
			  
#		else: # otherwise warning
		except IOError:
			print('No c-value group data for ' + can + ' available.')
			cind1, cind2 = 0, 0
			continue
		else:
			# check if there are any 0 and remove them
			indd = cind2>1e-6
			if (np.sum(indd)!=len(indd)): # not all cvals are above 0
				cind2 = cind2[indd]
				cind1 = cind1[indd]
				can_error.append(can)
				print(str(np.sum(~indd))+' cvals removed from '+can)  # print warning
				
			# do you want to plot the c-values?
			if plot_sol:
				plt.plot(cind1,'*-',cind2,'o-')
				plt.legend(['lasso','group lasso'])
				plt.title(str(can))
				plt.show()

		tt, pp = ttest(cind1,cind2)
#		
#		# save t and p values
		tval_vec[jj] = tt
		pval_vec[jj] = pp
		av_lasso[jj] = np.mean(cind1)
		av_grlasso[jj] = np.mean(cind2)
		std_lasso[jj] = np.std(cind1)
		std_grlasso[jj] = np.std(cind2)
#		c_lasso_vec.append(cind1)
#		c_grlasso_vec.append(cind2)

		# compute average congruence
		if os.path.isfile('./results/beta_single_nogroup_vec_'+can+'.npy'): # is there a file of this name?
			try:
				beta2 = np.load('./results/beta_single_vec_'+can+'.npy')
			except:
				print('fail', can)
			try:
				beta1 = np.load('./results/beta_single_nogroup_vec_'+can+'.npy')
			except:
				print('fail',can)
			# remove failed computations
			if (np.sum(indd)!=len(indd)): # not all cvals are above 0
				beta1 = beta1[:,indd,:]
				beta2 = beta2[:,indd,:]
				
			cong[jj], cong_std[jj], congvec = cong_av(beta1)
			cong_gr[jj], cong_std_gr[jj], congvec_gr = cong_av(beta2)
			groups = pd.read_csv('results/'+can+'_groups.csv').values[:,1:].flatten()
##			cong[jj], cong_std[jj] = pathways(beta1,groups)
			path_gr = pathways(beta2,groups)	
			npathways += path_gr
			ntests += np.sum(indd)
#			print(cong_gr)
			
			# significance of cong
			tt_cong[jj], pp_cong[jj] = ttest(congvec,congvec_gr)
		else: # otherwise warning
			print('No beta data for ' + can + ' available.')
 
# is the congruence difference significant
tt_c,pp_c = ttest(cong,cong_gr)
print('Congruence p-value: '+str(pp_c))

indsort = np.argsort(-npathways)
path = npathways[indsort]/ntests*100 # in percent
pnames_sort = pnames[indsort]
pnames_list = [str(pnames_sort[ii][0]) for ii in range(npath)]
#
#plt.bar(np.arange(npath),path,tick_label = pnames_list)

## plot most common pathways ever
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots()
howmany = 10
cmap = mtplt.cm.get_cmap('PiYG') # spectral / viridis
#color = cmap(jj/ng))
plt.barh(np.arange(howmany),np.flip(path[:howmany]), tick_label = np.flip(pnames_list[:howmany]),
			color = [cmap(ii/howmany) for ii in range(howmany)])
plt.xlabel('%')
#ax.set_xticklabels(pnames_list[:10], rotation = 90)#, ha="right")
plt.title('Most frequently selected pathways',fontsize=20)
plt.savefig('frequency.pdf',bbox_inches='tight')
plt.show()

# plot average cvals + tvals # sorted in descending order
plt.rcParams.update({'font.size': 18})
ind = np.argsort(-av_grlasso) # sorted indices
cancers_sort = [str.upper(cancers[ii]) for ii in ind]  # capitalize

fig, ax = plt.subplots()
#plt.plot(av_lasso[ind], '*-', av_grlasso[ind],'o-')
plt.errorbar(np.arange(NC),av_lasso[ind],std_lasso[ind], marker='o', capsize=3,color = 'cornflowerblue',ecolor = 'lightsteelblue') #linestyle='None',
plt.errorbar(np.arange(NC),av_grlasso[ind],std_grlasso[ind], marker='*', capsize=3,ecolor = 'peachpuff',color = 'lightsalmon') #linestyle='None',
plt.xticks(np.arange(NC), cancers_sort, rotation='vertical') # change tick labels to text
plt.title('Average c-values')
plt.legend(['lasso','group lasso'])

# annotate datapoints with their t-values
for i, txt in enumerate(tval_vec[ind]):
	if pval_vec[ind[i]] < 0.05:   # if the p-value is significant, add an asterisk to the correspondning t-value
		ax.annotate('*', (i,0.02+ av_grlasso[ind[i]]), fontsize=18)
	#else:
	   #ax.annotate(np.round(txt,2), (i, av_grlasso[ind[i]]), fontsize=10)
#plt.savefig('cvals_lasso_vs_grlasso.png',bbox_inches='tight')
plt.savefig('cvals_lasso_vs_grlasso.pdf',bbox_inches='tight')
plt.show()


# plot average congruence
fig, ax = plt.subplots()
#plt.plot(cong[ind], '*-', cong_gr[ind],'o-')
plt.errorbar(np.arange(NC),cong[ind],cong_std[ind], marker='o', capsize=3,color = 'cornflowerblue',ecolor = 'lightsteelblue')#,linestyle='None'
plt.errorbar(np.arange(NC),cong_gr[ind],cong_std_gr[ind], marker='*', capsize=3,ecolor = 'peachpuff',color = 'lightsalmon')#,linestyle='None'
plt.xticks(np.arange(NC), cancers_sort, rotation='vertical') # change tick labels to text
plt.title('Average congruence')
#plt.legend(['lasso','group lasso'],loc='lower left')
#plt.savefig('cong_lasso_vs_grlasso.png',bbox_inches='tight')
# annotate datapoints with their t-values
for i, txt in enumerate(pp_cong[ind]):
	if pp_cong[ind[i]] < 0.05:   # if the p-value is significant, add an asterisk to the correspondning t-value
		ax.annotate('*', (i,0.02+ max(cong[ind[i]],cong_gr[ind[i]])), fontsize=18)
plt.savefig('cong_lasso_vs_grlasso.pdf',bbox_inches='tight')
plt.show()

## plot a typical beta solution
indplot = 42  # which of the 100 runs to plot
plt.plot(beta1[0,indplot,:],'o')#,color='cornflowerblue')
plt.title('Lasso solution for '+str.upper(can))
plt.xlabel('gene index ordered by pathway')
#plt.ylabel('beta')
plt.ylabel(r'$\beta$')
#plt.savefig('beta_lasso_'+can+'.png',bbox_inches='tight')
plt.savefig('beta_lasso_'+can+'.pdf',bbox_inches='tight')
plt.show()

plt.plot(beta2[0,indplot,:],'o')#,color='cornflowerblue')
plt.title('Group lasso solution for '+str.upper(can))
plt.xlabel('gene index ordered by pathway')
#plt.ylabel('beta')
plt.ylabel(r'$\beta$')
#plt.savefig('beta_grlasso_'+can+'.png',bbox_inches='tight')
plt.savefig('beta_grlasso_'+can+'.pdf',bbox_inches='tight')
plt.show()
