#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 11:53:17 2021

@author: gabriela
"""
# OBS: the files are padded unnecessarily so that the format fits the one 
# required by Cox_GD_group_multitask.py
# The gene expressions (matrix x) are already loged, so one has to manually 
# toggle the corresponding line in the prepare_data function
# in Cox_GD_group_multitask.py

from __future__ import print_function
from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mtplt
from mpl_toolkits import mplot3d
import copy
import time
import sys
import os.path
import torch
import torch.optim as optim
import torch.optim.lr_scheduler

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds.\n" %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)

def get_data(cancer_type):
	x=pd.read_csv('surv_files/'+cancer_type+'_x.csv', header=None, low_memory=False).values
	gnames = x[0,1:].astype('str')		# gene names
	os_months=pd.read_csv('surv_files/'+cancer_type+'_os_months.csv',dtype='float32').values
	os_status=pd.read_csv('surv_files/'+cancer_type+'_os_status.csv',dtype='int').values
	pathway_names=pd.read_csv('surv_files/pathway_names.csv').values  # pathway names
	# remove first column = numberin line
	os_months = os_months[:,1:]
	os_status = os_status[:,1:]
	x = x[1:,1:].astype('float32')  # remove first row (names)
	pathway_names = pathway_names[:,1:].astype('str')
	# flatten into a list (not a list of single lists)
	os_months = os_months.flatten()
	os_status = os_status.flatten()
	indices = np.argsort(os_months)
	x_ordered = x[indices][:]
	os_status_ordered = os_status[indices] #np.array(os_status)[indices.astype(int)]
	os_months_ordered = os_months[indices] #np.array(os_months)[indices.astype(int)]
	return x_ordered,os_months_ordered,os_status_ordered, gnames, pathway_names
###############################################################################
def group_betas(beta):
	nb, npa = beta.shape  # nr of betas, nr of pathways
	gr = np.zeros((nb,1))
	for ii in range(npa):
		gr[np.abs(beta[:,ii]) > 1e-6] = ii+1
	if np.any(gr==0):
		raise Warning('Some betas not assigned.')
	return gr
###############################################################################
def npstnd(x_ordered):
	# Standardize
	meansx = np.mean(x_ordered, axis=0)
	stdsx = np.std(x_ordered, axis=0)
	stdin = (stdsx>0)   # only consider stdx > 0
	x_ordered[:,stdin] = (x_ordered[:,stdin] - meansx[stdin]) / stdsx[stdin]
	return x_ordered
###############################################################################
def clean(beta,gnames,tgenes):
#	first remove 0
	ind = []
	target = set(tgenes)
	for ii in range(len(gnames)):
		# don't even consider if all beta is 0
		if np.linalg.norm(beta[ii,:])<1e-6:
			ind.append(ii)
			continue
		elif gnames[ii].upper() in target:
			continue
#			print('yada')
		else:
			try:
				D[gnames[ii].upper()] in target
				continue
			except:
				print('Removing '+gnames[ii].upper())
#				beta = np.delete(beta,ii,0)
#				gnames = np.delete(gnames,ii)
				ind.append(ii)

	beta = np.delete(beta,ind,0)
	gnames = np.delete(gnames,ind)
	
	return beta, gnames
###############################################################################
def censor(score,sc):
	ind = np.argsort(score) # index of sorting
	D = dict()
	for count,value in enumerate(ind):
		D[count] = value  # fill dictionary
	score_sort = np.sort(score)
	sc_sort = sc[ind]
	score1 = np.zeros(len(score))
#	for ii in range(len(score)):
#		if sc_sort[ii]==True:
#			score_sort[ii] = score_sort[ii] + 100
#			score1[D[ii]] = score_sort[ii]
	for ii in range(len(score)):
		if sc_sort[ii]==True:
			ind_new = min(ii+int(np.round((len(score)-ii)*np.random.rand(1))),len(score)-2)  # new position
			value_new = (score_sort[ind_new]+score_sort[ind_new+1])/2 # average
#			score1[D[ii]] = score_sort[ii] + 100 #np.abs(np.median(score)*np.random.rand(1))
			score1[D[ii]] = value_new
		else:
			score1[D[ii]] = score_sort[ii]
	return score1
########################################################################################		
class can( object ):   # create class to hold cancer attributes
	def __init__( self ):
		self.name = str
		self.x = float   # ordered
		self.os_months = float  # ditto
		self.os_status = bool   # ditto
		self.gnames = str
		self.pnames = str
		self.groups = int
		self.n_samples = int
		self.n_genes = int
###############################################################################
### MAIN

switch = 2 # 1 or 2 cancer	
trim = False # trim to size of the significant genes	

pat1 = 300 # nr of toy patients
pat2 = 200
		
# rename obsolete genes
D = dict()
D['C9ORF91'] = 'TMEM268'
D['EIF2C3'] = 'AGO3'
D['KBTBD10'] = 'KLHL41'
D['KIAA0141'] = 'DELE1'
D['KIAA0247'] = 'SUSD6'
D['KIAA0564'] = 'VWA8'
D['KIAA1467'] = 'FAM234B'
D['RGAG4'] = 'RTL5'
D['TENC1'] = 'TNS2'
D['TMEM27'] = 'CLTRN'
D['C11ORF82'] = 'DDIAS'
D['C15ORF42'] = 'TICRR'
D['C16ORF59'] = 'TEDC2'
D['C1ORF135'] = 'AUNIP'
D['C1ORF38'] = 'THEMIS2'
D['FAM54A'] = 'MTFR2'
D['FAM64A'] = 'PIMREG'
D['GSG2'] = 'HASPIN'
D['MFI2'] = 'MELTF'
D['MLF1IP'] = 'CENPU'
D['TMEM48'] = 'NDC1'
D['C10ORF137'] = 'EDRF1'
#D['C10ORF2'] = 'TWNK' # doesn't exist
D['C17ORF90'] = 'OXLD1'
D['C21ORF59'] = 'CFAP298'
D['CCDC41'] = 'CEP83'
D['CECR5'] = 'HDHD5'
D['CHCHD8'] = 'COA4'
D['CTPS'] = 'CTPS1'
D['GUCY1A3'] = 'GUCY1A1'
#D['HN1'] = 'JPT1'  # same
D['KIAA0020'] = 'PUM3'
D['KIAA0664'] = 'CLUH'
D['KIAA0922'] = 'TMEM131L'
#D['LOC388796'] = 'SNHG17'
D['MINA'] = 'RIOX2'
D['MTERFD1'] = 'MTERF3'
D['NHP2L1'] = 'SNU13'
D['QTRTD1'] = 'QTRT2'
D['SDCCAG3'] = 'ENTR1'
D['SKIV2L2'] = 'MTREX'
D['SMCR7L'] = 'MIEF1'
D['TMEM5'] = 'RXYLT1'
D['WBSCR22'] = 'BUD23'
D['CD3EAP'] = 'POLR1G'
D['DARS'] = 'DARS1'
D['H2AFZ'] = 'H2AZ1'
D['HARS'] = 'HARS1'
D['HIST1H2AC'] = 'H2AC6'
D['IARS'] = 'IARS1'
D['KARS'] = 'KARS1'
D['NARS'] = 'NARS1'
#D['DSCAM'] = 'DSCAM' # not found
#D['NHLH2'] = 'NHLH2' # not found
#D['SCGB2A2'] = 'SCGB2A2' # not found
#D['TRIML2'] = 'TRIML2' # not found

# load betas
x = pd.read_csv('dataforG.csv', header=None, low_memory=False).values
pnames = x[0,1:].astype('str')		# pathway names
gnames = x[1:,0].astype('str')      # gene names
beta = x[1:,1:].astype('float32')   # remove first row (names)

if switch == 1:
	# load cancer matrix
	ca1 = can()
	ca1.name = 'coadread'   # save c name
	ca1.x,ca1.os_months,ca1.os_status,ca1.gnames,ca1.pnames = get_data(ca1.name)
	ca1.n_genes = ca1.x.shape[1]
	ca1.x = np.log(ca1.x+1)  # log
	ca1.groups = pd.read_csv('surv_files/'+ca1.name+'_groups.csv').values[:,1:].flatten()

	beta, gnames = clean(beta,gnames,ca1.gnames)

if switch == 2:
	# load cancer matrix
	ca2 = can()
	ca2.name = 'stad'   # save c name
	ca2.x,ca2.os_months,ca2.os_status,ca2.gnames,ca2.pnames = get_data(ca2.name)
	ca2.n_genes = ca2.x.shape[1]
	ca2.x = np.log(ca2.x+1)
	ca2.groups = pd.read_csv('surv_files/'+ca2.name+'_groups.csv').values[:,1:].flatten()

	beta, gnames = clean(beta,gnames,ca2.gnames)

# remove unused genes
if switch == 1 and trim == True:
	gr = []
	xcol1 = []
	gnamesnew = []
	beta1, beta21 = [], []
	for jj in range(len(pnames)):
		for ii in range(len(gnames)):
			if np.abs(beta[ii,jj])>1e-6: # if nonzero coeff
				gr = np.concatenate([gr, [jj]])  # add group tag
				if len(gnamesnew) == 0:
					gnamesnew = [str(gnames[ii])]
				else:
					gnamesnew.append(str(gnames[ii]))
				
				if jj == 0:
					beta1 = np.concatenate([beta1,[beta[ii,0]]])
#					beta1 = np.concatenate([beta1,[1]])
				else:
					beta1 = np.concatenate([beta1,[0]])
				if jj == 1:
					beta21 = np.concatenate([beta21,[beta[ii,1]]])
#					beta21 = np.concatenate([beta21,[1]])
				else:
					beta21 = np.concatenate([beta21,[0]])
					
				temp = np.where(np.char.upper(ca1.gnames) == np.char.upper(gnames[ii]))
				if len(temp[0])>0:
					a = temp[0][0]  # uppercase # unfold the tuple
				else:  # try the dictionary
					a = np.where(np.char.upper(ca1.gnames) == D[str(np.char.upper(gnames[ii]))])[0][0]
#				if ii==0:
				if len(xcol1) == 0: # no elements yet
					xcol1 = ca1.x[:,a]
				else:
					xcol1 = np.vstack((xcol1,ca1.x[:,a]))
	
	ca1.x = np.transpose(xcol1)
#	ca1.gnames = gnames	
	ca1.gnames = np.array(gnamesnew)
	ca1.n_genes = ca1.x.shape[1]
	ca1.groups = gr
	
if switch == 2 and trim == True:
	gr = []
	xcol2 = []
	gnamesnew = []
	beta22, beta3 = [], []
	for jj in range(len(pnames)):
		for ii in range(len(gnames)):
			if np.abs(beta[ii,jj])>1e-6: # if nonzero coeff
				gr = np.concatenate([gr, [jj]])  # add group tag
				if len(gnamesnew) == 0:
					gnamesnew = [str(gnames[ii])]
				else:
					gnamesnew.append(str(gnames[ii]))
				
				if jj == 1:
					beta22 = np.concatenate([beta22,[beta[ii,1]]])
#					beta22 = np.concatenate([beta22,[1]])
				else:
					beta22 = np.concatenate([beta22,[0]])
				if jj == 2:
					beta3 = np.concatenate([beta3,[beta[ii,2]]])
#					beta3 = np.concatenate([beta3,[1]])
				else:
					beta3 = np.concatenate([beta3,[0]])
					
				temp = np.where(np.char.upper(ca2.gnames) == np.char.upper(gnames[ii]))
				if len(temp[0])>0:
					a = temp[0][0]  # uppercase # unfold the tuple
				else:  # try the dictionary
					a = np.where(np.char.upper(ca2.gnames) == D[str(np.char.upper(gnames[ii]))])[0][0]
#				if ii==0:
				if len(xcol2) == 0: # no elements yet
					xcol2 = ca2.x[:,a]
				else:
					xcol2 = np.vstack((xcol2,ca2.x[:,a]))
	
	ca2.x = np.transpose(xcol2)
#	ca2.gnames = gnames	
	ca2.gnames = np.array(gnamesnew)
	ca2.n_genes = ca2.x.shape[1]
	ca2.groups = gr

##

if switch == 1:
	### Generate data with the same distributions
	tic()
	mu1 = np.mean(ca1.x, axis = 0)  # mean
	sig1 = np.cov(ca1.x, rowvar = False) # covariance matrix
#	sig1 = np.eye(ca1.n_genes)
	xx1 = np.abs(np.random.multivariate_normal(mu1,sig1,pat1))
	xx1 = npstnd(xx1)
	toc()

if switch == 2:
	tic()
	mu2 = np.mean(ca2.x, axis = 0)  # mean
	sig2 = np.cov(ca2.x, rowvar = False) # covariance matrix
#	sig2 = np.eye(ca2.n_genes)
	xx2 = np.abs(np.random.multivariate_normal(mu2,sig2,pat2))
	xx2 = npstnd(xx2)
	toc()

if switch == 1:
	if trim == False:
		beta1 = np.zeros(ca1.n_genes)
		beta21 = np.zeros(ca1.n_genes)
		for ii in range(len(gnames)):
			temp = np.where(np.char.upper(ca1.gnames) == np.char.upper(gnames[ii]))
			if len(temp[0])>0: # if the gene has been found
				a = temp[0][0]  # uppercase # unfold the tuple
			else:  # dictionary to translate the gene name
				a = np.where(np.char.upper(ca1.gnames) == D[str(np.char.upper(gnames[ii]))])[0][0]
			beta1[a] = beta[ii,0]
			beta21[a] = beta[ii,1]
	
	score1 = np.matmul(xx1,beta1+beta21)

if switch == 2:
	if trim == False:
		beta22 = np.zeros(ca2.n_genes)
		beta3 = np.zeros(ca2.n_genes)
		for ii in range(len(gnames)):
			temp = np.where(np.char.upper(ca2.gnames) == np.char.upper(gnames[ii]))
			if len(temp[0])>0: # if the gene has been found
				a = temp[0][0]  # uppercase # unfold the tuple
			else:  # dictionary to translate the gene name
				a = np.where(np.char.upper(ca2.gnames) == D[str(np.char.upper(gnames[ii]))])[0][0]
#			beta22[a] = beta[ii,1]
			beta22[a] = beta[ii,0]
			beta3[a] = beta[ii,2]
			
	score2 = np.matmul(xx2,beta22+beta3)

## all deceased
#status1 = np.ones((pat1+1,2))  # add junk
#status2 = np.ones((pat2+1,2))
	
#status1 = np.random.randint(2, size=(pat1+1,2))  # pateints randomly event/censored
#status2 = np.random.randint(2, size=(pat2+1,2))
#
# p = 10* proportion deceased, 0<p<10
p = 7	
status1 = np.array(list(map(int,np.random.randint(10,size=(pat1+1))<p)))
status1 = np.hstack((np.zeros((status1.shape[0],1)),np.expand_dims(status1,axis=1)))
status2 = np.array(list(map(int,np.random.randint(10,size=(pat2+1))<p)))
status2 = np.hstack((np.zeros((status2.shape[0],1)),np.expand_dims(status2,axis=1)))

# adjust censored patients: the time of censoring should be less than death
if switch == 1:
	sc = status1[1:,1] == 0  # which are censored
#	score1[sc] = score1[sc] + np.abs(np.median(score1)*np.random.rand(np.sum(sc))) # increase their score
	score1 = censor(score1,sc)

if switch == 2:
	sc = status2[1:,1] == 0  # which are censored
#	score2[sc] = score2[sc] + np.abs(np.median(score2)*np.random.rand(np.sum(sc)))
	score2 = censor(score2,sc)

if switch == 1:
	xx1 = np.vstack((ca1.gnames,xx1))  # x already loged! Needs to be toggled in Cox_GD_group_multitask.py
	xx1 = np.hstack((np.zeros((xx1.shape[0],1)),xx1))
	
	score1 = np.append(np.array((0)),-score1)
	score1 = np.hstack((np.zeros((score1.shape[0],1)),np.expand_dims(score1,axis=1)))
	
	groups1 = np.append(np.array((0)),ca1.groups)
	groups1 = np.hstack((np.zeros((groups1.shape[0],1)),np.expand_dims(groups1,axis=1)))
	
	np.savetxt('toy_'+ca1.name+'_x.csv',xx1, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca1.name+'_os_months.csv',score1, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca1.name+'_os_status.csv',status1, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca1.name+'_groups.csv',groups1, delimiter=",",  fmt="%s")

if switch == 2:
	xx2 = np.vstack((ca2.gnames,xx2))  # x already loged! Needs to be toggled in Cox_GD_group_multitask.py
	xx2 = np.hstack((np.zeros((xx2.shape[0],1)),xx2))
	
	score2 = np.append(np.array((0)),-score2)
	score2 = np.hstack((np.zeros((score2.shape[0],1)),np.expand_dims(score2,axis=1)))
	
	groups2 = np.append(np.array((0)),ca2.groups)
	groups2 = np.hstack((np.zeros((groups2.shape[0],1)),np.expand_dims(groups2,axis=1)))
	
	np.savetxt('toy_'+ca2.name+'_x.csv',xx2, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca2.name+'_os_months.csv',score2, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca2.name+'_os_status.csv',status2, delimiter=",",  fmt="%s")
	np.savetxt('toy_'+ca2.name+'_groups.csv',groups2, delimiter=",",  fmt="%s")
