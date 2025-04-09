#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    INITIALIZE MODEL

    Script to initialize the model by starting from a set of randomly choosen fractures found with MontpInv software. 

    Last modification: 25/02/12 - Lea

"""

import numpy as np
import random as rd
from os.path import join as pjoin
import os
import shutil


# ----- Prepare files and folders needed for each Markov chain -----

# Obtain the name of the current folder (which depends on the chain)
path_folder = os.getcwd()
name_new_folder = os.path.basename(os.getcwd())
print("***** Hej från Python initialize_fracmodels.py på " + name_new_folder + "! *****")

# Go to the forward model repository
os.chdir('../..')
path_folder = os.getcwd()

if not os.path.isdir('results_chains'):
    os.mkdir('results_chains')

os.chdir('results_chains')

# Check if directory already exists 
if os.path.isdir(name_new_folder):
        print(pjoin('WARNING: ', name_new_folder,' directory already exists! It is removed before starting the new simulation.'))
        shutil.rmtree(name_new_folder)

# Copy common folders between the chains to avoid overwritting and short-circuits
shutil.copytree(pjoin('..','src','input'),pjoin(name_new_folder,'input')) 
shutil.copytree(pjoin('..','src','results'),pjoin(name_new_folder,'results_fwd'))
os.chdir('..')

# ----- Initialize the model -----

# Load a list of fractures obtained with pattern search optimization for three different seeds (MATLAB version)
lib_frac = np.loadtxt(pjoin(path_folder,'src','lib_fractures.txt'))
nb_lib = np.shape(lib_frac)[0]

# Select random fractures from the fracture library file to initialize each Markov chain differently
nb_frac = 5
with open(pjoin(path_folder,'results_chains',name_new_folder,'input','testDFN','DFN.txt'), 'w') as f:
    f.write('DETERMINISTIC1\n')
    f.write(str(int(nb_frac))+"\n")
    for i in range (0,nb_frac):
        nb_i = rd.randint(0,nb_lib-1)
        frac_i = lib_frac[nb_i,:]
        frac = '\t'.join(str(x) for x in frac_i)
        f.write(frac+"\n")
    f.close()

# Save the initial model
path_data=pjoin(path_folder,'results_chains',name_new_folder,'input','testDFN','DFN.txt')
new_path_data = pjoin(path_folder,'results_chains',name_new_folder,'results_fwd','DFN_initial.txt')
shutil.copyfile(path_data,new_path_data)

