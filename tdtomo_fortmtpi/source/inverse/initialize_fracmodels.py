#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Program under the GNU General Public License <http://www.gnu.org/licenses/>.

"""
    INITIALIZE MODEL

    Script to initialize the model for each Markov chain by starting from a set of randomly choosen pathways. 

    Last modification: December 2025 - Léa Lelimouzin
"""

import numpy as np
import random as rd
from os.path import join as pjoin
import os
import shutil


# ---------- Get the name of the current folder (which depends on the chain) ----------

path_folder0 = os.getcwd()
name_new_folder = os.path.basename(os.getcwd())
print("***** Hej från Python initialize_fracmodels.py på " + name_new_folder + "! *****")


# ---------- Prepare forward model directory for each Markov chain ----------

os.chdir('../')
path_folder = os.getcwd()

# Copy common folders between the chains to avoid overwritting and short-circuits
shutil.copytree(pjoin('..','..','source','forward','ert_ddp2.5d','src','input'),pjoin(name_new_folder,'input')) 
shutil.copytree(pjoin('..','..','source','forward','ert_ddp2.5d','src','results'),pjoin(name_new_folder,'results_fwd'))
# If the C++ code is compiled as an executable, uncomment the following two lines:
#shutil.copytree(pjoin('..','src','electrical_resistivity'),pjoin(name_new_folder,'electrical_resistivity'))
#shutil.copytree(pjoin('..','src','output'),pjoin(name_new_folder,'output')) 


# ---------- Initialize the model ----------

# Load a list of fractures obtained with pattern search optimization for three different seeds (MATLAB version)
lib_frac = np.loadtxt(pjoin('..','lib_fractures.txt'))
nb_lib = np.shape(lib_frac)[0]

# Select random fractures from the fracture library file to initialize each Markov chain differently
nb_frac = 5
#rd.seed(10) # Set the seed for reproductibility. This line can be comented.
with open(pjoin(name_new_folder,'input','testDFN','DFN.txt'), 'w') as f:
    f.write('DETERMINISTIC1\n')
    f.write(str(int(nb_frac))+"\n")
    for i in range (0,nb_frac):
        nb_i = rd.randint(0,nb_lib-1)
        frac_i = lib_frac[nb_i,:]
        frac = '\t'.join(str(x) for x in frac_i)
        f.write(frac+"\n")
    f.close()


# ---------- Save the initial model ----------

path_data=pjoin(name_new_folder,'input','testDFN','DFN.txt')
new_path_data = pjoin(name_new_folder,'results_fwd','DFN_initial.txt')
shutil.copyfile(path_data,new_path_data)

os.chdir(path_folder0)
