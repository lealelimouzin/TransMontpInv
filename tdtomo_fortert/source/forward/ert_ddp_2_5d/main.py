#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import os
from os.path import join as pjoin
import fwdFunctions as fwd

'''
    MAIN

    Main script to run the forward of electrical resistivity data in a context of heterogeneous or fractured media.

    Last modification: 25/02/12 - Lea
'''


if __name__ == "__main__":


    # Obtain the name of the current folder (name depends on the chain)
    
    path_folder = os.getcwd()
    name_new_folder = os.path.basename(os.getcwd())
    # path_folder = '/home/leal/Desktop/Doctorat/Scripts/TransMontpInv/TransMontpInv/tdtomo_fortert/source/forward/ert_ddp_2_5d/results_chains/res_chain_1/'
    # name_new_folder = 'res_chain_1'
    print("***** Demat da Python main.py " + name_new_folder + " ! *****")

    # Set parameters
    expt_name_epm = 'testEPM'                                                # tabVar[0]
    expt_name_dfn = 'testDFN'                                                # tabVar[1]
    data_xx_yy = np.loadtxt(pjoin('..','..','fielddata','data_xxyy_R4_noratio.txt')) # tabVar[2]
    #data_xx_yy = np.loadtxt(pjoin('fielddata', 'data_xxyy_R4_noratio.txt'))  # tabVar[2]
    nb_electrodes = 128                                                      # tabVar[3]
    path_results = 'results_fwd'                                             # tabVar[4]
    ratioData = pd.read_csv(pjoin('..','..','..','..','..','bin','data.dat'), delimiter='\t', header=None,skiprows=1)
    #ratioData = pd.read_csv(pjoin('..', '..', '..', 'bin', 'data.dat'), delimiter='\t', header=None,skiprows=1)
    ratioData = ratioData.to_numpy()                                         # tabVar[5]
    tabVar = [expt_name_epm, expt_name_dfn, data_xx_yy, nb_electrodes, path_results, ratioData]

    # Forward EPM
    if not os.path.exists(pjoin(tabVar[4], "rho_a_EPM.txt")):
        rho_a_EPM = fwd.fwdEPM(tabVar)
        np.savetxt(pjoin(tabVar[4],"rho_a_EPM.txt"), rho_a_EPM, delimiter='\t')

    # Forward DFN
    rho_a_DFN = fwd.fwdDFN(tabVar)
    np.savetxt(pjoin(tabVar[4], 'rho_a_DFN.txt'),rho_a_DFN, delimiter='\t')

   