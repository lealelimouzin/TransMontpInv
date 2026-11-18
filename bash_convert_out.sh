#!/bin/sh

# Move to the correct folder
cd '/home/leal/Desktop/Doctorat/Data/Data_Farum/TransMontpInv_E6789/26-01-19_TERT_angle_tmp/' # name of the foler containing the output files

nb_chain=6 # number of chains run in parallel
path_chain2txt='/home/leal/Desktop/Doctorat/Scripts/TransMontpInv/TransMontpInv/tdtomo_fortert/bin/chain2txt' # path towards the executable
#path_chain2txt='/home/leal/Desktop/Doctorat/Scripts/TransMontpInv/TransMontpInv/tdtomo_fortert/bin/chain2txt_xzl'

# Convert files misfit
for n in $(seq 1 $nb_chain)
do
filename='misfit.out'
${path_chain2txt} ${filename}${n}
done

# Convert files ncells
for n in $(seq 1 $nb_chain)
do
filename='ncells.out'
${path_chain2txt} ${filename}${n}
done


# Convert files acceptance rate
for n in $(seq 1 $nb_chain)
do
filename='aratios.out'
${path_chain2txt} ${filename}${n}
done

# Convert files samples
for n in $(seq 1 $nb_chain)
do
filename='samples.out'
${path_chain2txt} ${filename}${n}
done

# Convert files sigmas
for n in $(seq 1 $nb_chain)
do
filename='sigmas.out'
${path_chain2txt} ${filename}${n}
done

# Convert files residuals
for n in $(seq 1 $nb_chain)
do
filename='residuals.out'
${path_chain2txt} ${filename}${n}
done

# Convert files temper
filename='temper.out'
${path_chain2txt} ${filename}

# Generate files related to the processing of the samples: WARNING: not implemented for TransMontpInv
#${procsamples}




