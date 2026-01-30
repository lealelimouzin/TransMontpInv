#!/bin/bash

# ********************************************************************* #
# Script to run multiple mksamples jobs on the ECDF cluster
# Copyright (C) 2017  Erica Galetti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ********************************************************************* #


##############################################################################################

# Initialise the environment
. /etc/profile
. /etc/profile.d/modules.sh

# Load the necessary modules
#module load intel
#module load openmpi/1.6.5
#module load igmm/mpi/gcc/mpich/3.1.4

# Set work and bin directories
PROJECT=/home/egaletti/erica/ForTERT
#PROJECT=/home/leal/Desktop/Doctorat/Scripts/transdimension/TERT/tdtomo_fortert
bin_dir=${PROJECT}/RUN_TOMO/bin_model2b_ptrand
out_dir=${PROJECT}/TOMO_RESULTS/model2b/OUT_ptrand

# Define some parameters
number_of_runs=10               # Number of runs in total
run_number=1                    # Start from this run number
number_of_chains=32             # How many parallel chains?
runtime=16:00:00                # Run time for one run of mksamples
sleep_time=10m                  # Sleep time between checks for job completion
max_fails=10                    # Maximum number of consecutive times for code failure

##############################################################################################

# Create results directory
mkdir -p ${out_dir}

# Move to bin directory
cd ${bin_dir}

### ------------------------- OUTPUT FILES ----------------------- ###
# Get the names of the output files from mksamples.in
samples_file=`awk '{ if (NR==11) print $1}' mksamples.in`
ncell_file=`awk '{ if (NR==12) print $1}' mksamples.in`
noise_file=`awk '{ if (NR==13) print $1}' mksamples.in`
misfits_file=`awk '{ if (NR==14) print $1}' mksamples.in`
acceptance_file=`awk '{ if (NR==15) print $1}' mksamples.in`
temperatures_file=`awk '{ if (NR==16) print $1}' mksamples.in`
residuals_file=`awk '{ if (NR==18) print $1}' mksamples.in`
echo
echo "Output files are:"
echo "   - ${samples_file}"
echo "   - ${ncell_file}"
echo "   - ${noise_file}"
echo "   - ${misfits_file}"
echo "   - ${acceptance_file}"
echo "   - ${temperatures_file}"
echo "   - ${residuals_file}"
echo
### -------------------------------------------------------------- ###


### ----------------------- PREPARE FUNCTION --------------------- ###
# This function checks if a job has finished executing
function checkjob {
        
        local job_code=$1
        qstat > job_status.txt
        grep "${job_code} " job_status.txt > line1.txt
        awk -v jn=`echo ${job_code}` '{if ($1==jn) print $0}' line1.txt > line.txt
        local job_done=0

        if [ -s line.txt ] ; then
                local job_status=`awk '{print $5}' line.txt`
        else
                #echo Job ${job_code} does not exist
                local job_done=1
        fi
        echo ${job_done}
}
### -------------------------------------------------------------- ###


### ------------ PREPARE SUBMISSION SCRIPT FOR MPI JOB ----------- ###
echo Preparing submission script...
echo

echo "#!/bin/sh" > jobsub.sh
echo >> jobsub.sh

echo "# Grid Engine Options" >> jobsub.sh
echo "#$ -N mksamples" >> jobsub.sh
echo "#$ -cwd" >> jobsub.sh
echo "#$ -l h_rt=${runtime}" >> jobsub.sh
echo "#$ -l h_vmem=2G" >> jobsub.sh
echo "#$ -pe mpi ${number_of_chains}" >> jobsub.sh
echo "#$ -R y" >> jobsub.sh

#echo "#$ -P geos_andrewcurtis" >> jobsub.sh
echo >> jobsub.sh

echo "# Initialise the environment" >> jobsub.sh
echo . /etc/profile >> jobsub.sh
echo . /etc/profile.d/modules.sh >> jobsub.sh
echo >> jobsub.sh

echo "# Load the intel and MPI modules" >> jobsub.sh
echo module load intel >> jobsub.sh
echo module load openmpi/1.6.5 >> jobsub.sh
echo module load igmm/mpi/gcc/mpich/3.1.4 >> jobsub.sh
echo >> jobsub.sh

echo "# Run executable" >> jobsub.sh
echo "mpiexec -np \$NSLOTS ./mksamples 1" >> jobsub.sh
### -------------------------------------------------------------- ###


### ----- RUN THE INVERSIONS AS SUCCESSIVE RUNS OF MKSAMPLES ----- ###
echo Running inversion...
echo
### -------------------------------------------------------------- ###


### ---------------------- SUBMIT MKSAMPLES ---------------------- ###
job_name=mksamples

# Loop over the total number of runs
count=1
while [ ${run_number} -le ${number_of_runs} ] ; do
        
        rm -f completed?.txt
        
        # Modify the submission script to include the run number
        awk -v jname=`echo ${job_name}` -v nm=`echo ${num}` -v nn=`echo ${run_number}` '
                { if(NR==4) print "#$ -N " jname ;
                else if(NR==21) print "mpiexec -np $NSLOTS ./" jname " " nn " ";
                else print $0 }
          ' jobsub.sh > jobsub1.sh
        mv jobsub1.sh jobsub.sh
        
        # Submit job
        if [ ${count} -le ${max_fails} ] ; then # if the number of failed attempts is <= 10
                qsub jobsub.sh > sub.txt
        else
                echo
                echo "Code crashed too many times on same iteration!!"
                exit
        fi
        more sub.txt
        
        # Get job number
        job_number=`awk '{print $3}' sub.txt`
        echo ${job_name} job number is ${job_number} - ${run_number}/${number_of_runs}
        
        # Check if the job has completed at intervals of sleep_time
        completed=0
        while [[ ${completed} == "0" ]] ; do
                sleep ${sleep_time}
                completed=`checkjob ${job_number}`
        done
        
        # Check if the job has completed successfully
        qdel ${job_number}
        grep "Job finished on " ${job_name}.o${job_number} > completed.txt
        
        if  [ -s completed.txt ] ; then # if the job has completed successfully, add output files to compressed folder before starting a new run
        
                # Create directory for outputs
                job_directory=${job_name}__run_${run_number}
                mkdir ${bin_dir}/${job_directory}
                zip_directory=${job_directory}.zip
                
                # Move/copy outputs into job_directory
                mv ${job_name}.*${job_number} ${job_directory}
                cp ${samples_file}* ${job_directory}
                cp ${ncell_file}* ${job_directory}
                cp ${noise_file}* ${job_directory}
                cp ${misfits_file}* ${job_directory}
                cp ${acceptance_file}* ${job_directory}
                cp ${temperatures_file}* ${job_directory}
                cp ${residuals_file}* ${job_directory}
                cp last.tmp* ${job_directory}
                cp parameters.dat parameters.txt ${job_directory}
                zip -r ${zip_directory} ${job_directory}
                mv ${zip_directory} ${out_dir}
                rm -r ${job_directory}
                
                # Echo message and proceed to next iteration
                echo "Inversion ${run_number}/${number_of_runs} OK (job number ${job_number})"
                run_number=`echo ${run_number} + 1 | bc`
                count=1
                
        else # if the job crashed or got interrupted, echo message and maybe try again
                
                echo "Inversion ${run_number}/${number_of_runs} FAILED (job number ${job_number})"
                
                # If the code crashed because there is a mismatch in the number of lines of output files, stop
                grep "Lengths of output files are inconsistent" ${job_name}.o${job_number} > completed1.txt
                
                # If the code stopped because the desired chain length has been reached, stop
                grep "Maximum number of samples " ${job_name}.o${job_number} > completed2.txt
                
                # If the code stopped for any other error, stop
                grep "TERMINATING PROGRAM" ${job_name}.o${job_number} > completed3.txt
                
                # Create directory for outputs
                job_directory=${job_name}__run_${run_number}__FAIL
                mkdir ${bin_dir}/${job_directory}
                zip_directory=${job_directory}.zip
                
                # Move/copy outputs into job_directory
                mv ${job_name}.*${job_number} ${job_directory}
                zip -r ${zip_directory} ${job_directory}
                mv ${zip_directory} ${out_dir}
                rm -r ${job_directory}
                
                if  [ -s completed1.txt ] ; then
                        echo
                        echo "Code crashed because lengths of output files are inconsistent!!"
                        exit
                elif [ -s completed2.txt ] ; then
                        echo
                        echo "Code stopped because desired number of samples has been reached!!"
                        exit
                elif [ -s completed3.txt ] ; then
                        echo
                        echo "Code stopped because of an error!! Check file ${job_name}.*${job_number}"
                        echo "in compressed folder ${out_dir}/${zip_directory}!!"
                        exit
                fi
                
                count=`echo ${count} + 1 | bc`
                
        fi
done
rm completed?.txt line1.txt line.txt job_status.txt
### -------------------------------------------------------------- ###

