# TransMontpInv

## Description

TransMontpInv is a software to invert electrical resistivity data based on a discrete dual porosity approach. It is a combination of the already existent codes from Galetti & Curtis (2018) and the forward model adapted from Caballero Sanz et al. (2017). The two codes have been merged and adapted to support the geometry of thin pathways such as fractures.

## Table of Contents

- [Contributions](#contributions)
- [Prepare the working environment to run the code](#prepare-the-working-environment-to-run-the-code)
- [Set files](#set-files)
- [Compile](#compile)

## Contributions

### Main code: origins

The main code has been developed by Erica Galetti under the name `TERT - Transdimensional Electrical Resistivity Tomography`. The original codes are deposited following this link: . 

### Main code: modifications

Galetti & Curtis (2018) codes have been modified to support thin pathways structures.

Contacts:
+ lea
+ cedric
+ delphine

### Forward model

### Subroutine to clip the fractures

The subroutine `box_segment_clip_2d` in `inverse\chains.f90` has been developed by John Bukardt and is extracted from his package `GEOMETRY`. It clips the lines (here the fractures) that are outside the boundaries of a rectangular domain. The computer code is under the GNU LGPL license. The complete library and details on the project are available on the web site of the [geometry project](https://people.math.sc.edu/Burkardt/f_src/geometry/geometry.html).

## Prepare the working environment to run the code

Create a Python environment to locally installed the required packages.

An option is to use miniconda as we did. Here are the steps to create a new environment named "tert" and install the packages:

- Open a terminal window
- Install [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions)
- Create a new environment: `conda create -n tert`
- Enter the environment: `conda activate tert`
- Install the packages: ̀`conda install --yes --file requirements.txt` (the package installer "pip" can also be used instead: `pip install -r requirements.txt`)
- To close the environment: `conda deactivate`

Another option to create a Python environment can be to use Python "venv" package:

- Open a terminal window
- Check if "venv" is already installed: `find venv`
- If not, install "venv": `pip install venv`
- Go to the ert_ddp_2_5d folder in the project: `cd tdtomo_fortert/source/forward/ert_ddp_2_5d`
- Create a virtual environment: `python3 -m venv env_python`
- Activate the environment: `source env_python/bin/activate`
- Install the dependencies in the virtual environment: `pip install -r requirements.txt`
- To close the environment: `source deactivate`

## Set files

In the file `tdtomo_fortert/source/main/rjtomo/chain2txt.f90`, change line 23 ("opt" variable) to indicate if there is inversion of angle/aperture.



## Compile

Once the working environment is created, the code needs to be compiled in the environment from the terminal:

- To compile the C++ forward code as a library:
    - Go to: `cd tdtomo_fortert/source/forward/ert_ddp_2_5d/src/electrical_resistivity/Release`
    - Make sure everything is cleaned: `make clean`
    - Compile: `make shared`
    - In `initialize_fracmodels.py`, comment lines 45 and 46
    - In `main.py` change optCompile to `1`
- To compile the C++ forward code as an executable (text files):
    - Go to: `cd tdtomo_fortert/source/forward/ert_ddp_2_5d/src/electrical_resistivity/Release`
    - Make sure everything is cleaned: `make clean`
    - Compile: `make`
    - In `initialize_fracmodels.py`, uncomment lines 45 and 46
    - In `main.py` change optCompile to `2`
- To compile the Fortran code:
    - Go to the compile directory: `cd tdtomo_fortert/source/compile`
    - Make sure everything is cleaned: `make clean`
    - Compile: `make` 

If using a computing cluster, one can execute the bash script "bash_tert" to run the code. If not, go to the "bin" directory and enter the command line: `mpirun -n nbChains ./mksamples nbJob` where instead of "nbChains" write the number of chains launch in parallel and replace "nbJob" with the number of job.

## Test call Python in Fortran

compile: gfortran test.f90 -o test
run: ./test

If mpi, compile: mpif90 test.f90 -o test
run: mpirun -n 4 ./test
mpirun -n 2 ./mksamples 1 > screenDisplay.txt

If problem with mpi.h file, compile as followed: gfortran -I/usr/lib/x86_64-linux-gnu/openmpi/include test.f90 -o test

## Fortran

By default, the naming initially used by Galletti and Curtis (2017) is kept in the scripts. The correspondances to the forward model ert_ddp_2_5d is indicated in comment in the concerned files.

## Test run with ert_ddp_2_5d forward model

Test with 3 fractures for a start, and 3 markov chains.

## Alpha component adapted to fracture

alpha_component = l_proposed^2 / (2*birth_change^2)

instead of (rho_proposed - rho_current)^2 / (2*birth_change^2)

## Temperring

If more than 4 Markov chains are used with two of them having a temperature greater than 1°C, the sequence of temperatures tested is calculated automatically with a linspace subroutine from 2°C to 8°C (here the max indicated in mkasamples.in) with a step corresponding to the number of temperature greater than 1°C.

## Noise in the input data

The noise is estimated with the propagation of uncertainties formula based on the standard deviation given for each measurement and adapted to a ratio implementation.

We implemented the noise option "noise_frac" corresponding to the already existent "noise_lambda" to which the standard deviation of the homogeneous matrix (no fracture) is added. As only fractures are modelled with ert\_ddp\_2.5d model, the standard deviation enables to correct the effect of the matrix present in the measured data.

## Python post-processing

Simplifications have been done in the scripts.

## Run on the cluster

change path in mksamples.in

In source/compile/Makefile, change:

- lib_dir_local = /usr/lib64
- inc = ../../../COMMON/source/tools/superlu/SRC

Make a link to libsuperlu: ln -s /usr/lib64/libsuperlu.so.5.1 libsuperlu.so

## Convert unformatted files to text files

To convert the results files from .out to .txt, use the function `chain2txt.f90` which executable is in the `bin` directory. If the optimisation of the angle or aperture is selected, dimension of some variables need to be changed in the code. To do so, uncomment all the lines having the note "! Inversion of angle/aperture" and comment all the lines with the note "! NO inversion of angle/aperture".

## Common errors and how to handle them

### "PERFORM" not found

First, check if the library `PERFORM.so` has been genrated in the folder `tdtomo_fortert/source/forward/ert_ddp_2_5d/src/electrical_resistivty/Release`. If not, follow the steps of the next item `Error in "Parameters.cpp"`. Then try to compile again.

If the error still appears check python version (here, python3.13 was used on the local computer, and 3.10 on the cluster). Depending on the version, some errors may occur with python thread state. Check the file "PERFORM.d", the version should be the same as the one used in your Anaconda environment or the version installed on your computer. To check python version, use the command: `python3 -V` (ubuntu).

### Error in "Parameters.cpp"

If the compilation fails and an error occurs on the parameters, comment `#include <filesystem>` (write `//` in front of the line) at the beginning of `InputOut/Parameters.cpp`.

### Alpha delayed set to Infinity

It may happens to get alpha_delayed equal to -Infinity (TT3t=-Infinity in `alpha_calculation.f90`). This occurs when the likelihood delayed is equal or higher than the likelihood proposed, i.e. the second try leads to a worse likelihood (delayed) than the first try (proposed). Setting alpha to -Infinty is a way to ensure the rejection of the model.

### keyError (still unsolved)

A keyError may occur when loading the C++ results in the post-processing python functions using the libraries PERFORM.so and PyBind11. It means that a given data is not found in the dictionnary. This error seems to be linked with eigen structures.
