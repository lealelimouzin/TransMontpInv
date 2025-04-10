# TERT - Transdimensional Electrical Resistivity Tomography #

### Description ###

This repository contains part of the source code necessary to perform non-linear electrical resistivity tomography with the rj-McMC algorithm.

It must be downloaded together with the companion repository [TDTOMO\_COMMON](https://bitbucket.org/glterica/tdtomo_common).

### How do I get set up? ###

You can download the latest version of the complete source code as a pre-packaged zip file by clicking [here](http://glterica.pythonanywhere.com/software/download/tert).

Alternatively, if you wish to obtain the source codes by downloading/cloning the repository, follow these steps: 

+ create a directory on your machine which will contain both this repository and the companion repository [TDTOMO\_COMMON](https://bitbucket.org/glterica/tdtomo_common)
+ download/clone this repository
+ download/clone the companion repository [TDTOMO\_COMMON](https://bitbucket.org/glterica/tdtomo_common), which contains the rest of the necessary source codes
+ rename the folder of the companion repository to "COMMON" to ensure all links in this repository work properly

The user guide in this repository ([```user_guide.pdf```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/user_guide.pdf)) provides instructions on how to compile the Fortran codes using either the Intel or GNU compilers, as well as a description of each program.
Example datasets are included in the [```examples```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/examples) folder.

### Suggested reading ###

To learn about the theory behind non-linear and transdimensional tomography, we suggest reading the following papers:

+ Galetti E, Curtis A (2017) Transdimensional Electrical Resistivity Tomography. *Journal of Geophysical Research: Solid Earth*, submitted. The submitted copy is included in the [```papers```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/papers) folder.
+ Galetti E, Curtis A (2016) Transdimensional Monte Carlo Electrical Resistivity Tomography. *Near Surface Geoscience 2016, 22nd European Meeting of Environmental and Engineering Geophysics*. [DOI: 10.3997/2214-4609.201602017](http://dx.doi.org/10.3997/2214-4609.201602017). A copy is included in the [```papers```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/papers) folder.
+ Galetti E, Curtis A, Meles GA, Baptie B (2015) Uncertainty loops in travel-time tomography from nonlinear wave physics. *Physical Review Letters* 114: 148501. [DOI: 10.1103/PhysRevLett.114.148501](http://dx.doi.org/10.1103/PhysRevLett.114.148501). A copy is included in the [```papers```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/papers) folder.
+ Bodin T, Sambridge M, Rawlinson N, Arroucau P (2012) Transdimensional tomography with unknown data noise. *Geophysical Journal International* 189(3): 1536-1556. [DOI: 10.1111/j.1365-246X.2012.05414.x](http://dx.doi.org/10.1111/j.1365-246X.2012.05414.x)
+ Bodin T, Sambridge M (2009) Seismic tomography with the reversible jump algorithm. *Geophysical Journal International* 178(3): 1411-1436. [DOI: 10.1111/j.1365-246X.2009.04226.x](http://dx.doi.org/10.1111/j.1365-246X.2009.04226.x)

### Contact ###

If you would like to collaborate or contribute to this code, for questions and to report errors/bugs in the code or examples, you can contact:

+ Erica Galetti: [erica.galetti@ed.ac.uk](mailto:erica.galetti@ed.ac.uk)
+ Andrew Curtis: [Andrew.Curtis@ed.ac.uk](mailto:Andrew.Curtis@ed.ac.uk)
+ Both: [td.tomo.codes@gmail.com](mailto:td.tomo.codes@gmail.com)

***

### To-do list ###

+ Add 3D modeller (in progress)
+ Add option to use a custom random seed
+ Add option to save single samples in ```procsamples```
+ Clean up call to noise calculation routines from ```mksamples.f90```
+ Add option to average over neighbouring pixels for maximum-a-posteriori model in ```procsamples```

### Change log ###

--

***

### License ###

This program is licensed under [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0-standalone.html). A copy of the license is available within this repository by clicking [here](https://bitbucket.org/glterica/tdtomo_fortert/src/master/LICENSE.md).

### Credits ###

The programs in the TERT package make use of the following libraries:

+ [Bessel programs](http://jean-pierre.moreau.pagesperso-orange.fr/f_bessel.html)
+ [BLAS](http://www.netlib.org/blas/)
+ [LAPACK](http://www.netlib.org/lapack/)
+ [mt19937](https://jblevins.org/mirror/amiller/mt19937.f90)
+ [ORDERPACK 2.0](http://www.fortran-2000.com/rank/index.html)
+ [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)

The 2.5D forward modeller in [```fw2_5d```](https://bitbucket.org/glterica/tdtomo_fortert/src/master/source/forward/fw2_5d) is Erica Galetti's Fortran translation of the Matlab package [FW2\_5D by Pidlisecky & Knight (2008)](https://doi.org/10.1016/j.cageo.2008.04.001).

