# README #

### Description ###

This repository contains part of the source code necessary to perform non-linear tomography with the rj-McMC algorithm. 

It must be downloaded together with the companion repositories:

+ [TDTOMO\_ForTERT](https://bitbucket.org/glterica/tdtomo_fortert) for electrical resistivity tomography
+ [TDTOMO\_ForTSWT](https://bitbucket.org/glterica/tdtomo_fortswt) for seismic traveltime tomography

### How do I get set up? ###

You can download the latest version of the complete source code for electrical resistivity or seismic traveltime tomography as pre-packaged zip files by clicking [here for reistivity tomography (TERT)](http://glterica.pythonanywhere.com/software/download/tert), and [here for seismic tomography (TSWT)](http://glterica.pythonanywhere.com/software/download/tswt).

Alternatively, if you wish to obtain the source codes by downloading/cloning the repositories, follow these steps:

+ create a directory on your machine which will contain both this repository and either companion repository
+ download/clone this repository
+ download/clone either companion repository [TDTOMO\_ForTERT](https://bitbucket.org/glterica/tdtomo_fortert) or [TDTOMO\_ForTSWT](https://bitbucket.org/glterica/tdtomo_fortswt), which contain the rest of the necessary source codes
+ rename the folder of this repository to "COMMON" to ensure all links in the companion repositories work properly

Detailed instructions on how to compile and execute the codes are provided in the ```user_guide.pdf``` files in the companion repositories. Example datasets for each case are included in the ```examples``` folder of each companion repository.

### Suggested reading ###

Visit the companion repositories for a list of papers describing the theory and applications of non-linear transdimensional tomography.

### Contact ###

If you would like to collaborate or contribute to this code, for questions and to report errors/bugs in the code or examples, you can contact:

+ Erica Galetti: [erica.galetti@ed.ac.uk](mailto:erica.galetti@ed.ac.uk)
+ Andrew Curtis: [Andrew.Curtis@ed.ac.uk](mailto:Andrew.Curtis@ed.ac.uk)
+ Both: [td.tomo.codes@gmail.com](mailto:td.tomo.codes@gmail.com)

***

### To-do list ###

+ Add option to use a custom random seed
+ Add option to save single samples in ```procsamples```
+ Add option to average over neighbouring pixels for maximum-a-posteriori model in ```procsamples```

### Changelog ###

--

***

### License ###

This program is licensed under [GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0-standalone.html). A copy of the license is available within this repository by clicking [here](https://bitbucket.org/glterica/tdtomo_common/src/master/LICENSE.md).

### Credits ###

The programs in this package make use of the following libraries:

+ [Bessel programs](http://jean-pierre.moreau.pagesperso-orange.fr/f_bessel.html)
+ [BLAS](http://www.netlib.org/blas/)
+ [LAPACK](http://www.netlib.org/lapack/)
+ [mt19937](https://jblevins.org/mirror/amiller/mt19937.f90)
+ [ORDERPACK 2.0](http://www.fortran-2000.com/rank/index.html)
+ [SuperLU](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/)
