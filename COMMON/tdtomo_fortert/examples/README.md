# README #

This folder contains two example datasets that can be used to test the correctness of your compiled programs.

Each example directory is structured as follows:

```
<example>
|-- data
    |-- data.dat
    |-- electrodes.dat
|-- inversion
    |-- forward.in
    |-- mksamples.in
    |-- procsamples.in
|-- results
    |-- plots
        |-- ...
    |-- r2
        |-- ...
    |-- plotting.sh
|-- <synthetic_model>.vtx
|-- <synthetic_model>_app.vtx
```

where:

+ ```data``` contains data files for the inversion
+ ```inversion``` contains the parameter files for inversion and processing
+ ```results/plots``` contains plots from the TERT inversion results 
+ ```results/r2``` contains results from the same dataset inverted using [R2](http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2.htm)
+ ```plotting.sh``` contains the plotting parameters used to obtain the plots in ```results/plots```
+ ```<synthetic_model>.vtx``` is the synthetic model used to generate the dataset
+ ```<synthetic_model>_app.vtx``` is the apparent resistivity pseudosection corresponding to the synthetic model

The inversion for each example was carried out by running the following:

```
$ mpiexec -np 6 mksamples 1
$ mpiexec -np 6 mksamples 2
$ mpiexec -np 6 mksamples 3
$ mpiexec -np 6 mksamples 4
$ mpiexec -np 6 mksamples 5
$ mpiexec -np 6 mksamples 6
$ ./procsamples
```

where ```mksamples``` and ```procsamples``` were compiled using the GNU compilers.

***

#### Note that the results from these example inversions are only meant to be used to check that your installation is running correcly, hence they were obtained by running *very short* Markov chains (only 6000 samples!) --> the results are not particularly good!! ####

#### When the code is used for real studies, longer Markov chains (order 10^6 samples) should be run in order to move past the burn-in period, reach convergence, and ensure that a large enough ensemble of samples is obtained. ####

