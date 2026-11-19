# TransMontpInv

Trandimensional inversion of electrical data in heterogeneous media.

## Description

TransMontpInv is a software to invert electrical resistivity data based on a discrete dual porosity approach. It is a combination of the already existent codes from Erica Galetti and Andrew Curtis ([TERT](https://bitbucket.org/glterica/tdtomo_fortert/src/master/)) and the forward model adapted from Caballero Sanz et al. (2017) (ERT_DDP_2.5D). The two codes have been merged and adapted to support the geometry of thin pathways such as fractures.

A user-guide is provided in the projetc package as well as a copy of the user-guides of TERT and ERT\_DDP\_2.5D.

## Project architecture

Part of the architecture of TransMontpInv is detailed in the following chart. Files that need to be updated by the user depending on their application sites have been highlighted using comments.

```text
.
├── TransMontpInv
│   ├── COMMON                                      # TERT source code
|   |── tdtomo_fortmtpi                             # Adaptation of tdtomo_fortert from TERT
│   │   ├── bin
|   |   |   ├── data
|   |   |   |   ├── Farum
|   |   |   |   └── Synthetic
|   |   |   ├── results
|   |   |   ├── data.dat                            # To be updated by the user
|   |   |   ├── electrodes.dat                      # To be updated by the user
|   |   |   └── ...
│   │   ├── examples
│   │   ├── lib
│   │   ├── papers
│   │   ├── source
|   |   |   ├── compile
│   │   │   ├── cptfiles
│   │   │   ├── forward
|   |   |   |   ├── ert_ddp2_5d
|   |   |   |   |   ├── fielddata                   # To be updated by the user
|   |   |   |   |   ├── post_processing
|   |   |   |   |   |    ├── farum_elecConfig.py    # To be updated by the user
|   |   |   |   |   |    └── ...
|   |   |   |   |   ├── results_chains
|   |   |   |   |   ├── src
|   |   |   |   |   |    ├── electrical_resistivity # C++ codes of ERT_DDP_2.5D
|   |   |   |   |   |    ├── input                  # To be updated by the user
|   |   |   |   |   |    └── ...
|   |   |   |   |   ├── fwdFunctions.py
|   |   |   |   |   ├── initialize_fracmodels.py    # To be updated by the user
|   |   |   |   |   └── main.py                     # To be updated by the user
|   |   |   |   ├── fw2_5d
|   |   |   |   └── forward.f90
│   │   │   ├── input_output
│   │   │   ├── inverse
│   │   │   ├── main
|   |   |   |   |   ├── other
|   |   |   |   |   └── rjtomo
|   |   |   |   |   |    ├── chain2txt.f90          # To be updated by the user
|   |   |   |   |   |    ├── mksamples.in           # To be updated by the user
|   |   |   |   |   |    ├── mksamples.f90
|   |   |   |   |   |    └── ...
│   │   │   ├── plotting
│   │   │   ├── random
│   │   │   ├── running
│   │   │   └── tools
├── batch_compilation
├── batch_transinv
├── LICENSE
├── README.md
├── requirements.txt
├── user_guide_ERT_DDP_2_5D.pdf
├── user_guide_TERT.pdf
├── user_guide_TransMontpInv.pdf
└── version
```

## License

TransMontpInv is under the GNU GPL license.

## Contacts

For questions, errors report, or contributions to this code, you can contact:

+ Léa Lelimouzin: lea.lelimouzin@gmail.fr