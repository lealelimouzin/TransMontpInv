//============================================================================
// Name        : PERFORM.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : C++ code for solute transport in heterogeneous fractured porous media
//============================================================================

#include "InputOutput/Parameters.h"
#include "DomainDefinition/Domain.h"
#include "DFNModel/NetworkMeshes.h"
#include "InputOutput/Results.h"
#include "Utilitaries/Storage/Structures.h"
#include "DFNModel/DFNComputation.h"
#include "DFNModel/DFNGeneration.h"
#include "Utilitaries/Visualization/VisuStructures.h"
#include "Utilitaries/Visualization/DisplayResults.h"
#include "EPMModel/EPMDefinition.h"
#include "EllipseComputation/EllipseComputation.h"
#include "DDPModel/DDPModel.h"
#include "DDPModel/DDPStudies.h"
#include "Utilitaries/Math/Operations.h"
#include "DDPModel/ExchangeFunctions.h"
#include "FourierModel/FourierModelFunctions.h"

#include <iostream>
#include <string>

#include <Eigen/IterativeLinearSolvers>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <map>

using Eigen::SparseMatrix;
using namespace std;
namespace py = pybind11;

typedef Eigen::Triplet<double> T;

std::map<int, Eigen::SparseMatrix<double>> perform() {

	clock_t t_begin=clock(),t_end;

	//cout << "PERFORM begin" << endl;

	bool singularity=false,bc_on_border=false;
	// 1. Parameter and domain definition reading
	// 1.1. Parameter reading
	Parameters param(singularity);
	// 1.2. Domain definition
	Domain domain(param.domain_param.Lx,param.domain_param.Ly);
	// 1.3. Display parameters
	double min_visu=-1,max_visu=-1;
    // 1.4. Output results
    std::map<int, Eigen::SparseMatrix<double> > res_map;

	// 2. Simulation
	// 2.1. DFN generation
	if (param.simu_param.simu_option==DFN_GENERATION){
		DFNVisu DFN_visu=DFNGeneration(param,domain);

	// Compute the backbone
		/*NetworkMeshes net_mesh1=DFN_visu.net_mesh,backbone;
		backbone=net_mesh1.return_backbone(param,EPSILON);
		DFN_visu.net_mesh=backbone;*/

		// Cut the domain
		/*NetworkMeshes net_mesh2=DFN_visu.net_mesh;
		//DFNCutting(net_mesh2,26.6);
		//DFNCutting(net_mesh2,25);
		DFN_visu.net_mesh=net_mesh2;*/

		// DFN writing
		Results results(param.files_param.output_path+param.files_param.output_file,DFN_visu.net_mesh);
		// TODO remove display function
		//DisplayResults(argc,argv,DFN_visu);
	}

	else if (param.simu_param.simu_option==FOURIER_POTENTIAL_COMPUTATION){
		// Equivalent Porous Media
		if (param.simu_param.model_option==EPM){
			// Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			// Run the simulations
			if (source_term_positions.size()>0){
			    //FourierPotentialEPMSourceTerm(param,domain,singularity,bc_on_border,source_term_positions);
			    res_map=FourierPotentialEPMSourceTerm(param,domain,singularity,bc_on_border,source_term_positions);
			 }
			else if (source_term_positions.size()==0){FourierPotentialEPM(param,domain,singularity);}
			else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}
		}
		// Discrete-Dual-Porosity model for fractured porous domains
		else if (param.simu_param.model_option==DDP){
			// Read the source term positions
			map<int,pair<double,double> > source_term_positions=ReadSourceTermsCoordinates(param);
			// Run the simulations
			if (source_term_positions.size()>0){
				//FourierPotentialDDPSourceTerm(param,domain,singularity,source_term_positions);
				res_map=FourierPotentialDDPSourceTerm(param,domain,singularity,source_term_positions);
				//cout << "after FourierPotentialDDPSourceTerm" << endl;
			}
			else if (source_term_positions.size()==0){FourierPotentialDDP(param,domain,singularity,bc_on_border);}
			else{cout << "WARNING in main (PERFORM.cpp): option not implemented" << endl;}

		}
		else{cout << "WARNING in main (PERFORM.cpp): solving method not implemented" << endl;}
	}



	// 2.2. Potential computation and visualization
	else if (param.simu_param.simu_option==POTENTIAL_COMPUTATION){
		// Discrete Fracture Network representation
		if (param.simu_param.model_option==DFN){
			//DFNVisu fracture_network=PotentialComputationDFN(param,domain);
			//DisplayResults(argc,argv,fracture_network);
		}
		// Equivalent Porous Media
		else if (param.simu_param.model_option==EPM){
			EPMSystem epd_model(param,domain);
			ublas_matrix EPM_Potential=epd_model.PotentialComputation();
			Results results(param.files_param.output_path+param.files_param.output_file,EPM_Potential);
		}
		// Equivalent Porous Media with a projected DFN
		else if (param.simu_param.model_option==EPMDFN){
			NetworkMeshes net_mesh(domain,param,param.simu_param.simu_option);
			//EPMVisu EPM_Potential=PotentialComputationEPMDFN(param,domain,net_mesh,min_visu,max_visu);
			//DisplayResults(argc,argv,EPM_Potential);
		}
		// Discrete Dual Porosity representation
		else if (param.simu_param.model_option==DDP){
			PotentialComputationDDP(param,domain,min_visu,max_visu);
			//DisplayResults(argc,argv,EPMDFN_Potential,"Matrix");	// "Matrix/Current" to display electric potential or current
		}
		else{cout << "WARNING in main (PERFORM.cpp): solving method not implemented" << endl;}
	}

	// 2.3. Permeability ellipse computation for DFN domain
	else if (param.simu_param.simu_option==PERMEABILITY_ELLIPSE){
		EllipseComputation(param,domain,param.simu_param.model_option);
	}

	else if (param.simu_param.simu_option==CONDUCTIVITY_PLOT){
		if (param.simu_param.model_option==DDP){
			ConductivityPlotDDP(param,domain);
		}
		else if (param.simu_param.model_option==EPMDFN){
			ConductivityPlotEPMDFN(param,domain);
		}
	}

	// 2.4. Study of the flow exchanged between fractures and matrix at the matrix-block scale
	else if (param.simu_param.simu_option==EXCHANGE_STUDY){
		/*int nb_div=6,nb_seed=1; // nb_div:number of times that we divide the domain and nb_seed:number of generated domain
		ExchangeStudy(param,domain,nb_div,nb_seed);*/
		//ExchangeStudy2(param,domain);
	}
	// 2.5. Comparison of the potential with DDP and EPMDFN
	else if (param.simu_param.simu_option==POTENTIAL_STUDY){
		int nb_div=param.simu_param.DDP_param.nb_div,nb_seed=param.simu_param.DDP_param.nb_seed;
		PotentialStudy(param,domain,nb_div,nb_seed);
	}
	else if (param.simu_param.simu_option=="TEST"){
		cout << "TEST" << endl;
		int size1=1e5,size2=size1;
		ublas_matrix matrix_syst(size1,size2);
		//print_matrix(matrix_syst);
	}
	else{cout << "WARNING in main (PERFORM.cpp): simulation option not implemented" << endl;}

	//cout << "PERFORM end" << endl;

	t_end=clock();
	//cout << "CPU Time = " << (t_end-t_begin)/CLOCKS_PER_SEC << endl;

    /*// Test eigen matrix
    int n = 5;
    SparseMatrix<double> A(n,n);
    std::vector<T> coefficients;

	for (int i=0;i<n;i++){
	    for (int j=0;j<n;j++){
            	coefficients.push_back(T(i,j,10.));
                        }
                }
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    std::map<int,Eigen::SparseMatrix<double>> map_matrices;
    map_matrices[0] = A;
	return map_matrices;*/
	return res_map;
}

PYBIND11_MODULE(PERFORM, m) {
    m.doc() = "pybind11 electrical resistivity"; // optional module docstring
    m.def("perform", &perform, py::return_value_policy::reference);
	//m.def("perform", &perform, py::return_value_policy::reference_internal);
	//m.def("perform", &perform, py::return_value_policy::move);
}