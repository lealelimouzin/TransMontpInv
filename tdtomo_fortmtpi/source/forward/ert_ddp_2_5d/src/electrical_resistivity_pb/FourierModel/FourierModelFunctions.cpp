/*
 * FourierModelFunctions.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#include "FourierModelFunctions.h"
#include "../EPMModel/EPMDefinition.h"
#include "../InputOutput/Results.h"
#include "../DDPModel/DDPModel.h"
#include "../DDPModel/SystLinearFunctions.h"
#include "../DFNModel/DFNGeneration.h"
#include <omp.h>

using namespace std;

typedef Eigen::Triplet<double> T;

std::map<int, Eigen::SparseMatrix<double> > FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border,
                                   map<int, pair<double, double> > source_term_positions) { // 25-12 map version
/*std::vector<Eigen::SparseMatrix<double> > FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border,
                                   map<int, pair<double, double> > source_term_positions) {*/  // 25-12 vector version
//std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
//         Eigen::aligned_allocator<std::pair<const int, Eigen::SparseMatrix<double>> > >
//         FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border,
//                                   map<int, pair<double, double> > source_term_positions) {
/*void FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border,
                                   map<int, pair<double, double> > source_term_positions) {*/
    // Parameters
    std::map<int, Eigen::SparseMatrix<double> > res_matrices; // 25-12 map version
//    std::vector<Eigen::SparseMatrix<double> > res_matrices(source_term_positions.size()); // 25-12 vector version
    //cout << "EPM: res_matrices size 1 = " << res_matrices.size() << endl;
//std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
    //     Eigen::aligned_allocator<std::pair<const int, Eigen::SparseMatrix<double>> > > res_matrices; // modif lea 12/05/25

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    bool fourier = true;
    int elec_nb;
    double source_position_x, source_position_y, source_value = 1.0;
    // Initial domain and system
    clock_t t_begin=clock(),t_end;
    EPMSystem epd_model_init(param), epd_model;
    epd_model_init.DetermineBasicLinearSystem();
    t_end=clock();
	//cout << "initial domain and system Time = " << (t_end-t_begin)/(double)CLOCKS_PER_SEC << endl;
    // Loop on the different electrode positions <-> different experiments
    for (map<int, pair<double, double> >::iterator it1 = source_term_positions.begin();
         it1 != source_term_positions.end(); it1++) {
        // Definition of the considered electrode position
        elec_nb = it1->first;
        //cout << "Experiment " << elec_nb << endl;
        source_position_x = it1->second.first, source_position_y = it1->second.second;
        if (param.simu_param.source_terms.size() == 0) {
            param.simu_param.source_terms.push_back(SourceTerm(source_position_x, source_position_y, source_value));
        } else { param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y, source_value); };
        // Each experiment is run for several wave values
        //map<double,ublas_matrix> Potential_Matrices;
        map<double, ublas_vector> Potential_Vector;
        for (vector<double>::iterator it2 = w_vect.begin(); it2 != w_vect.end(); it2++) {
            param.simu_param.EPM_param.w_current = *it2;
            // update the EPM model and linear system
            epd_model = epd_model_init;
            epd_model.w = param.simu_param.EPM_param.w_current;
            // Fourier terms
            epd_model.UpdateLinearSystemFourier();
            // Define the boundary conditions and source
            SourceTermsEPM source_terms_epm_;
            epd_model.bc_def_epm = BoundaryConditionsDefinition(domain, epd_model.Nx, epd_model.Ny,
                                                                param.simu_param.bc_map, param.simu_param.source_terms,
                                                                source_terms_epm_, fourier, epd_model.w);
            epd_model.source_terms_epm = source_terms_epm_;
            // Update the linear system with the boundary conditions and source terms
            epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);
            //Potential_Matrices[*it1]=epd_model.PotentialComputation();
            clock_t t_begin=clock(),t_end;
            Potential_Vector[*it2] = LinearSystemSolving(epd_model.matrix_syst, epd_model.vector_syst);
            t_end=clock();
	       // cout << "system lineaire par elec par w Time = " << (t_end-t_begin)/(double)CLOCKS_PER_SEC << endl;
        }
        // Write the results
        /*stringstream ss;
        string str_cpt;
        ss << elec_nb;
        str_cpt = ss.str();
        string output_file = "Results" + str_cpt;
        output_file = output_file + ".txt";
        //Results results(param.files_param.output_path+output_file,Potential_Matrices);
        Results results(param.files_param.output_path + output_file, Potential_Vector, epd_model.Nx, epd_model.Ny);*/

        // Store results in the map of eigen matrices to be used then in the matlab or python post-processing
        // (required transposition and translation of the matrix which contains the electrical potential distribution
        // defined from the equivalent vector)
        //Eigen::SparseMatrix<double> A(w_vect.size()*epd_model.Ny, epd_model.Nx);
        clock_t t_begin=clock(),t_end;
        Eigen::SparseMatrix<double> A(w_vect.size()*epd_model.Nx, epd_model.Ny);
        ublas_vector current_vector;
        std::vector<T> coefficients;
        int cpt_w=0;
        for (map<double,ublas_vector>::iterator it=Potential_Vector.begin();it!=Potential_Vector.end();it++){
                current_vector=it->second;
                for (int i=0;i<epd_model.Nx;i++){
                    for (int j=0;j<epd_model.Ny;j++){
                            //coefficients.push_back(T(j+epd_model.Ny*cpt_w,epd_model.Nx-1-i,current_vector(return_index(i,j,epd_model.Ny))));
                            //coefficients.push_back(T(j+epd_model.Ny*cpt_w,(epd_model.Nx-1)-i,current_vector(return_index(i,j,epd_model.Ny))));
                            //coefficients.push_back(T(i+epd_model.Nx*cpt_w,j,current_vector(return_index(i,j,epd_model.Ny))));
                            coefficients.push_back(T(j+epd_model.Ny*cpt_w,i,current_vector(return_index(i,(epd_model.Ny-1)-j,epd_model.Ny))));
                            //cout << i << " " << j << " " << current_vector(return_index(i,(epd_model.Ny-1)-j,epd_model.Ny)) << endl;
                            /*A.setFromTriplets(coefficients.begin(), coefficients.end());
                            res_matrices[elec_nb]=A;*/
                        }
                }
                cpt_w++;
        t_end=clock();
	    //cout << "remplissage eigen par electrode Time = " << (t_end-t_begin)/(double)CLOCKS_PER_SEC << endl;
        }
        A.setFromTriplets(coefficients.begin(), coefficients.end());
        //cout << "test1" << endl;
        res_matrices[elec_nb-1]=A;
        //cout << "test2" << endl;
    }
    //cout << "EPM: res_matrices size 2 = " << res_matrices.size() << endl;
    return res_matrices;
}

void FourierPotentialEPMSourceTerm(Parameters param, Domain domain, bool singularity, bool bc_on_border) {

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    if (param.simu_param.source_terms.size() != 1) {
        cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
    }

    double source_position_x = param.simu_param.source_terms[0].x_position;
    double source_position_y = param.simu_param.source_terms[0].y_position;
    double source_value = param.simu_param.source_terms[0].source_value;
    double Lsimu = domain.domain_size_x();
    double source_spacing = param.simu_param.source_terms[0].source_spacing;
    int cpt = 0;
    bool continue_loop = true, fourier = true;

    EPMSystem epd_model_init(param), epd_model;
    epd_model_init.DetermineBasicLinearSystem();

    while (continue_loop) {
        cpt++;
        //cout << "Experiment " << cpt << endl;

        param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y, source_value);

        map<double, ublas_matrix> Potential_Matrices;
        for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
            //cout << "w=" << *it1 << endl;
            param.simu_param.EPM_param.w_current = *it1;

            // update the EPM model and linear system
            epd_model = epd_model_init;
            epd_model.w = param.simu_param.EPM_param.w_current;
            // Fourier terms
            epd_model.UpdateLinearSystemFourier();
            // Define the boundary conditions and source
            SourceTermsEPM source_terms_epm_;
            epd_model.bc_def_epm = BoundaryConditionsDefinition(domain, epd_model.Nx, epd_model.Ny,
                                                                param.simu_param.bc_map, param.simu_param.source_terms,
                                                                source_terms_epm_, fourier, epd_model.w);
            epd_model.source_terms_epm = source_terms_epm_;
            // Update the linear system with the boundary conditions and source terms
            epd_model.UpdateLinearSystemBoundaryConditionsAndSourceTerms(fourier);

            Potential_Matrices[*it1] = epd_model.PotentialComputation();
        }

        stringstream ss;
        string str_cpt;
        ss << cpt;
        str_cpt = ss.str();

        string output_file = "EPM_Xu" + str_cpt;
        output_file = output_file + ".txt";

        Results results(param.files_param.output_path + output_file, Potential_Matrices);

        source_position_x += source_spacing;
        if (source_position_x < Lsimu && source_spacing != -1) { continue_loop = true; }
        else { continue_loop = false; }

    }
}

void FourierPotentialEPM(Parameters param, Domain domain, bool singularity) {

    vector<double> w_vect = param.simu_param.EPM_param.w_vect;
    map<double, ublas_matrix> Potential_Matrices;
    for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {

        param.simu_param.EPM_param.w_current = *it1;
        EPMSystem epd_model(param, domain, true, singularity);
        Potential_Matrices[*it1] = epd_model.PotentialComputation();
    }

    Results results(param.files_param.output_path + param.files_param.output_file, Potential_Matrices);
}

void FourierPotentialDDP(Parameters param, Domain domain, bool singularity, bool bc_on_border) {
    // 1.1. Variables
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    map<double, ublas_matrix> Potential_Matrices;
    map<double, ublas_vector> Potential_Fractures;
    ublas_vector DFN_potential;
    ublas_matrix matrix_potential;
    // 1.2. Fracture network generation
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    bool fourier = true;
    clock_t t1 = clock();
    DDPModel ddp_model(param, domain, net_mesh, fourier, singularity);
    clock_t t2 = clock();
    cout << "time for generating the fracture network = t2-t1" << endl;
    // 2. Simulation for each wave number
    for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
        // 2.1. Update the current wave number
        param.simu_param.DDP_param.w_current = *it1;
        // 2.2. Define the boundary conditions for ddp model (depend on w)
        ddp_model.DDPModelFourierBC(param, domain, net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map);
        // 2.3. Compute electrical potential
        cout << "w=" << ddp_model.current_w << endl;
        ddp_model.ComputePotentialFourier(DFN_potential, matrix_potential, bc_on_border);
        // 2.4. Store the computed values
        Potential_Matrices[*it1] = matrix_potential;
        Potential_Fractures[*it1] = DFN_potential;
    }
    Results results(param.files_param.output_path + param.files_param.output_file, Potential_Matrices);
    Results results1(param.files_param.output_path + param.files_param.output_file1, Potential_Fractures);
}

std::map<int, Eigen::SparseMatrix<double> > FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity,
                                   map<int, pair<double, double> > source_term_positions) {
/*std::vector<Eigen::SparseMatrix<double> > FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity,
                                   map<int, pair<double, double> > source_term_positions) {*/ // 25-12 vector version
/*std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
         Eigen::aligned_allocator<std::pair<const int,
         Eigen::SparseMatrix<double>> > > FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity,
         map<int, pair<double, double> > source_term_positions) {*/
/*void FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity,
                                  map<int, pair<double, double> > source_term_positions) {*/
    // 1. Global variables and fracture network generation
    // 1.1. Variables
    std::map<int, Eigen::SparseMatrix<double> > res_matrices; // 25-12 map version
	//std::vector<Eigen::SparseMatrix<double> > res_matrices(source_term_positions.size()); // 25-12 vector version
	//std::vector<int> res_index_elec;
    //std::map<int, Eigen::SparseMatrix<double>, std::less<int>,
    //     Eigen::aligned_allocator<std::pair<const int, Eigen::SparseMatrix<double>> > > res_matrices; // modif lea 12/05/25
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    // 1.2 Fracture network generation
    clock_t t1 = clock();
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    clock_t t2 = clock();
    //cout << "Computation time for generating the fracture network = " << (t2-t1)/CLOCKS_PER_SEC << " seconds" << endl;
    bool fourier = true;
    // 1.3. Variables for the source term positions
    double source_value = 1.0, source_position_x, source_position_y;
    int nb_elec;
    // 2. Simulations
    // 2.1. Basic linear system (defined the coefficients depending on the porous domain properties and which do not depend on the Fourier variable, boundary conditions, and source terms)
    DDPModel ddp_model_init(param, domain, net_mesh, fourier, singularity);
    ddp_model_init.DetermineBasicLinearSystemDDP();
    t1 = clock();
    // 2.2. Loop for the simulations
    //cout << "omp_get_max_threads() = " << omp_get_max_threads() << endl;
    map<int, pair<double, double> >::iterator it1;
    vector<double>::iterator it2;
    map<double, ublas_matrix> Potential_Matrices;
    map<double, ublas_vector> Potential_Fractures;
    ublas_vector DFN_potential;
    ublas_matrix matrix_potential;
    string str_cpt, output_file_matrix;
    Results results;
    BoundaryConditionsEPM bc_def_epm = ddp_model_init.bc_def_epm, bc_def_epm_sing = ddp_model_init.bc_def_epm_sing;
    BoundaryConditionsDFN bc_map = ddp_model_init.bc_map;
    int Nx = ddp_model_init.Nx, Ny = ddp_model_init.Ny, nb_nodes_fract = ddp_model_init.nb_nodes_fract;
    double delta_x = ddp_model_init.delta_x, delta_y = ddp_model_init.delta_y;
    ublas_matrix matrix_syst;
    ublas_vector vector_syst;
    ublas_matrix porous_cond = ddp_model_init.porous_cond;
    SourceTermsEPM source_terms_epm = ddp_model_init.source_terms_epm, source_terms_epm_sing = ddp_model_init.source_terms_epm_sing;
    SubNetworkMap subnetwork_map = ddp_model_init.subnetwork_map;
/*#pragma omp parallel
    {
//#pragma omp single nowait
//#pragma omp single
 //      {
            for (it1 = source_term_positions.begin(); it1 != source_term_positions.end(); it1++) {
//#pragma omp task firstprivate(res_matrices, it1, it2, nb_elec, source_position_x, source_position_y, param, Potential_Matrices, Potential_Fractures, DFN_potential, matrix_potential, str_cpt, output_file_matrix, results, matrix_syst, vector_syst, bc_def_epm, bc_map, bc_def_epm_sing, source_terms_epm, source_terms_epm_sing, subnetwork_map)
#pragma omp task firstprivate(it1, it2, nb_elec, source_position_x, source_position_y, param, Potential_Matrices, Potential_Fractures, DFN_potential, matrix_potential, str_cpt, output_file_matrix, results, matrix_syst, vector_syst, bc_def_epm, bc_map, bc_def_epm_sing, source_terms_epm, source_terms_epm_sing, subnetwork_map)
                {*/
    #pragma omp parallel private(it1, it2, nb_elec, source_position_x, source_position_y, param, Potential_Matrices, Potential_Fractures, DFN_potential, matrix_potential, str_cpt, output_file_matrix, results, matrix_syst, vector_syst, bc_def_epm, bc_map, bc_def_epm_sing, source_terms_epm, source_terms_epm_sing, subnetwork_map) 
    {
        #pragma omp single nowait
        {
        for (it1 = source_term_positions.begin(); it1 != source_term_positions.end(); it1++){         
                    nb_elec = it1->first;
                    source_position_x = it1->second.first;
                    source_position_y = it1->second.second;
                    //printf("Experiment %d traité par le thread %d \n", nb_elec, omp_get_thread_num());
                    //2.1. Variables
                    if (param.simu_param.source_terms.size() == 0) {
                        param.simu_param.source_terms.push_back(
                                SourceTerm(source_position_x, source_position_y, source_value));
                    }
                    else {
                        param.simu_param.source_terms[0] = SourceTerm(source_position_x, source_position_y,
                                                                      source_value);
                    }
		    
                    map<double, ublas_matrix> Potential_Matrices;
                    map<double, ublas_vector> Potential_Fractures;
                    ublas_vector DFN_potential;
                    ublas_matrix matrix_potential;
                    // 2.2. Simulation for each wave number
                    for (it2 = w_vect.begin(); it2 != w_vect.end(); it2++) {
                        // 2.1. Update the current wave number and linear system
                        param.simu_param.DDP_param.w_current = *it2;
                        matrix_syst = ddp_model_init.matrix_syst;
                        vector_syst = ddp_model_init.vector_syst;
                        //cout << matrix_syst(0,0) << " ";

                        // 2.2. Define the boundary conditions for ddp model (depend on w)
                        DDPModelFourierBC(singularity, param, domain, net_mesh.domain, net_mesh.border_map,
                                          net_mesh.nodes_map, source_terms_epm_sing, source_terms_epm, bc_map,
                                          bc_def_epm_sing, bc_def_epm);
                        // 2.3. Update the linear system with the coefficients which depend on the Fourier variables, boundary conditions, and source terms
                        UpdateLinearSystemFourier(Nx, Ny, delta_x, delta_y, nb_nodes_fract, *it2, bc_def_epm, bc_map,
                                                  source_terms_epm, porous_cond, subnetwork_map, singularity,
                                                  matrix_syst, vector_syst);
                        //cout << matrix_syst(0,0) << " ";
                        // 2.4. Compute electrical potential
                        EvaluatePotential(nb_nodes_fract, Nx, Ny, matrix_syst, vector_syst, DFN_potential,
                                          matrix_potential);
                        // 2.5. Store the computed values
                        //cout << matrix_potential(0,0) << " ";
                        Potential_Matrices[*it2] = matrix_potential;
                        Potential_Fractures[*it2] = DFN_potential;
                        //matrix_syst.resize(0,0,false);vector_syst.clear();
                    }
                    //cout << endl;
                    // 2.3. Write the results
                    /*str_cpt = to_string(nb_elec);
                    output_file_matrix = "Matrix_Potential" + str_cpt;
                    output_file_matrix = output_file_matrix + ".txt";
                    results = Results(param.files_param.output_path + output_file_matrix, Potential_Matrices);*/

                    // 2.4. Free memory
                    DFN_potential.clear();

                    // Store results in the map of eigen matrices to be used then in the matlab or python post-processing
                    // (required transposition and translation of the matrix which contains the electrical potential distribution
                    // defined from the equivalent vector)
                    Eigen::SparseMatrix<double> A(w_vect.size()*Nx, Ny);
                    //ublas_vector current_vector;
                    ublas_matrix current_matrix;
                    std::vector<T> coefficients;
                    int cpt_w=0;
                    for (map<double,ublas_matrix>::iterator it=Potential_Matrices.begin();it!=Potential_Matrices.end();it++){
                            current_matrix=it->second;
                            for (int i=0;i<Nx;i++){
                                for (int j=0;j<Ny;j++){
                                	//coefficients.push_back(T(j+epd_model.Ny*cpt_w,i,current_vector(return_index(i,(epd_model.Ny-1)-j,epd_model.Ny))));
                                	coefficients.push_back(T(j+Ny*cpt_w,i,current_matrix(i,(Ny-1)-j)));
                                	//cout << i << " " << j << " " << current_matrix(i,(Ny-1)-j) << endl;
                                }
                            //if (int i=1){return res_matrices;}
                            }
                            cpt_w++;
                    }
                    A.setFromTriplets(coefficients.begin(), coefficients.end());
                    res_matrices[nb_elec-1]=A;
                    if (param.simu_param.debug_mode){
                        cout << "DEBUG: res_matrices calculated for simulation number " << nb_elec-1 << endl;
                        cout << res_matrices[nb_elec-1] << endl;
                    }
                    //res_matrices.push_back(A);
					//res_index_elec.push_back(nb_elec);
                }
            }
    }
    t2 = clock();
    //cout << "Computation time for simulating all the experiments = "
     //    << (t2 - t1) / CLOCKS_PER_SEC / omp_get_max_threads() << " seconds" << endl;
    //exit(0);
    if (param.simu_param.debug_mode){
        cout << "DEBUG: print the first res_matrices at the end of the Fourrier DDP function" << endl;
        cout << res_matrices[0] << endl;
    }
    return res_matrices;
}

void FourierPotentialDDPSourceTerm(Parameters param, Domain domain, bool singularity) {
    // 1. Global variables and fracture network generation
    // 1.1. Variables
    vector<double> w_vect = param.simu_param.DDP_param.w_vect;
    if (param.simu_param.source_terms.size() != 1) {
        cout << "WARNING in main (PERFORM.cpp): configuration not implemented" << endl;
    }
    // 1.2 Fracture network generation
    NetworkMeshes net_mesh(domain, param, param.simu_param.simu_option);
    /*cout << "NETWORK RESTRICTRED TO ITS BACKBONE AND UPPER VALUES" << endl;
    // Compute the backbone
    net_mesh=net_mesh.return_backbone(param,EPSILON);
    // Cut the domain
    cout << "NETWORK RESTRICTRED TO THE UPPER FRACTURES" << endl;
    DFNCutting(net_mesh,25);*/
    bool fourier = true;
    // 2. Simulation for each source position
    // 2.1. Parameters for the simulations
    // standard parameters and variables
    double source_value = param.simu_param.source_terms[0].source_value;
    double source_position_y = param.simu_param.source_terms[0].y_position;
    double Lsimu = domain.domain_size_x();
    int cpt = 0;
    // definition of the x-coordinate of the source term
    double source_spacing = param.simu_param.source_terms[0].source_spacing;
    vector<double> vect_source_pos_x, vect_values_I;
    ublas_matrix elec_positions;
    // several positions for the source terms defined in a text file (electrical resistivity survey)
    if (source_spacing == -1) { vect_source_pos_x = ReadSourceTermsCoordX(param); }
        // one position for the source terms
    else if (source_spacing == -2) {
        double source_position_x = param.simu_param.source_terms[0].x_position;
        vect_source_pos_x.push_back(source_position_x);
    } else {
        double source_position_x = param.simu_param.source_terms[0].x_position;
        vect_source_pos_x.push_back(source_position_x);
        //while (source_position_x<=Lsimu-4*source_spacing){
        while (source_position_x < Lsimu - source_spacing) {
            source_position_x += source_spacing;
            vect_source_pos_x.push_back(source_position_x);
        }
    }
    // Basic linear system (defined the coefficients depending on the porous domain properties and which do not depend on the Fourier variable, boundary conditions, and source terms)
    DDPModel ddp_model_init(param, domain, net_mesh, fourier, singularity), ddp_model;
    ddp_model_init.DetermineBasicLinearSystemDDP();
    ublas_vector potential;
    // 2.2. Loop for the simulations
    for (vector<double>::iterator it1 = vect_source_pos_x.begin(); it1 != vect_source_pos_x.end(); it1++) {
        cpt++;
        //cout << "Experiment " << cpt << " and x_coord=" << *it1 << endl;
        //2.1. Variables
        param.simu_param.source_terms[0] = SourceTerm(*it1, source_position_y, source_value);
        map<double, ublas_matrix> Potential_Matrices;
        map<double, ublas_vector> Potential_Fractures;
        ublas_vector DFN_potential;
        ublas_matrix matrix_potential;
        // 2.2. Simulation for each wave number
        for (vector<double>::iterator it1 = w_vect.begin(); it1 != w_vect.end(); it1++) {
            ddp_model = ddp_model_init;
            // 2.1. Update the current wave number
            param.simu_param.DDP_param.w_current = *it1;
            ddp_model.current_w = *it1;
            //cout << "w=" << ddp_model.current_w << endl;
            // 2.2. Define the boundary conditions for ddp model (depend on w)
            ddp_model.DDPModelFourierBC(param, domain, net_mesh.domain, net_mesh.border_map, net_mesh.nodes_map);
            // 2.3. Update the linear system with the coefficients which depend on the Fourier variables, boundary conditions, and source terms
            ddp_model.UpdateLinearSystemFourier();
            // 2.4. Compute electrical potential
            ddp_model.EvaluatePotential(DFN_potential, matrix_potential);
            // 2.4. Store the computed values
            Potential_Matrices[*it1] = matrix_potential;
            Potential_Fractures[*it1] = DFN_potential;
        }
        stringstream ss;
        string str_cpt;
        ss << cpt;
        str_cpt = ss.str();
        string output_file_matrix = "Matrix_Potential" + str_cpt;
        output_file_matrix = output_file_matrix + ".txt";
        string output_file_fract = "Fracture_Potential" + str_cpt;
        output_file_fract = output_file_fract + ".txt";
        Results results(param.files_param.output_path + output_file_matrix, Potential_Matrices);
        Results results1(param.files_param.output_path + output_file_fract, Potential_Fractures);
    }
}

// Function to read the x-coordinates of the source term
vector<double> ReadSourceTermsCoordX(Parameters param) {
    vector<double> vect_source_pos_x;
    //string file_name_coord=param.files_param.input_path+param.files_param.file_name_source_xcoord;
    string file_name_coord = param.files_param.input_path + param.files_param.file_name_source_positions;
    ifstream file(file_name_coord.c_str());
    if (file.is_open()) {
        int Nb_coord;
        double x_coord;
        file >> Nb_coord;
        for (int i = 0; i < Nb_coord; i++) {
            file >> x_coord;
            vect_source_pos_x.push_back(x_coord);
        }
    } else {
        cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
    }
    return vect_source_pos_x;
}

// Function to read the source term positions
map<int, pair<double, double> > ReadSourceTermsCoordinates(Parameters param) {
    map<int, pair<double, double> > source_positions;
    string file_name_coord = param.files_param.input_path + param.files_param.file_name_source_positions;
    ifstream file(file_name_coord.c_str());
    int elec_nb;
    double coord_x, coord_y;
    if (file.is_open()) {
        while (!file.eof()) {
            file >> elec_nb >> coord_x >> coord_y;
            source_positions[elec_nb] = make_pair(coord_x, coord_y);
        }
    } else {
        cout << "WARNING in ReadSourceTermsCoordX (FourierModelFunctions.cpp): file not found" << endl;
    }
    return source_positions;
}
