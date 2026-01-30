/*
 * FourierModelFunctions.h
 *
 *  Created on: Apr 25, 2013
 *      Author: delphineroubinet
 */

#ifndef FOURIERMODELFUNCTIONS_H_
#define FOURIERMODELFUNCTIONS_H_

#include "../InputOutput/Parameters.h"
#include "../DomainDefinition/Domain.h"

#include<Eigen/IterativeLinearSolvers>
//#include </home/leal/Desktop/Doctorat/Scripts/MontpInv_Python/transdimension/refonte_matlab_python/montp_inv/Library/Eigen/IterativeLinearSolvers>

void FourierPotentialEPM(Parameters,Domain,bool);
//void FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >);
std::vector<Eigen::SparseMatrix<double> > FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >); // vector version
//std::map<int, Eigen::SparseMatrix<double> > FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >); // map version
//std::map<int, Eigen::SparseMatrix<double>, std::less<int>,Eigen::aligned_allocator<std::pair<const int,
//    Eigen::SparseMatrix<double>> > > FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool,std::map<int,std::pair<double,double> >);
void FourierPotentialEPMSourceTerm(Parameters,Domain,bool,bool);
//void FourierPotentialDDPSourceTerm(Parameters,Domain,bool,std::map<int,std::pair<double,double> >);
std::vector<Eigen::SparseMatrix<double> > FourierPotentialDDPSourceTerm(Parameters,Domain,bool,std::map<int,std::pair<double,double> >); // vector version
//std::map<int, Eigen::SparseMatrix<double> > FourierPotentialDDPSourceTerm(Parameters,Domain,bool,std::map<int,std::pair<double,double> >); // map version
//std::map<int, Eigen::SparseMatrix<double>, std::less<int>,Eigen::aligned_allocator<std::pair<const int,
//    Eigen::SparseMatrix<double>> > > FourierPotentialDDPSourceTerm(Parameters,Domain,bool,std::map<int,std::pair<double,double> >);
void FourierPotentialDDPSourceTerm(Parameters,Domain,bool);
void FourierPotentialDDP(Parameters,Domain,bool,bool bc_on_border=false);

std::vector<double> ReadSourceTermsCoordX(Parameters);
std::map<int,std::pair<double,double> > ReadSourceTermsCoordinates(Parameters);


#endif /* FOURIERMODELFUNCTIONS_H_ */
