/*
 * LinearSystem.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: delphineroubinet
 */

#include "LinearSystem.h"

#include <iostream>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace ublas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;
using namespace std;

typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::column_major, 0,boost::numeric::ublas::unbounded_array<int>, boost::numeric::ublas::unbounded_array<double> > ublas_matrix_syst;

// return the solution x of the linear system Ax = b (A = mat_syst, b = vect_syst and x = vect_sol)
ublas_vector LinearSystemSolving(ublas_matrix mat_syst, ublas_vector vect_syst){

	// 1. Variables definition
	int Ny = vect_syst.size();
	ublas::vector<double> X(Ny);

	ublas_matrix_syst mat_syst_current=mat_syst;

	// function implemented for 2D matrix
	if ((unsigned)mat_syst.size2()==vect_syst.size()){
		// 2. System solving
		umf::symbolic_type<double> Symbolic;
		umf::numeric_type<double> Numeric;
		/*umf::symbolic (mat_syst, Symbolic);
		umf::numeric (mat_syst, Symbolic, Numeric);
		umf::solve (mat_syst, X, vect_syst, Numeric);*/

		umf::symbolic (mat_syst_current, Symbolic);
		umf::numeric (mat_syst_current, Symbolic, Numeric);
		umf::solve (mat_syst_current, X, vect_syst, Numeric);

	}
	else{cout << "WARNING in LinearSystemSolving (LinearSystemUtilitaries.cpp): system not well defined" << endl;}

	return X;
}



