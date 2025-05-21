# ----- (A) SELECT COMPILER ----- #
# F90 = ifort
# F90mpi = mpiifort
# CC = icc
F90 = gfortran
F90mpi = mpif90
CC = gcc
OPTS = -O3

# ----- (B) CREATE DIRECTORIES ----- #
bin_dir = ../../bin
gmt_dir = ${bin_dir}/gmtplot
lib_dir = ../../lib
lib_dir_local = /usr/lib64
inc = ../../../COMMON/source/tools/superlu/SRC

# ----- (C) COMPILE! ----- #
all  : dir clean tools tomo fwd main

dir	::
	mkdir -p ${bin_dir}
	mkdir -p ${gmt_dir}
	mkdir -p ${lib_dir}

tools	::
	${F90} ${OPTS} -c ../tools/debug.f90
	${F90} ${OPTS} -c ../tools/orderpack/unirnk.f90
	${F90} ${OPTS} -c ../tools/orderpack/mrgrnk.f90
	${F90} ${OPTS} -c ../tools/coordinate_to_sparse_column.f90
	${F90} ${OPTS} -c ../tools/coordinate_to_sparse_row.f90
	${F90} ${OPTS} -c ../tools/linspace.f90
	${F90} ${OPTS} -c ../tools/bessel/tbessk.f90
	${F90} ${OPTS} -c ../tools/string_conversion.f90
	ar -r ${lib_dir}/libtools.a *.o
	rm *.o
	${F90} ${OPTS} -c ../tools/lapack/*.f
	ar -r ${lib_dir}/liblapack.a *.o
	rm *.o
	${F90} ${OPTS} -c ../random/MT/mt19937.f90
	${F90} ${OPTS} -c ../random/random_generator.f90
	ar -r ${lib_dir}/librandom.a *.o
	rm *.o

fwd	::
	${F90} ${OPTS} -c ../forward/fw2_5d/modelling_2_5D.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/get_data.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/read_grid2D.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_divergence2D.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_gradient2D.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/get_fourier_coefficients.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/optimize_k_g.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/interpolate_electrodes.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_current_term.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_potential_term.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/average_conductivity_on_face.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_DS.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_A.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/calculate_data.f90
	${F90} ${OPTS} -c ../forward/fw2_5d/create_resistivity_grid2D.f90
	${F90} ${OPTS} -c ../forward/forward.f90
	ar -r ${lib_dir}/libfwd.a *.o
	rm *.o


tomo	::
	${F90} ${OPTS} -c ../input_output/acquisition.f90
	${F90} ${OPTS} -c ../inverse/voronoi_operations.f90
	${F90} ${OPTS} -c ../inverse/chains.f90
	${F90} ${OPTS} -c ../input_output/input_files.f90
	${F90} ${OPTS} -c ../input_output/check_output_files.f90
	${F90} ${OPTS} -c ../inverse/initialize_models.f90
	${F90} ${OPTS} -c ../input_output/write_last_model.f90
	${F90} ${OPTS} -c ../inverse/choose_step.f90
	${F90} ${OPTS} -c ../inverse/noise_calculation.f90
	${F90} ${OPTS} -c ../inverse/alpha_calculation.f90
	${F90} ${OPTS} -c ../inverse/chain_update.f90
	${F90} ${OPTS} -c ../input_output/write_chains_to_file.f90
	${F90} ${OPTS} -c ../input_output/inversion_parameters.f90
	${F90} ${OPTS} -c ../inverse/print_status.f90
	${F90} ${OPTS} -c ../inverse/tempering.f90
	${F90} ${OPTS} -c ../inverse/processing.f90
	${F90} ${OPTS} -c ../input_output/write_debug_info.f90
	${F90} ${OPTS} -c ../input_output/write_results_to_file.f90
	${F90} ${OPTS} -c ../input_output/write_vtx.f90
	ar -r ${lib_dir}/libtomo.a *.o
	rm *.o

main	::
	${CC}  ${OPTS} -O3 -g -DAdd_ -I${inc} -L${lib_dir_local} -c ../tools/superlu/c_fortran_dgssv.c
	${F90} ${OPTS} -o ${bin_dir}/chain2txt ../main/rjtomo/chain2txt.f90 
	${F90} ${OPTS} -o ${gmt_dir}/value_3D_post ../main/rjtomo/value_3D_post.f90 -L${lib_dir} -ltomo -ltools -lrandom 
	${F90} ${OPTS} -o ${gmt_dir}/interpolate_2D ../main/other/interpolate_2D.f90 -L${lib_dir} -ltomo -ltools -lrandom
	${F90} ${OPTS} -o ${bin_dir}/add_noise ../main/other/add_noise.f90 -L${lib_dir} -ltools -lrandom
	${F90} ${OPTS} -o ${bin_dir}/forward_2_5d ../main/rjtomo/forward_2_5d.f90 c_fortran_dgssv.o -I${inc} -L${lib_dir_local} -L${lib_dir} -lsuperlu -lfwd -ltomo -ltools  -lblas -llapack
	${F90mpi} ${OPTS} -o ${bin_dir}/mksamples ../main/rjtomo/mksamples.f90 c_fortran_dgssv.o -L${lib_dir} -lfwd -ltools -ltomo -lsuperlu -lblas -llapack -lrandom 
	${F90} ${OPTS} -o ${bin_dir}/procsamples ../main/rjtomo/procsamples.f90 c_fortran_dgssv.o -I${inc} -L${lib_dir} -L${lib_dir_local} -lfwd -ltomo -ltools -lsuperlu -lblas -llapack -lrandom
	cp ../main/rjtomo/value_3D_post.in ${gmt_dir}
	cp ../main/other/interpolate_2D.in ${gmt_dir}
	cp ../main/other/add_noise.in ${bin_dir}
	cp ../main/rjtomo/forward_2_5d.in ${bin_dir}
	cp ../main/rjtomo/mksamples.in ${bin_dir}
	cp ../main/rjtomo/procsamples.in ${bin_dir}
	cp ../running/jobsub.sh ${bin_dir}
	cp ../running/run_tert.sh ${bin_dir}

clean	::
	rm -f *.mod ${lib_dir}/libfwd.a ${lib_dir}/libtools.a ${lib_dir}/libtomo.a ${lib_dir}/liblapack.a ${lib_dir}/librandom.a 