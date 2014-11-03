LAPACK_LIB =-L/Users/noldor/physics/lapack-3.1.1/ -llapack
BLAS_LIB =-L/Users/noldor/physics/BLAS/ -lblas

ARPACK_LIB = -L/Users/noldor/physics/ARPACK/ -larpack

FORTRAN_LIB = -L/opt/local/lib/ -lf95 -lgfortran

CAMB_LIB = -L/Users/noldor/physics/CAMB_110419/camb/ -lcamb
CAMB_INC = -I/Users/noldor/physics/CAMB_110419/camb/

CFITSIO_LIB = -L/Users/noldor/physics/cfitsio/ -lcfitsio

WMAP_LIB = -L/Users/noldor/physics/WMAP7likelihood/ -lwmap7 -lpthread
WMAP_INC = -I/Users/noldor/physics/WMAP7likelihood/

R_PATH = /Users/noldor/physics/lib/

LIBRARIES = $(LAPACK_LIB) $(BLAS_LIB) $(ARPACK_LIB) $(FORTRAN_LIB)

INCLUDE =

WMAP_LIBRARIES = $(LIBRARIES) $(WMAP_LIB) $(CAMB_LIB) $(CFITSIO_LIB)

WMAP_INCLUDE = $(INCLUDE) $(CAMB_INC) $(WMAP_INC)

#do not use these compilers with omp
gg = g++ -Wno-write-strings -O3
ff = gfortran -O3

#use these compilers if you are going to use omp
#gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
#ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

containers.o: containers.cpp containers.h
	$(gg) -c containers.cpp

goto_tools.o: goto_tools.h goto_tools.cpp containers.o
	$(gg) -c goto_tools.cpp containers.o $(LIBRARIES) $(INCLUDE)

test_containers: containers.o test_containers.cpp goto_tools.o
	$(gg) -o test_containers test_containers.cpp containers.o goto_tools.o

kd.o: kd.cpp kd.h goto_tools.o containers.o
	$(gg) -c kd.cpp goto_tools.o containers.o $(LIBRARIES) $(INCLUDE) \
         -Wno-deprecated

box.o: box.cpp box.h goto_tools.o kd.o
	$(gg) -c box.cpp

test_box: test_box.cpp box.o
	$(gg) -o test_box test_box.cpp goto_tools.o containers.o box.o

test_kd: test_kd.cpp kd.o
	$(gg) -o test_kd test_kd.cpp containers.o goto_tools.o kd.o

chisq.o: goto_tools.o chisq.h chisq.cpp containers.o kd.o
	$(gg) -c chisq.cpp goto_tools.o containers.o kd.o $(LIBRARIES) $(INCLUDE)

eigen_wrapper.o: goto_tools.o eigen_wrapper.cpp eigen_wrapper.h containers.o
	$(gg) -c eigen_wrapper.cpp goto_tools.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

aps_extractor.o: aps_extractor.h aps_extractor.cpp goto_tools.o
	$(gg) -c aps_extractor.cpp

aps_extract: aps_extraction_runner.cpp aps_extractor.o kd.o
	$(gg) -o aps_extract aps_extraction_runner.cpp \
	containers.o goto_tools.o aps_extractor.o kd.o

extract: mcmc_extraction_runner.cpp mcmc_extractor.o eigen_wrapper.o kd.o
	$(gg) -o extract mcmc_extraction_runner.cpp containers.o kd.o \
	goto_tools.o mcmc_extractor.o eigen_wrapper.o kde.o \
        $(LIBRARIES) $(INCLUDE)

test_eigen: test_eigen.cpp eigen_wrapper.o
	$(gg) -o test_eigen test_eigen.cpp goto_tools.o containers.o \
	eigen_wrapper.o $(LIBRARIES) $(INCLUDE) -Wno-deprecated

gaussian_process.o: gaussian_process.cpp gaussian_process.h kd.o goto_tools.o \
eigen_wrapper.o containers.o box.o
	$(gg) -c gaussian_process.cpp goto_tools.o kd.o eigen_wrapper.o \
	containers.o \
	$(LIBRARIES) $(INCLUDE)

gp_wrapper.o: gp_wrapper.h gp_wrapper.cpp chisq.o gaussian_process.o
	$(gg) -c gp_wrapper.cpp

simplex.o: gp_wrapper.o simplex.h simplex.cpp
	$(gg) -c simplex.cpp

node.o: node.cpp node.h gp_wrapper.o
	$(gg) -c node.cpp

aps.o: aps.h aps.cpp kd.o goto_tools.o eigen_wrapper.o gaussian_process.o \
chisq.o containers.o gp_wrapper.o node.o box.o
	$(gg) -c aps.cpp goto_tools.o \
	eigen_wrapper.o kd.o gaussian_process.o \
	chisq.o containers.o gp_wrapper.o node.o $(LIBRARIES) \
	$(INCLUDE) -Wno-deprecated

s_control: s_curve_control.cpp chisq.o
	$(gg) -o s_control s_curve_control.cpp containers.o goto_tools.o \
	kd.o chisq.o \
	$(LIBRARIES) $(INCLUDE)

ellipse: aps_runner_ellipses.cpp aps.o chisq.o
	$(gg) -o ellipse aps_runner_ellipses.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	chisq.o aps.o gp_wrapper.o node.o box.o \
	$(LIBRARIES) $(INCLUDE)

s_curve: aps_runner_s_curve.cpp aps.o chisq.o
	$(gg) -o s_curve aps_runner_s_curve.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	gp_wrapper.o chisq.o aps.o node.o box.o \
	$(LIBRARIES) $(INCLUDE)

coverage: aps_s_curve_coverage.cpp aps.o chisq.o
	$(gg) -o coverage aps_s_curve_coverage.cpp \
	goto_tools.o containers.o kd.o eigen_wrapper.o gaussian_process.o \
	gp_wrapper.o chisq.o aps.o node.o box.o \
	$(LIBRARIES) $(INCLUDE)

s_curve_analysis: s_curve_analyzer.cpp chisq.o aps_extractor.o
	$(gg) -o s_curve_analysis s_curve_analyzer.cpp \
	goto_tools.o containers.o kd.o aps_extractor.o \
	chisq.o $(LIBRARIES) $(INCLUDE)


s_curve_multinest_analysis: s_curve_multinest_analyzer.cpp chisq.o aps_extractor.o
	$(gg) -o s_curve_multinest_analysis s_curve_multinest_analyzer.cpp \
	goto_tools.o containers.o kd.o aps_extractor.o \
	chisq.o $(LIBRARIES) $(INCLUDE)

s_curve_mcmc_analysis: s_curve_mcmc_analyzer.cpp chisq.o aps_extractor.o
	$(gg) -o s_curve_mcmc_analysis s_curve_mcmc_analyzer.cpp \
	goto_tools.o containers.o kd.o aps_extractor.o \
	chisq.o $(LIBRARIES) $(INCLUDE)



all:
	make test_containers
	make test_kd
	make test_eigen
	make test_box
	make ellipse
	make aps_extract
	make s_control
	make s_curve
	make s_curve_analysis
	make s_curve_mcmc_analysis
	make s_curve_multinest_analysis
	make coverage
clean:
	rm *.o test_containers test_kd test_box test_eigen ellipse \
	aps_extract s_curve s_control s_curve_analysis \
	s_curve_mcmc_analysis s_curve_mcmc_analysis coverage
