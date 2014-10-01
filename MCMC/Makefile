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

MAPSDIR = /Users/noldor/physics/nonparametric_cosmology/code_sfd/GitDirectory/ModifiedAPS/

MAPSOBJ = $(MAPSDIR)containers.o $(MAPSDIR)goto_tools.o $(MAPSDIR)kd.o $(MAPSDIR)eigen_wrapper.o \
$(MAPSDIR)chisq.o

#do not use these compilers with omp
#gg = g++ -Wno-write-strings -O3
#ff = gfortran -O3

#use these compilers if you are going to use omp
gg = /opt/local/bin/g++-mp-4.8 -Wno-write-strings -O3 -fopenmp -DUSE_OPENMP -g
ff = /opt/local/bin/gfortran-mp-4.8 -O3 -g

mcmc_extractor.o: mcmc_extractor.h mcmc_extractor.cpp
	$(gg) -c mcmc_extractor.cpp

mcmc.o: mcmc.cpp mcmc.h mcmc_extractor.o
	$(gg) -c mcmc.cpp

kde.o: kde.h kde.cpp
	$(gg) -c kde.cpp

test: mcmc_test.cpp mcmc.o
	$(gg) -o test mcmc_test.cpp \
	$(MAPSOBJ) \
	mcmc.o mcmc_extractor.o \
	$(LIBRARIES) $(INCLUDE)

scurve: mcmc.o mcmc_s_curve.cpp
	$(gg) -o scurve mcmc_s_curve.cpp \
	$(MAPSOBJ) \
	mcmc.o mcmc_extractor.o \
	$(LIBRARIES) $(INCLUDE)

scurveBatch: mcmc.o mcmc_s_curve_batch.cpp
	$(gg) -o scurveBatch mcmc_s_curve_batch.cpp \
	$(MAPSOBJ) \
	mcmc.o mcmc_extractor.o \
	$(LIBRARIES) $(INCLUDE)



scurve_extract: mcmc_extractor.o mcmc_s_curve_extractor.cpp kde.o
	$(gg) -o scurve_extract mcmc_s_curve_extractor.cpp \
	$(MAPSOBJ) \
	mcmc_extractor.o kde.o \
	$(LIBRARIES) $(INCLUDE)

all:
	make mcmc.o
	make test
	make scurve
	make scurveBatch
	make scurve_extract
	
clean:
	rm *.o test scurve scurveBatch scurve_extract
