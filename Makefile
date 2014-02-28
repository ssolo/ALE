#Please change to reflect your Bio++ and Boost installation:
bpp_DIR= /Users/ssolo/newest_bpp/
boost_DIR=/usr/local/
#works:
# for Bio++ v2.0.3 (maybe works for 2.0.x)
# for Boost v1.51 (should work for other versions) 

CC=g++ -pipe
mCC=mpic++

#works: 
#on Ubuntu 
#g++ -v
#gcc version 4.6.3 (Ubuntu/Linaro 4.6.3-1ubuntu5)  
#
#and OS X 10.7
#g++ -v
#gcc version 4.8.0 20120930 (experimental) (GCC) 
#ALSO
#llvm-g++ -v 
#gcc version 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)
#from http://hpc.sourceforge.net
#
#Bio++ did not complie/link with clang.

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall -std=gnu++11
OMP_FLAGS = -fopenmp
DEV_FLAGS =  -g -Wall -std=gnu++11 

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
  export OSTYPE
endif


ifeq ($(OSTYPE),darwin)
	#pre Mountaion Lion compatibility
	CC=llvm-g++ -mmacosx-version-min=10.6 
	CC=g++
	bpp_libs = $(bpp_DIR)lib/libbpp-core.a $(bpp_DIR)lib/libbpp-phyl.a $(bpp_DIR)lib/libbpp-seq.a 
	OSSTRING=OSX
else
	bpp_libs = -lbpp-core -lbpp-seq -lbpp-phyl
	OSSTRING=LINUX

endif

STATIC= libexODT.a
STATIC_OMP= libexODT_omp.a

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) 
LINK = $(DYNAMIC) 
INCLUDE = -I$(boost_DIR)include -I$(bpp_DIR)include  
MPI_INCLUDE=/usr/local/include/openmpi #-I/usr/lib/openmpi/include/ 
MPI_LINK= -lboost_mpi -lboost_serialization 

ALE.o: ALE.h ALE.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o ALE.o ALE.cpp

ALE_util.o: ALE.h ALE_util.h ALE_util.cpp Makefile
	 $(CC) $(FLAGS) $(INCLUDE)  -c -o ALE_util.o ALE_util.cpp

exODT.o: ALE.h exODT.h exODT.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o exODT.o exODT.cpp

model.o: ALE.h exODT.h model.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o model.o model.cpp

model_omp.o: ALE.h exODT.h model_omp.cpp Makefile 
	$(CC) $(FLAGS) $(OMP_FLAGS) $(INCLUDE) -c -o model_omp.o model_omp.cpp

traceback.o: ALE.h exODT.h traceback.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback.o traceback.cpp

sample.o: ALE.h exODT.h sample.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample.o sample.cpp	


libexODT.a: ALE.o ALE_util.o exODT.o model.o traceback.o sample.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model.o traceback.o  ALE_util.o  sample.o

libexODT_omp.a: ALE.o ALE_util.o exODT.o model_omp.o traceback.o sample.o Makefile
	ar rcs libexODT_omp.a ALE.o exODT.o model_omp.o traceback.o  ALE_util.o  sample.o

ALEml:	libexODT.a ALEml.cpp Makefile
	$(CC) ALEml.cpp -o ALEml $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEml_omp:	libexODT_omp.a ALEml.cpp Makefile
	$(CC) ALEml.cpp -o ALEml_omp $(FLAGS) $(OMP_FLAGS) $(INCLUDE) $(STATIC_OMP) $(LINK)

ALEsample:	libexODT.a ALEsample.cpp Makefile
	$(CC) ALEsample.cpp -o ALEsample $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEsample_omp:	libexODT_omp.a ALEsample.cpp Makefile
	$(CC) ALEsample.cpp -o ALEsample_omp $(FLAGS) $(OMP_FLAGS) $(INCLUDE) $(STATIC_OMP) $(LINK)

ALEobserve:	libexODT.a ALEobserve.cpp Makefile
	$(CC) ALEobserve.cpp -o ALEobserve $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEcount:	libexODT.a ALEcount.cpp Makefile
	$(CC) ALEcount.cpp -o ALEcount $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

test_omp:	libexODT_omp.a test.cpp Makefile
	$(CC) test.cpp -o test_omp $(FLAGS)  $(OMP_FLAGS) $(INCLUDE) $(STATIC_OMP) $(LINK)

test:	libexODT.a test.cpp Makefile
	$(CC) test.cpp -o test $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALE_tutorial:	libexODT.a ALE_tutorial.cpp Makefile
	$(CC) ALE_tutorial.cpp -o ALE_tutorial $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


test_simpleML:	libexODT.a test_simpleML.cpp Makefile
	$(CC) test_simpleML.cpp -o test_simpleML $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

simulation: simulation.cpp Makefile
	$(CC) simulation.cpp -o simulation $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


bin:  ALEobserve ALEml ALEml_omp ALEsample ALEsample_omp
	mv ALEobserve binaries/ALEobserve_$(OSSTRING)
	mv ALEml binaries/ALEml_$(OSSTRING)
	mv ALEsample binaries/ALEsample_$(OSSTRING)
	mv ALEml_omp binaries/ALEml_omp_$(OSSTRING)
	mv ALEsample_omp binaries/ALEsample_omp_$(OSSTRING)

mpi_tree.o: ALE.h exODT.h mpi_tree.h mpi_tree.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o mpi_tree.o mpi_tree.cpp

mpi_ml:	libexODT.a mpi_tree.o mpi_ml.cpp Makefile
	$(mCC) mpi_ml.cpp -o mpi_ml $(FLAGS) $(INCLUDE) $(STATIC) mpi_tree.o $(LINK) $(MPI_LINK) 
