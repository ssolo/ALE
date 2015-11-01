#Please change to reflect your Bio++ and Boost installation:
bpp_DIR= /home/ssolo/newest_bpp/
boost_DIR=/home/ssolo/
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

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall -std=c++0x
OMP_FLAGS = -fopenmp
DEV_FLAGS =  -g -Wall -std=c++0x 

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
	MPI_INCLUDE=/usr/local/include/openmpi #-I/usr/lib/openmpi/include/ 
	bpp_DIR= /Users/ssolo/newest_bpp/

else
	bpp_libs = -lbpp-core -lbpp-seq -lbpp-phyl
	OSSTRING=LINUX
	MPI_INCLUDE=-I/usr/lib/openmpi/include/ 	
endif

STATIC_LEGACY= libexODT_legacy.a
STATIC_SCALED= libexODT_scaled.a
STATIC_OMP= libexODT_omp.a
STATIC= libexODT.a

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) 
LINK = $(DYNAMIC) 
INCLUDE = -I$(boost_DIR)include -I$(bpp_DIR)include -I/usr/include/mpich2/
MPI_LINK= -L/home/ssolo/lib -lboost_mpi -lboost_serialization 

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

model_qvec.o: ALE.h exODT.h model_qvec.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE) -c -o model_qvec.o model_qvec.cpp

traceback_qvec.o: ALE.h exODT.h traceback_qvec.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback_qvec.o traceback_qvec.cpp

sample_qvec.o: ALE.h exODT.h sample_qvec.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample_qvec.o sample_qvec.cpp	

model_scaled.o: ALE.h exODT.h model_scaled.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE) -c -o model_scaled.o model_scaled.cpp

traceback_scaled.o: ALE.h exODT.h traceback_scaled.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback_scaled.o traceback_scaled.cpp

sample_scaled.o: ALE.h exODT.h sample_scaled.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample_scaled.o sample_scaled.cpp	

undated.o: ALE.h exODT.h undated.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o undated.o undated.cpp	

libexODT.a: ALE.o ALE_util.o exODT.o model_scaled.o traceback_scaled.o sample_scaled.o undated.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model_scaled.o traceback_scaled.o  ALE_util.o  sample_scaled.o undated.o

libexODT_scaled.a: ALE.o ALE_util.o exODT.o model_scaled.o traceback_scaled.o sample_scaled.o Makefile
	ar rcs libexODT_scaled.a ALE.o exODT.o model_scaled.o traceback_scaled.o  ALE_util.o  sample_scaled.o

libexODT_qvec.a: ALE.o ALE_util.o exODT.o model_qvec.o traceback_qvec.o sample_qvec.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model_qvec.o traceback_qvec.o  ALE_util.o  sample_qvec.o

libexODT_omp.a: ALE.o ALE_util.o exODT.o model_omp.o traceback_qvec.o sample_qvec.o Makefile
	ar rcs libexODT_omp.a ALE.o exODT.o model_omp.o traceback_qvec.o  ALE_util.o  sample_qvec.o

libexODT_legacy.a: ALE.o ALE_util.o exODT.o model.o traceback.o sample.o Makefile
	ar rcs libexODT_legacy.a ALE.o exODT.o model.o traceback.o  ALE_util.o  sample.o

ALEml:	libexODT.a ALEml.cpp Makefile
	$(CC) ALEml.cpp -o ALEml $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEml_undated:	libexODT.a ALEml_undated.cpp Makefile
	$(CC) ALEml_undated.cpp -o ALEml_undated $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEml-verbose:	libexODT.a ALEml-verbose.cpp Makefile
	$(CC) ALEml-verbose.cpp -o ALEml-verbose $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEml_scaled:	libexODT_scaled.a ALEml_scaled.cpp Makefile
	$(CC) ALEml_scaled.cpp -o ALEml_scaled $(FLAGS) $(INCLUDE) $(STATIC_SCALED) $(LINK)

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

test_qvec:	libexODT_qvec.a test.cpp Makefile
	$(CC) test.cpp -o test_qvec $(FLAGS) $(INCLUDE) $(STATIC_QVEC) $(LINK)

test:	libexODT.a test.cpp Makefile
	$(CC) test.cpp -o test $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

times:	libexODT.a times.cpp Makefile
	$(CC) times.cpp -o times $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

times_undated:	libexODT.a times_undated.cpp Makefile
	$(CC) times_undated.cpp -o times_undated $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


test_scaled:	libexODT_scaled.a test_scaled.cpp Makefile
	$(CC) test_scaled.cpp -o test_scaled $(FLAGS) $(INCLUDE) $(STATIC_SCALED)  $(LINK)

bw:	libexODT.a bw.cpp Makefile
	$(CC) bw.cpp -o bw $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

computeALEcomplexity: libexODT.a computeALEcomplexity.cpp Makefile
	$(CC) computeALEcomplexity.cpp -o computeALEcomplexity $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

mlsampler:	libexODT.a mlsampler.cpp Makefile
	$(CC) mlsampler.cpp -o mlsampler $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

mlresampler:	libexODT.a mlresampler.cpp Makefile
	$(CC) mlresampler.cpp -o mlresampler $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

mlresampler_undated:	libexODT.a mlresampler_undated.cpp Makefile
	$(CC) mlresampler_undated.cpp -o mlresampler_undated $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

mlsample:	libexODT.a mlsample.cpp Makefile
	$(CC) mlsample.cpp -o mlsample $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

mlsample_undated:	libexODT.a mlsample_undated.cpp Makefile
	$(CC) mlsample_undated.cpp -o mlsample_undated $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

wALE_ml_sample:	libexODT.a wALE_ml_sample.cpp Makefile
	$(CC) wALE_ml_sample.cpp -o wALE_ml_sample $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

wALE_ml_sample_undated:	libexODT.a wALE_ml_sample_undated.cpp Makefile
	$(CC) wALE_ml_sample_undated.cpp -o wALE_ml_sample_undated $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

summary:	libexODT.a summary.cpp Makefile
	$(CC) summary.cpp -o summary $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


ALE_tutorial:	libexODT.a ALE_tutorial.cpp Makefile
	$(CC) ALE_tutorial.cpp -o ALE_tutorial $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


test_simpleML:	libexODT.a test_simpleML.cpp Makefile
	$(CC) test_simpleML.cpp -o test_simpleML $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

simulation: simulation.cpp Makefile
	$(CC) simulation.cpp -o simulation $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

wol_tree: wol_tree.cpp Makefile
	$(CC) wol_tree.cpp -o wol_tree $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

wol_host: wol_host.cpp Makefile
	$(CC) wol_host.cpp -o wol_host $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

wol_paras: wol_paras.cpp Makefile
	$(CC) wol_paras.cpp -o wol_paras $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


simulation_cp: simulation_cp.cpp Makefile
	$(CC) simulation_cp.cpp -o simulation_cp $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

GEsim: GEsim.cpp Makefile
	$(CC) GEsim.cpp -o GEsim $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

pln: pln.cpp Makefile
	$(CC) pln.cpp -o pln $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

RF: RF.cpp Makefile
	$(CC) RF.cpp -o RF $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ls: ls.cpp Makefile
	$(CC) ls.cpp -o ls $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

simulateSpAndGeneTrees: libexODT.a simulateSpAndGeneTrees.cpp  exODT_sim.cpp Makefile
	$(CC)  simulateSpAndGeneTrees.cpp  exODT_sim.cpp -o simulateSpAndGeneTrees  $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)


bin:  ALEobserve ALEml ALEml_omp ALEsample ALEsample_omp
	mv ALEobserve binaries/ALEobserve_$(OSSTRING)
	mv ALEml binaries/ALEml_$(OSSTRING)
	mv ALEsample binaries/ALEsample_$(OSSTRING)
	mv ALEml_omp binaries/ALEml_omp_$(OSSTRING)
	mv ALEsample_omp binaries/ALEsample_omp_$(OSSTRING)

mpi_tree.o: ALE.h exODT.h mpi_tree.h mpi_tree.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o mpi_tree.o mpi_tree.cpp

mpi_ml:	libexODT.a mpi_tree.o mpi_ml.cpp Makefile
	$(mCC) mpi_ml.cpp -o mpi_ml $(FLAGS) $(INCLUDE)   $(LINK) $(MPI_LINK) mpi_tree.o $(STATIC)
mpi_ml_undated:	libexODT.a mpi_tree.o mpi_ml_undated.cpp Makefile
	$(mCC) mpi_ml_undated.cpp -o mpi_ml_undated $(FLAGS) $(INCLUDE)   $(LINK) $(MPI_LINK) mpi_tree.o $(STATIC)
mpi_ml-bw_undated:	libexODT.a mpi_tree.o mpi_ml-bw_undated.cpp Makefile
	$(mCC) mpi_ml-bw_undated.cpp -o mpi_ml-bw_undated $(FLAGS) $(INCLUDE)   $(LINK) $(MPI_LINK) mpi_tree.o $(STATIC)

mpi_S_ml:	libexODT.a mpi_tree.o mpi_S_ml.cpp Makefile
	$(mCC) mpi_S_ml.cpp -o mpi_S_ml $(FLAGS) $(INCLUDE)   $(LINK) $(MPI_LINK) mpi_tree.o $(STATIC)

#WTF 
#mpic++ mpi_ml.cpp -o mpi_ml -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall -std=gnu++11 -I/usr/include -I/home/ssolo/newest_bpp/include -I/usr/include/mpich2/   -L. -L/home/ssolo/newest_bpp/lib -lbpp-core -lbpp-seq -lbpp-phyl   -L/home/ssolo/lib -lboost_mpi -lboost_serialization  mpi_tree.o libexODT.a
#mpic++ mpi_ml.cpp -o mpi_ml -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall -std=gnu++11 -I/usr/include -I/home/ssolo/newest_bpp/include -I/usr/include/mpich2/   -L. -L/home/ssolo/newest_bpp/lib -lbpp-core -lbpp-seq -lbpp-phyl   -L/usr/lib -lboost_mpi -lboost_serialization  mpi_tree.o libexODT.a
