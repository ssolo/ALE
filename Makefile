
CC=g++ -pipe
CC=g++

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall
DEV_FLAGS =  -g -Wall -lprofiler

boost_DIR=/usr/local/
#bpp_DIR= $(HOME)/newest_bpp/
bpp_DIR= /usr/local/
boost_libs= -lboost_mpi -lboost_serialization 


MPI_INCLUDE=-I/usr/lib/openmpi/include/ 

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
  export OSTYPE
endif

ifeq ($(OSTYPE),darwin)
	bpp_libs = $(bpp_DIR)lib/libbpp-core.a $(bpp_DIR)lib/libbpp-phyl.a $(bpp_DIR)lib/libbpp-seq.a 
	boost_libs = 
	STATIC = libexODT.a $(boost_DIR)lib/libboost_serialization.a  $(boost_DIR)lib/libboost_mpi.a	
else
	bpp_libs = -lbpp-core -lbpp-seq -lbpp-phyl
	boost_libs = -lboost_mpi -lboost_serialization 
	#boost_DIR= $(HOME)/
	STATIC= libexODT.a	
endif

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) $(boost_libs) 
LINK = $(DYNAMIC) #-lexODT
INCLUDE = -I$(boost_DIR)include -I$(bpp_DIR)include  $(MPI_INCLUDE)


ALE.o: ALE.h ALE.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o ALE.o ALE.cpp

ALE_util.o: ALE.h ALE_util.h ALE_util.cpp Makefile
	 $(CC) $(FLAGS) $(INCLUDE)  -c -o ALE_util.o ALE_util.cpp

exODT.o: ALE.h exODT.h exODT.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o exODT.o exODT.cpp

model.o: ALE.h exODT.h model.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o model.o model.cpp

traceback.o: ALE.h exODT.h traceback.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o traceback.o traceback.cpp

sample.o: ALE.h exODT.h sample.cpp Makefile 
	$(CC) $(FLAGS) $(INCLUDE)  -c -o sample.o sample.cpp	


libexODT.a: ALE.o ALE_util.o exODT.o model.o traceback.o sample.o Makefile
	ar rcs libexODT.a ALE.o exODT.o model.o traceback.o  ALE_util.o  sample.o

ALEml:	libexODT.a ALEml.cpp Makefile
	$(CC) ALEml.cpp -o ALEml $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEsample:	libexODT.a ALEsample.cpp Makefile
	$(CC) ALEsample.cpp -o ALEsample $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

ALEobsorve:	libexODT.a ALEobsorve.cpp Makefile
	$(CC) ALEobsorve.cpp -o ALEobsorve $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

