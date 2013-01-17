#Please change to reflect your Bio++ and Boost installation:
bpp_DIR= /usr/local/
boost_DIR=/usr/local/
#works:
# for Bio++ v2.0.3 (maybe works for 2.0.x)
# for Boost v1.51 (should work for other versions) 

CC=g++ -pipe

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

FLAGS = -O3  -fmerge-all-constants -funroll-loops -DNDEBUG -Wall
#for debuging..
DEV_FLAGS =  -g -Wall -lprofiler

ifndef OSTYPE
  OSTYPE = $(shell uname -s|awk '{print tolower($$0)}')
  export OSTYPE
endif

ifeq ($(OSTYPE),darwin)
	#CC=llvm-g++
	bpp_libs = $(bpp_DIR)lib/libbpp-core.a $(bpp_DIR)lib/libbpp-phyl.a $(bpp_DIR)lib/libbpp-seq.a 
	STATIC = libexODT.a 
else
	bpp_libs = -lbpp-core -lbpp-seq -lbpp-phyl
	STATIC= libexODT.a	
endif

DYNAMIC =  -L. -L$(bpp_DIR)lib $(bpp_libs) 
LINK = $(DYNAMIC) 
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

ALEobserve:	libexODT.a ALEobserve.cpp Makefile
	$(CC) ALEobserve.cpp -o ALEobserve $(FLAGS) $(INCLUDE) $(STATIC) $(LINK)

