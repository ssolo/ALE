# Installing ALE

There are two ways to use the ALE suite of programs. The simplest approach is using [Docker](#using-docker), the other is [building it from source](#building-from-source) which uses cmake and requires installing several dependencies.

## Using Docker

### Install Docker

To use ALE using Docker, you need to install Docker first.
* Under Mac OS X, this is done with a DMG file that can be downloaded from [there](https://docs.docker.com/docker-for-mac/install/#install-and-run-docker-for-mac).
* Under Windows, you can download the executable file from [there](https://docs.docker.com/docker-for-windows/install/)
* Under Linux, you need to choose your particular flavor on the left of [this page](https://docs.docker.com/engine/installation/) and then follow the instructions.

### Test Docker

Once Docker has been installed, you can use ALE as follows. As described in `README.md`, you need to run `ALEobserve` then `ALEml` or `ALEmcmc_undated` to get reconciled gene trees, and this is explained below.

Let's assume you have a file named `geneFamily.treelist` containing a distribution of gene trees, and a species tree named `species_tree.newick`.
We will launch the program from the folder containing those data, which means that we will launch the command from the folder `$PWD`. `$PWD` is the usual environmental variable that UNIX-type operating systems use to point to the Present Working Directory.
If you want to use the example data to test the software, this means your `$PWD` will be `/absolute_path/ALE/example_data`, your species tree will be `S.tree`, and your gene tree file will be either `HBG745965_real.1.treelist.txt` or `HBG745965_real.2.treelist.txt`.

The way the docker image is as follows: you can run all the programs of the ALE suite, and only have to precede the command you would use to launch them with `docker run -v $PWD:$PWD  -w $PWD boussau/alesuite `.

For instance, we provide below what the typical pipeline would look like.

* The first command to use is:
```sh
docker run -v $PWD:$PWD  -w $PWD boussau/alesuite ALEobserve $PWD/geneFamily.treelist
```
* The second command to use could be:
```sh
docker run -v $PWD:$PWD  -w $PWD boussau/alesuite ALEml_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```
or
```sh
docker run -v $PWD:$PWD  -w $PWD boussau/alesuite ALEmcmc_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```
for instance.

Those commands will produce the output files describing scenarios of gene family evolution including events of gene transfer, duplication and loss.

## Building from source

Installation from source requires:
* cmake
* a C++ compiler (e.g. g++ or clang >=3.8.0)
* the Bio++ libraries bpp-core, bpp-seq, and bpp-phyl
* the Boost C++ libraries (serialization and mpi)

### Installing dependencies

#### Basic dependencies

**Ubuntu 16.10**:
```sh
sudo apt-get install git cmake gcc g++
```

**openSUSE 42.2**:
```sh
sudo zypper install git cmake gcc gcc-c++
```

**CentOS 7**:
```sh
yum install git cmake gcc gcc-c++
```

On CentOS and Fedora based systems the MPI (e.g. OpenMPI, MPICH, etc.) libraries and binaries are not included in the `$PATH` environmental variable by default. You have to use the `environment-modules` package to load them manually.

First install the package:
```sh
yum install environment-modules
```

Then **logout** from your current shell and **login** again or open a new shell.

Test that the `module` command is available and working:
```sh
module avail
```

This will output something like:
```
----------------------------- /us/share/Modules/modulefiles -----------------------------
dot         module-git  module-info modules      null      use.own
```

**Mac OS X**:
At the time of this writing, the version of clang embarked in Mac OS X does not support OpenMP, so we suggest to install gcc. Assuming [homebrew](https://brew.sh/) has been installed on the mac already:
```sh
brew install gcc
brew install git
brew install cmake
```


#### Boost libraries

**Ubuntu 16.10**:                                                                                           
```sh                                                                                                         
sudo apt-get install libboost-dev libboost-serialization-dev libboost-mpi-dev
```                                                                                                           

**openSUSE 42.2**:                                                                                          
```sh                                                                                                         
sudo zypper install boost-devel libboost_mpi1_61_0 libboost_serialization_61_0                                                                 
```

**CentOS 7**:
```sh
yum install boost-devel boost-openmpi-devel boost-serialization openmpi-devel
```

**Mac OS X**:

* download boost, a version between 1.55 and 1.63 (included): https://sourceforge.net/projects/boost/files/boost/  
* unarchive it
* install it :
```sh
cd boost_directory
./bootstrap.sh --with-libraries=mpi,serialization
./b2
sudo ./b2 install
```

#### Bio++ libraries


**Mac OS X**:
For Mac OS X, we recommend installing from source, making sure to compile with the gcc installed previously with brew, adding the following options to the cmake command lines below:`  -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/7.1.0/bin/g++-7`


##### Building from source

First create a directory, where we can build the libraries.

```sh
mkdir bpp
cd bpp
```

Then clone the git repositories:
```sh
git clone https://github.com/BioPP/bpp-core
git clone https://github.com/BioPP/bpp-seq
git clone https://github.com/BioPP/bpp-phyl
```

ALE is only compatible up to BPP version 2.4.1, to select the proper version issue the following commands:
```sh
cd bpp-core
git checkout tags/v2.4.1 -b version2.4.1
cd ../bpp-seq
git checkout tags/v2.4.1 -b version2.4.1
cd ../bpp-phyl
git checkout tags/v2.4.1 -b version2.4.1
cd ..
```

Afterwards create the build directories:
```sh
mkdir bpp-core-build
mkdir bpp-seq-build
mkdir bpp-phyl-build
```

Finally build and install the libraries.

The default installation locations are `${CMAKE_INSTALL_PREFIX}/include` and `${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}`. `$LIB_SUFFIX` is empty by default.

To install locally (e.g. to your home directory) just replace the `CMAKE_INSTALL_PREFIX` switch with the desired path. In this case you don't need sudo rights to `make install`.

Depending on the g++ version a modification in `bpp-core/src/Bpp/Graph/GlobalGraph.cpp` may be necessary: add the line `#include <limits>` at line 45.


```sh
cd bpp-core-build/
cmake ../bpp-core -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install

cd ../bpp-seq-build/
cmake ../bpp-seq -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install

cd ../bpp-phyl-build/
cmake ../bpp-phyl -DCMAKE_INSTALL_PREFIX=/usr/ -DBUILD_TESTING=FALSE
make
sudo make install
```

Don't forget to step back two directories.
```sh
cd ../..
```

##### Precompiled from repository

The precompiled libraries don't always work (see Troubleshooting), our recommendation is to build Bio++ from source.

**Ubuntu 16.10**:                                                                                           
```sh                                                                                                         
sudo apt-get install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```                                                                                                           

**openSUSE 42.2**:

For openSUSE the Bio++ libraries are in a user's repository which needs to be added to Zypper:

```
sudo zypper addrepo http://download.opensuse.org/repositories/home:/jdutheil:/Bio++2.2.0/openSUSE_Factory/home:jdutheil:Bio++2.2.0.repo
sudo zypper refresh
```

After that you can install the needed Bio++ libraries:
```
sudo zypper install libbpp-core-devel libbpp-phyl-devel libbpp-seq-devel libbpp-seq-omics-devel libbpp-phyl-omics-devel
```

### Compiling ALE

Clone the git repository:
```
git clone https://github.com/ssolo/ALE.git
```

It's recommended to create a folder named `build` in the `ALE` folder:

```sh
cd ALE
mkdir build
cd build
```


**On CentOS 7**:
Before executing `cmake` load the `mpi` module:
```
module load mpi
```
Failing to do so `cmake` will find neither any MPI library nor the `boost_mpi` library. 

Remember that in case you're going to use ALE with MPI you have to execute the same command above beforehand. 

Then run cmake:
```sh
cmake .. -DBUILD_STATIC=OFF
make
```

**Any other GNU/Linux**:
Then run cmake:

```sh
cmake ..
make
```

Using make with the option "-j 4" paralellizes compilation across 4 processes and speeds it up. Those commands will produce executable files in the folder ALE/build/bin.


**Mac OS X**:
For Mac OS X, you need to make sure to compile with the gcc installed previously with brew, adding the following options to the cmake command line, i.e. :
```
cmake .. -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/7.1.0/bin/g++-7
make
```

## Troubleshooting  

### CMake doesn't find Boost libraries

Provided that the Boost libraries are installed, this error can mean that CMake just can't find any _static_ Boost libraries.

Delete the contents of the build directory:
```
rm -rf build/*
```
and try to build ALE with the `BUILD_STATIC` switch set to `OFF`:

```sh
cmake .. -DBUILD_STATIC=OFF
make
```

### `undefined reference to` Bio++ function(s)/variable(s)

The linkage is broken, if you use prebuilt Bio++ libraries you should delete (and purge) the packages and compile them from source.

If you have built Bio++ from source you could first try to set the paths explicitly:
```sh
cmake .. -DCMAKE_LIBRARY_PATH=/path/to/lib -DCMAKE_INCLUDE_PATH=/path/to/include
make
```

In addition, you need to make sure that you are compiling ALE with the same compiler (e.g. same version of g++, or same version of clang) you used for compiling Bio++.

If this does not work and you have multiple Bio++ library instances installed you can try to remove every instance except one and then try to build ALE.

If it is still not working a fresh start usually helps: remove all Bio++ libraries and rebuild them, then try to build ALE.

### Default C++ compiler does not support one of the required libraries

After installing a compiler which has support, you can tell cmake to use it:

```sh
cmake .. -DCMAKE_C_COMPILER=/path/to/bin/c-compiler -DCMAKE_CXX_COMPILER=/path/to/bin/c++-compiler
make
```

### Could not find MPI_CXX

If you get an error like:
```
Could not find MPI_CXX (missing: MPI_CXX_LIBRARIES MPI_CXX_INCLUDE_PATH)
```

Then you either don't have any MPI implementation (e.g. OpenMPI or MPICH) installed or if you're using Red Hat-based OS, then you may have forgot to load the MPI module with:
```
module load mpi
```


### Eigen not found

Create symlink and recompile:
```
ln -s /usr/include/eigen3/Eigen /usr/include/Eigen
```
