# Installing ALE

There are two ways to use the ALE suite of programs. The simplest approach is using Docker, the other is building it form source which uses cmake and requires installing several dependencies.

## Using Docker

### Install Docker

To use ALE using Docker, you need to install Docker first.
* Under Mac OS X, this is done with a DMG file that can be downloaded from [there](https://docs.docker.com/docker-for-mac/install/#install-and-run-docker-for-mac).
* Under Windows, you can download the executable file from [there](https://docs.docker.com/docker-for-windows/install/)
* Under Linux, you need to choose your particular flavor on the left of [this page](https://docs.docker.com/engine/installation/) and then follow the instructions.

### Test Docker 

Once Docker has been installed, you can use ALE as follows. As described in README.md, you need to run ALEobserve then ALEml or ALEmcmc_undated to get reconciled gene trees, and this is explained below.

Let's assume you have a file named "geneFamily.treelist" containing a distribution of gene trees, and a species tree named species_tree.newick.
We will launch the program from the folder containing those data, which means that we will launch the command from the folder $PWD. $PWD is the usual environmental variable that UNIX-type operating systems use to point to the Present Working Directory.
If you want to use the example data to test the software, this means your $PWD will be /absolute_path/ALE/example_data, your species tree will be S.tree, and your gene tree file will be either HBG745965_real.1.treelist.txt or HBG745965_real.2.treelist.txt.

The way the docker image is as follows: you can run all the programs of the ALE suite, and only have to precede the command you would use to launch them with "docker run -v $PWD:$PWD alesuite ".

For instance, we provide below what the typical pipeline would look like.

- The first command to use is:
```sh
docker run -v $PWD:$PWD alesuite ALEobserve $PWD/geneFamily.treelist
```

- The second command to use could be:
```sh
docker run -v $PWD:$PWD alesuite ALEml_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```

or

```sh 
docker run -v $PWD:$PWD alesuite ALEmcmc_undated $PWD/species_tree.newick $PWD/geneFamily.treelist.ale
```

for instance.

Running ALEobserve then one of those two commands will produce the output files describing scenarios of gene family evolution including events of gene transfer, duplication and loss.

## Building from source

Installation from source requires:
* cmake
* a C++ compiler (e.g. g++ or clang)
* the Bio++ libraries bpp-core, bpp-seq, and bpp-phyl
* the Boost C++ libraries (serialization and mpi)

### Installing dependencies

#### Basic dependencies

** Ubuntu 16.10 **:
```sh
sudo apt-get install git cmake gcc g++
```

** openSUSE 42.2 **:
```sh
sudo zypper install git cmake gcc gcc-c++
```

#### Boost libraries

** Ubuntu 16.10 **:                                                                                           
```sh                                                                                                         
sudo apt-get install libboost-dev libboost-serialization-dev libboost-mpi-dev
```                                                                                                           
                                                                                                              
** openSUSE 42.2 **:                                                                                          
```sh                                                                                                         
sudo zypper install boost-devel libboost_mpi1_61_0 libboost_serialization_61_0                                                                     
```

#### Bio++ libraries

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

Afterwards create the build directories:
```sh
mkdir bpp-core-build
mkdir bpp-phyl-build
mkdir bpp-seq-build
```

Finally build and install the libraries.

The default installation locations are `${CMAKE_INSTALL_PREFIX}/include` and `${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}`. `$LIB_SUFFIX` is empty by default.

To install locally (e.g. to your home directory) just replace the `CMAKE_INSTALL_PREFIX` switch with the desired path. In this case you don't need sudo rights to `make install`.

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

** Ubuntu 16.10 **:                                                                                           
```sh                                                                                                         
sudo apt-get install libbpp-core-dev libbpp-phyl-dev libbpp-seq-dev libbpp-seq-omics-dev
```                                                                                                           
                                                                                                              
** openSUSE 42.2 **:

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

We advise that you create a folder named `build` in the ALE folder:

```sh
mkdir build
cd build
```

Then run cmake:

```sh
cmake ..
make
```

Using make with the option "-j 4" paralellizes compilation across 4 processes and speeds it up. Those commands will produce executable files in the folder ALE/build/bin.

## Troubleshooting  

### Cmake doesn't find boost libraries

Provided the Boost libraries are installed this error can mean that CMake just can't find any static Boost libraries.

Delete the contents of the build directory and try to build ALE with the `BUILD_STATIC` switch set to `OFF`:

```sh
cmake .. -DBUILD_STATIC=OFF
make
```

### `undefined reference to` a Bio++ function/variable

The linkage is broken, if you use prebuilt Bio++ libraries you should delete (and purge) the packages and compile them from source.

If you built Bio++ from source you could first try to set the explicitly:
```sh
cmake .. -DCMAKE_LIBRARY_PATH=/path/to/lib -DCMAKE_INCLUDE_PATH=/path/to/include
make
```

If does not work and you have multiple Bio++ library instances installed you could try to remove every instance except one and then try to build ALE.

If still not working a fresh start usually helps: remove all Bio++ libraries and rebuild them, then try to build ALE.



