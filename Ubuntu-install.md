# ALE installation HowTo for Ubuntu Linux systems

This guide is tested on a fresh installation of Ubuntu 16.10.

## 1 Install dependecies

### 1.1 Basic dependencies

```
sudo apt-get install git cmake
```

### 1.2 Boost
```
sudo apt-get install libboost-dev libboost-serialization-dev libboost-mpi-dev
```
### 1.3 Compile Bio++ from source

We can use the newest, development version of Bio++.

First create a directory, where we can build the libraries.

```
mkdir bpp
cd bpp
```

Then clone the git repositories:
```
git clone https://github.com/BioPP/bpp-core
git clone https://github.com/BioPP/bpp-seq
git clone https://github.com/BioPP/bpp-phyl
```

Afterwards create the build directories:
```
mkdir bpp-core-build
mkdir bpp-phyl-build
mkdir bpp-seq-build
```

Finally build and install the libraries.

The default installation locations are `${CMAKE_INSTALL_PREFIX}/include` and `${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}`.

To install locally (e.g. to your home directory) just replace the `CMAKE_INSTALL_PREFIX` switch with the desired path. In this case you do not need sudo rights to `make install`.

```
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
```
cd ../..
```

## 2 Clone ALE repository

```
git clone https://github.com/ssolo/ALE.git
```

## 3 Build ALE

Create a build directory and enter it:
```
mkdir ALE-build
cd ALE-build
```

### 3.1 With globally installed Bio++ libraries

```
cmake ../ALE
make
```

### 3.2 With locally installed Bio++ libraries
Suppose, you have following directory structure in your home directory:
```
$ ls $HOME/local
include
lib
```

In case you have more complex directory structure you can specify multiple paths separated by `;` (semicolon).

```
cmake ../ALE -DCMAKE_LIBRARY_PATH=$HOME/local
make
```

## 4 Running ALE

If everything went well, you will find the binaries in the `ALE-build/bin` directory.
