# 1 Install dependecies

## 1.1 Basic dependencies

```
sudo zypper install git cmake gcc gcc-c++
```

## 1.2 MPI
```
sudo zypper install openmpi-devel openmpi-devel-static lam lam-devel 
```

## 1.3 Boost
```
sudo zypper install boost-devel 
```

Following dependencies include version number, you may need to adjust them:
```
sudo zypper libboost_mpi1_61_0 libboost_serialization_61_0
```

## 1.4a Compile Bio++ from source

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

And finally build and install the libraries:
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

If you just copy-paste all the commands, don't forget to step back two directories.
```
cd ../..
```

## 1.4b (You can leave out this section) Prebuilt Bio++

There is an openSUSE repository `jdutheil` with Bio++ 2.2.0 libraries prebuilt, but these seems to be wrongly built, because we get 
`undefined reference to` errors for functions of Bio++ library itself.

You may give it a try, but our recommendation is to build the Bio++ libraries yourself.

```
sudo zypper addrepo http://download.opensuse.org/repositories/home:/jdutheil:/Bio++2.2.0/openSUSE_Factory/home:jdutheil:Bio++2.2.0.repo
sudo zypper refresh
```

After that you can install the needed Bio++ libraries:
```
sudo zypper install libbpp-core-devel libbpp-phyl-devel 
libbpp-seq-devel libbpp-seq-omics-devel 
libbpp-phyl-omics-devel
```

Later on, at the build step you will most likely get some `undefined reference to` errors. If so, try to build the Bio++ libraries from source.

# 2 Clone ALE repository

```
git clone https://github.com/ssolo/ALE.git
```

# 3 Build ALE

Create a build directory and enter it:
```
mkdir ALE-build
cd ALE-build
```

For some reason ( [maybe because packaging guidelines](https://en.opensuse.org/openSUSE:Packaging_guidelines#Static_Libraries) ) the openSUSE Boost library can not be statically linked, so we have to use the `-DBUILD_STATIC=OFF` switch.

```
cmake ../ALE -DBUILD_STATIC=OFF
make
```

If everything went well, you will find the binaries in the `bin` directory.
