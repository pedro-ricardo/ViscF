
## AMReX Dependency

To install AMReX software framework, clone it with git
```bash
git clone https://github.com/AMReX-Codes/amrex.git
```

Than enter it and change the **branch** to the release you want
```bash
cd amrex
git checkout 20.10
```

There is a bug in the CMake configuration in the latest version (`20.10`). To fix it go to the file **Tools/CMake/FindPETSc.cmake** in line 40 and comment out:
```cmake
# set(CMAKE_Fortran_FLAGS ${PETSC_INCLUDE_DIRS})
```

### Cmake install

Create a build dir in the AMReX root and go into it
```bash
mkdir build
cd build
```

Get the libraries install path for `MPI`, `HYPRE`, `HDF5` and `PETSC`.
- The **Hypre** path is usualy the same as **PETSc**
- The **HDF5** lib has to be build with the flags `--enable-cxx`, `--enable-parallel` and `--enable-unsupported`

Configure CMake with flags, adding your instalation folders for the libraries:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/amrex \
-DDIM=3 -DENABLE_EB=yes -DENABLE_PARTICLES=yes \
-DENABLE_FORTRAN=yes -DENABLE_FORTRAN_INTERFACES=yes \
-DENABLE_MPI=yes -DMPI_ROOT=/path/to/mpich \
-DENABLE_HYPRE=yes -DHYPRE_ROOT=/path/to/petsc \
-DENABLE_HDF5=yes -DHDF5_ROOT=/path/to/hdf5 \
-DENABLE_PETSC=yes -DPETSC_DIR=/path/to/petsc \
-DBUILD_SHARED_LIBS=true -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=true ..
```

Now do `make -j` and `make install`
