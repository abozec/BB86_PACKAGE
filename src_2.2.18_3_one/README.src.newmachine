src_2.0.01_16/README.src.newmachine:

The Makefile sources ../config/$(ARCH)_$(TYPE) where ARCH defines exactly 
what machine architecture to target and TYPE is the parallelization 
strategy and precision (one, omp, mpi, mpisr, shmem).  The make process 
is automated by the script Make.com, which should be used instead of 
directly invoking the make command.

The source code directory name should end with _${TYPE}, where ${TYPE}
is the parallelization type (one,omp,mpi,ompi,shmem).  The script
Make.com should be edited to define ${ARCH} appropriately for the 
machine.  The executable is then created by the command:

    ./Make.com >& Make.log

In order for this to work, the file config/${ARCH}_${TYPE} must exist
and must contain the machine-specific parts of Makefile (see README.config).
