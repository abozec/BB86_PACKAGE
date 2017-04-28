#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./Make_esmf.com >& Make_esmf.log
#
# --- make esmf (ESMF HYCOM component) with TYPE=esmf.
# --- this directory's name must be src_*_esmf.
# --- assumes dimensions.h is correct for TYPE=esmf (i.e. for mpi).
#
# --- set ARCH to the correct value for this machine.
# --- ARCH that start with A are for ARCTIC patch regions
#
#setenv ARCH alphaL
#setenv ARCH alpha
#setenv ARCH amd64
#setenv ARCH intel
#setenv ARCH o2k
#setenv ARCH sp3
#setenv ARCH sp4
#setenv ARCH sp5
#setenv ARCH sun64
#setenv ARCH sun
#setenv ARCH t3e
#setenv ARCH xt3
#
setenv ARCH sp4
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
#
if     ($TYPE != "esmf") then
  echo "TYPE must be esmf to invoke esmf make target"
  exit 1
endif
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
endif
#
# --- esmf needs additional environment variables.
#
if ($TYPE == "esmf") then
  switch ($ARCH)
  case 'sp4':
    setenv BEI_HOME /site/BEI
    setenv ESMF_DIR ${BEI_HOME}/esmf/2.2.2r
    breaksw
  case 'o2k':
    setenv BEI_HOME /usr/local/usp/BEI
    setenv ESMF_DIR ${BEI_HOME}/esmf/2.2.2r
    breaksw
  case 'xt3':
    setenv BEI_HOME /usr/local/usp/BEI
    setenv ESMF_DIR ${BEI_HOME}/esmf/2.2.2r
    breaksw
  default:
    echo "TYPE = esmf  needs BEI_HOME and ESMF_DIR"
    exit (1)
  endsw
endif
#
# --- some machines require gmake
#
#gmake ARCH=$ARCH TYPE=$TYPE esmf
make ARCH=$ARCH TYPE=$TYPE esmf
