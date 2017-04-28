#!/bin/csh
#
set echo
#
# --- run HYCOM through CPP.
#
#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./CPP.com >& CPP.log
#
# --- Build a source code set for this ARCH and TYPE only.
# --- assumes dimensions.h is correct for $TYPE.
#
#setenv ARCH alphaL
#setenv ARCH alpha
#setenv ARCH intel
#setenv ARCH o2k
#setenv ARCH sp3GPFS
#setenv ARCH sp3
#setenv ARCH sun64
#setenv ARCH sun
#setenv ARCH t3e
#
setenv ARCH sp3GPFS
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
else
  make hycom ARCH=$ARCH TYPE=$TYPE
endif
