#
set echo
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
# --- make all with TYPE from this directory's name (src_*_$TYPE/TEST).
# --- set ARCH to the correct value for this machine.
# --- assumes dimensions.h is correct for $TYPE.
#
setenv ARCH sp3GPFS
#
cd ..
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
cd TEST
#
if (! -e ../../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
else
  cp ../mod*.[om]* .
  foreach t ( xca xcl xcs xct zaio )
    make $t ARCH=$ARCH TYPE=$TYPE
  end
endif
