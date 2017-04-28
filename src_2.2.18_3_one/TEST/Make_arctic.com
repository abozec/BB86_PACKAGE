#
set echo
cd $cwd
#
# --- Usage:  ./Make_arctic.com >& Make_arcticlog
#
# --- make arctic with TYPE from this directory's name (src_*_$TYPE/TEST).
# --- set ARCH to the correct value for this machine.
# --- assumes dimensions.h is correct for $TYPE.
#
setenv ARCH Asp4
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
  foreach t ( xct_arctic )
    make $t ARCH=$ARCH TYPE=$TYPE
  end
endif
