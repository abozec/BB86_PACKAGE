#!/bin/csh
#
set echo
set time = 1
set timestamp
#
# --- test xctilr, single cpu
#
touch   regional.depth.a
/bin/rm regional.depth.a
ln -s ../../topo/depth_GLBa10.0_06.a regional.depth.a
./test_xct_arctic
