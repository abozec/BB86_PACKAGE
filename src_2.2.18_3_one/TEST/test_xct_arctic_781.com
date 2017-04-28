#!/bin/csh
#
#@ job_name         = test_xct_arctic_781
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = parallel
#@ network.MPI      = csss,not_shared,US
#@ environment      = MP_EUILIB=us     
#@ node             = 98
#@ total_tasks      = 781
#@ node_usage       = not_shared
#@ resources        = ConsumableCpus(1) ConsumableMemory(500mb)
#@ wall_clock_limit = 0:30:00
#@ account_no       = NRLSSC1J
#@ class            = block2
#@ queue
#
set echo
set time = 1
set timestamp
#
# --- test xctilr, single cpu
#
touch   regional.depth.a
/bin/rm regional.depth.a
ln -s ../../topo/depth_GLBa0.08_05.a           regional.depth.a
cp    ../../topo/partit/depth_GLBa0.08_05.781 patch.input
#
#   AIX
#
    setenv MP_SHARED_MEMORY     yes
    setenv MP_SINGLE_THREAD     yes
    setenv MP_EAGER_LIMIT       65536
#   list where the MPI job will run
    env MP_LABELIO=YES poe hostname
    poe ./test_xct_arctic
