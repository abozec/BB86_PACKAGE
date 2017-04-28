#!/bin/csh
#
#@ job_name         = test_all
#@ output           = $(job_name).log
#@ error            = $(job_name).log
#@ restart          = yes
#@ job_type         = parallel
#@ network.MPI      = css0,not_shared,US
#@ environment      = MP_EUILIB=us     
#@ node             = 1
#@ total_tasks      = 3
#@ node_usage       = not_shared
#@ wall_clock_limit = 0:15:00
#@ account_no       = NRLSS018
#@ class            = batch
#@ queue
#
set echo
set time = 1
set timestamp
#
# --- run all test cases.
#
cd ~/hycom/ATLa2.00/src_2.0.01_16_ompi/TEST
#
foreach f ( test_xca test_xcl test_xcs test_xct test_zaio )
  touch   ${f}.log
  /bin/rm ${f}.log
  csh ${f}.com >& ${f}.log
end
