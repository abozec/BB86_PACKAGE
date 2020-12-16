#!/bin/csh

set echo
## CONFIGURATION BOX 20kmx20km
setenv R BB86   ## config name
setenv T 01     ## topography version (see name of your file..)
setenv E 010    ## experiment number (ex: 010)
setenv E1 `echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`   ## (ex: 01.0)

## time to run
setenv day1  1 
setenv day2  360
## restart run?
setenv restart -   ## - means no restart (+ means restart)

## get all the files in the same directory
# P : path where you are going to launch the run
# D : path where you are actually running 
setenv P /Net/yucatan/abozec/BB86_PACKAGE/expt_${E1}/
setenv D ${P}/data

mkdir -p ${D}

## get the grid and depth files
touch   ${D}/regional.grid.[ab] ${D}/regional.depth.[ab]
/bin/rm ${D}/regional.grid.[ab] ${D}/regional.depth.[ab]

/bin/cp ${P}/../topo/regional.grid.${R}.a  ${D}/regional.grid.a
/bin/cp ${P}/../topo/regional.grid.${R}.b  ${D}/regional.grid.b
/bin/cp ${P}/../topo/depth_${R}_${T}.a     ${D}/regional.depth.a
/bin/cp ${P}/../topo/depth_${R}_${T}.b     ${D}/regional.depth.b

## get the initial conditions
touch   ${D}/relax.temp.[ab] ${D}/relax.saln.[ab] ${D}/relax.intf.[ab]
/bin/rm ${D}/relax.temp.[ab] ${D}/relax.saln.[ab] ${D}/relax.intf.[ab]

#potential temperature
/bin/cp ${P}/../relax/${E}/relax_tem_${R}.a ${D}/relax.temp.a
/bin/cp ${P}/../relax/${E}/relax_tem_${R}.b ${D}/relax.temp.b

#salinity
/bin/cp ${P}/../relax/${E}/relax_sal_${R}.a ${D}/relax.saln.a
/bin/cp ${P}/../relax/${E}/relax_sal_${R}.b ${D}/relax.saln.b

#interface depth
/bin/cp ${P}/../relax/${E}/relax_int_${R}.a ${D}/relax.intf.a
/bin/cp ${P}/../relax/${E}/relax_int_${R}.b ${D}/relax.intf.b

# relaxation rmu
/bin/cp ${P}/../relax/${E}/relax_rmu.a ${D}/relax.rmu.a
/bin/cp ${P}/../relax/${E}/relax_rmu.b ${D}/relax.rmu.b


## to not get into trouble with some flags ..
touch ${D}/relax.weird 

## get the forcing
touch   ${D}/forcing.tauewd.[ab] ${D}/forcing.taunwd.[ab]
/bin/rm ${D}/forcing.tauewd.[ab] ${D}/forcing.taunwd.[ab]

#eastward windstress
/bin/cp ${P}/../force/forcing.tauewd.${R}.a ${D}/forcing.tauewd.a
/bin/cp ${P}/../force/forcing.tauewd.${R}.b ${D}/forcing.tauewd.b

#northward windstress
/bin/cp ${P}/../force/forcing.taunwd.${R}.a ${D}/forcing.taunwd.a
/bin/cp ${P}/../force/forcing.taunwd.${R}.b ${D}/forcing.taunwd.b

## get the executable
/bin/cp ${P}/../src_2.3.01_one/hycom ${D}/hycom

## get rid of the old files 
/bin/rm ${D}/restart_out* ${D}/ovrtn_out 


## Parameter file
/bin/cp ${P}/blkdat.input ${D}/blkdat.input

## length of the run
cat > limits <<EOF
${restart}${day1}  ${day2}
EOF
/bin/cp ${P}/limits ${D}/limits


## let''s run
cd ${D}
./hycom

## Put archives in OUTPUT directory
mkdir  -p ${D}/output
/bin/mv ${D}/archv* ${D}/output/


echo 'You are Done !'


