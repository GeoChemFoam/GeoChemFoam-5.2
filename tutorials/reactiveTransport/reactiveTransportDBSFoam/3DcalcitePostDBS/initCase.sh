#!/bin/bash

###### USERS INPUT ############################################################

#Define boudnary tupe
boundary_type="flow_rate"

#Define flow rate
flowRate=3.5e-10

#fluid properties
Visc=2.61e-6
Diff=5e-9

model='grayscale'

#Reaction constants
kreac=8.9125e-4 
scoeff=2
rhos=2710
Mws=100
cinlet=0.0126


#Kozeny-Carman constant
kf=1.8e12

#### END OF USER INPUT #######################################################

if [[ $boundary_type != 'flow_rate' ]];then
        echo "ERROR: only flow_rate boundary type is implemented for this solver"
        exit
fi

pressureDrop=0

if [[ $model == 'phases' ]];then
        echo "ERROR: only grayscale segmentation is implemented for this solver"
        exit
fi

echo -e "set flow and transport properties"
cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/rho_s/$rhos/g" constant/transportProperties
sed -i "s/Mw_s/$Mws/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

cp constant/thermoPhysicalProperties1 constant/thermoPhysicalProperties
sed -i "s/Diff/$Diff/g" constant/thermoPhysicalProperties
sed -i "s/s_coeff/$scoeff/g" constant/thermoPhysicalProperties
sed -i "s/k_reac/$kreac/g" constant/thermoPhysicalProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f 0/eps ] && grep -q "frontAndBack" 0/eps; } || \
   { [ -f processor0/0/eps ] && grep -q "frontAndBack" processor0/0/eps; }; then
   dimension="2D"
else dimension="3D"

fi

echo -e "create U and p"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then

        srun python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ $cinlet
    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ $cinlet
    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
    python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop
    python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ $cinlet
fi

echo "Case initialised"
