#!/bin/bash

###### USERS INPUT ############################################################

#direction of strain
direction=1
#confining pressure (Pa)
pconf=2e7
#strain displacement in axis (m)
epsilon1=1e-6
#Solid density
rhos=2650
#Young modulus and poisson ratio
E=8e10
nu=0.15
#fluid kinematic pressurei (p/rho) (m2/s2) 
fluidPressure=2e6

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/mechanicalProperties1 constant/mechanicalProperties
sed -i "s/rho_s/$rhos/g" constant/mechanicalProperties
sed -i "s/Young_Modulus/$E/g" constant/mechanicalProperties
sed -i "s/Poisson/$nu/g" constant/mechanicalProperties
sed -i "s/fluid_pressure/$fluidPressure/g" constant/mechanicalProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f 0/eps ] && grep -q "frontAndBack" 0/eps; } || \
   { [ -f processor0/0/eps ] && grep -q "frontAndBack" processor0/0/eps; }; then
   dimension="2D"
else dimension="3D"

fi

isSolidWalls=0

echo -e "create displacement D"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createDisplacement.py $dimension $NPX $NPY $NPZ $direction $pconf $epsilon1 $isSolidWalls $fluidPressure

    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDisplacement.py $dimension $NPX $NPY $NPZ $direction $pconf $epsilon1 $isSolidWalls $fluidPressure
    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createDisplacement.py $dimension $NPX $NPY $NPZ $direction $pconf $epsilon1 $isSolidWalls $fluidPressure
fi
             
echo -e "Case initialised"

