#!/bin/bash

###### USERS INPUT ############################################################

#Initial Temperature and boundary condition (dimensionless)
T0=1.0
Tin=0.0
wall_boundary_type="fixedValue" # This value will be T0
#wall_boundary_type="zeroGradient"

#Diffusion model
diffModel="harmonic"
#diffModel="arithmetic"
#fluid properties
#Heat conductivity (kW/m/K)
kappa_f=0.04 
kappa_s=3.3

#Density (kg/m3)
rho_f=300
rho_s=2600

#Heat capacity (kJ/kg/K)
gamma_f=3000
gamma_s=700

#### END OF USER INPUT #######################################################

echo -e "set flow and transport properties"
cp constant/transportProperties2 constant/transportProperties
sed -i "s/diffModel/$diffModel/g" constant/transportProperties
sed -i "s/kappa_s/$kappa_s/g" constant/transportProperties
sed -i "s/kappa_f/$kappa_f/g" constant/transportProperties
sed -i "s/rho_f/$rho_f/g" constant/transportProperties
sed -i "s/rho_s/$rho_s/g" constant/transportProperties
sed -i "s/gamma_f/$gamma_f/g" constant/transportProperties
sed -i "s/gamma_s/$gamma_s/g" constant/transportProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f 0/eps ] && grep -q "frontAndBack" 0/eps; } || \
   { [ -f processor0/0/eps ] && grep -q "frontAndBack" processor0/0/eps; }; then
   dimension="2D"
else dimension="3D"
fi

echo -e "create T"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"
    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"
    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ $T0 $Tin $wall_boundary_type "zeroGradient"
fi

echo -e "Case initialised"


