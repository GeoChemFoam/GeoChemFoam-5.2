#!/bin/bash

###### USERS INPUT ############################################################

#Define boundary type
#boundary_type="flow_rate"
boundary_type="pressure_drop"
#Enter value (m^3/s) if boundary_type is flow_rate, otherwise put 0
flowRate=0 #4.2e-9
#Enter value (Pa/rho) if boundary_type is pressure-drop, otherwise put 0 
pressureDrop=2.41e-1 
#fluid kinematic viscosity (m2/s)
Visc=1e-6

#model='phases'

#define the labels of the phases
#phases=(0,255)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0
#micro_por=(0.0001,1)

#define the permeability of each label (note: grain phase should be < 1e-20, pore should be > 1e6)
#micro_k=(1e-16,1e30)

#Kozeny-Carman constant
kf=1e12

#### END OF USER INPUT #######################################################

if [[ $model == 'phases' ]];then
	kf=0
fi

echo -e "set flow and transport properties"
cp constant/transportProperties1 constant/transportProperties
sed -i "s/Visc/$Visc/g" constant/transportProperties
sed -i "s/k_f/$kf/g" constant/transportProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f constant/polyMesh/boundary ] && grep -q "frontAndBack" constant/polyMesh/boundary; } || \
   { [ -f processor0/constant/polyMesh/boundary ] && grep -q "frontAndBack" processor0/constant/polyMesh/boundary; }; then
   dimension="2D"
else dimension="3D"

fi

echo -e "create U and p"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        if [[ "$model" == "phases" ]];then
          echo "Calculate micro permeability"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createKinvPhases.py $dimension $micro_por $micro_k $NPX $NPY $NPZ
        fi

        srun python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop

    else
        if [[ "$model" == "phases" ]];then
          echo "Calculate micro permeability"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createKinvPhases.py $dimension $micro_por $micro_k $NPX $NPY $NPZ
        fi

        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop
    fi
else
    if [[ "$model" == "phases" ]];then
          echo "Calculate micro permeability"
          python $GCFOAM_DIR/applications/utilities/pyTools/createKinvPhases.py $dimension $micro_por $micro_k $NPX $NPY $NPZ
    fi

    python $GCFOAM_DIR/applications/utilities/pyTools/createU.py $dimension $NPX $NPY $NPZ $boundary_type $flowRate
    python $GCFOAM_DIR/applications/utilities/pyTools/createP.py $dimension $NPX $NPY $NPZ $boundary_type $pressureDrop
fi
