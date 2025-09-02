#!/bin/bash

###### USERS INPUT ############################################################

#Define boudnary tupe
#boundary_type="flow_rate"
boundary_type="pressure_drop"
#Define flow rate
flowRate=0 #4.2e-9
pressureDrop=0.00005 #0.1
#fluid properties
Visc=1e-6 #8.2e-8

model='phases'

#define the labels of the phases
phases=(0,1,2)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=(0.0001,0.4,1)

#define the permeability of each label (note: solid phase should be < 1e-20, pore should be > 1e6)
micro_k=(1e-16,9.876543209876547e-10,1e30)

#Kozeny-Carman constant
#kf=0

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
