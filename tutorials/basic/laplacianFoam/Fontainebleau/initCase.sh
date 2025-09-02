#!/bin/bash

###### USERS INPUT ############################################################

#Diffusion coefficient (m^2/s)
Diff=1e-6

model='molDiff' ##Deff=Diff

#model='tortuosity'
## D=Diff/tau

#define the labels of the phases
#phases=(1,2,3)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
#micro_por=(1.0,0.445167,0.0001)

##tortuosity iof micropores
#tau=(1,2.636,0)

#### END OF USER INPUT #######################################################

cp constant/transportProperties1 constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/transportProperties

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
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ 0.0 1.0 "zeroGradient" "fixedValue" 
	if [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro effective diffusion coefficient"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
        fi
    else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ 0.0 1.0 "zeroGradient" "fixedValue" 
	if [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro effective diffusion coefficient"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
        fi
    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createT.py $dimension $NPX $NPY $NPZ 0.0 1.0 "zeroGradient" "fixedValue" 
    if [[ "$model" == "tortuosity" ]];then
      echo "Calculate micro effective diffusion coefficient"
      python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
    fi
fi

echo -e "Case initialised"


