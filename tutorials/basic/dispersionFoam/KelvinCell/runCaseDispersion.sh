#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=10000
WriteTime=10000
RunTimeStep=1

#### END OF USER INPUT #######################################################

cp system/controlDict2 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTime/$WriteTime/g" system/controlDict
sed -i "s/RunTimeStep/$RunTimeStep/g" system/controlDict

cp system/fvSolution2 system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # Run dispersionFoam in parallel
    echo -e "Run dispersionFoam in parallel on $NP processors"
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun dispersionFoam -parallel  > dispersionFoamD.out
    else
        mpirun -np $NP dispersionFoam -parallel  > dispersionFoamD.out
    fi
else
    echo -e "Run dispersionFoam"
    dispersionFoam > dispersionFoamD.out
fi
