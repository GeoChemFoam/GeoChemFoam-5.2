#!/bin/bash

###### USERS INPUT ############################################################

## Define the total number of iterations of the simulation
TotalTime=1000

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
sed -i "s/TotalTime/$TotalTime/g" system/controlDict
sed -i "s/WriteTimestep/$TotalTime/g" system/controlDict
sed -i "s/runTimestep/1/g" system/controlDict

cp system/fvSolution1 system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    # Run solidDisplacementFoam in parallel
    echo -e "Run solidDisplacementFoam in parallel on $NP processors"
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun solidDisplacementElasticVOSFoam -parallel  > solidDisplacementFoamStress.out
    else
        mpirun -np $NP solidDisplacementElasticVOSFoam -parallel  > solidDisplacementFoamStress.out
    fi
else
    echo -e "Run solidDisplacementFoam"
    solidDisplacementElasticVOSFoam > solidDisplacementFoamStress.out
fi

echo -e "Note: Please check the last line of solidDisplacementFoamStress.out to confirm the displacement field has converged. If it has not, change TotalTime and/or the D tolerance  in system/fvSolution and re-run the script" 
