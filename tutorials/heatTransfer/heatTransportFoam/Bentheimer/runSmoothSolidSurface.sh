#!/bin/bash

###### USERS INPUT ############################################################

#Smoothing parameters: smooth surface when image has artifical roughness created ny segmentation to avoid error when using adaptive mesh
nSmooth=0
cSmooth=0

#### END OF USER INPUT #######################################################

cp system/fvSolution1 system/fvSolution
sed -i "s/nSmooth/$nSmooth/g" system/fvSolution
sed -i "s/cSmooth/$cSmooth/g" system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
    echo -e "smoothSolidSurface in parallel on $NP processors"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        srun smoothSolidSurface -parallel  > smoothSolidSurface.out
    else
        mpiexec -np $NP smoothSolidSurface -parallel  > smoothSolidSurface.out
    fi
else
    echo -e "smoothSolidSurface"
    smoothSolidSurface > smoothSolidSurface.out
fi




