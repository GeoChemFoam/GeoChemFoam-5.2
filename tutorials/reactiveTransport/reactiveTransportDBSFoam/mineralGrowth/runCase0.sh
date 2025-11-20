#!/bin/bash

###### USERS INPUT ############################################################

#### END OF USER INPUT #######################################################

cp system/controlDict1 system/controlDict
cp system/fvSolution1 system/fvSolution

sed -i "s/TotalTime/1e-06/g" system/controlDict
sed -i "s/WriteTimestep/1e-06/g" system/controlDict
sed -i "s/initTimestep/1e-06/g" system/controlDict
sed -i "s/maxTimestep/1e-06/g" system/controlDict
sed -i "s/NCORR/5000/g" system/fvSolution

if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    echo -e "run reactiveTransportDBSFoam for 1e-06 sec in parallel on $NP processors"
    mpirun -np $NP reactiveTransportDBSFoam -parallel > reactiveTransportDBSFoam0.out

    echo -e "move 1e-06 to 0"
    mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/moveResultTo0.py "U" "p" "phi" "A" "B" "R"

else
    echo -e "run reactiveTransportDBSFoam for 1e-06"
    reactiveTransportDBSFoam > reactiveTransportDBSFoam0.out

    cp -f 1e-06/A 0/.
    cp -f 1e-06/B 0/.
    cp -f 1e-06/U 0/.
    cp -f 1e-06/p 0/.
    cp -f 1e-06/phi 0/.
    cp -f 1e-06/R 0/.
    rm -rf 1e-06
fi

echo -e "Note: Please check the last line of reactiveTransportDBSFoam0.out to confirm the equation have converged. If it has not, re-run script and/or change the tolerance and residual controls in system/fvSolution" 


