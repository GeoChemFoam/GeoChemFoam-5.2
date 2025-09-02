#!/bin/bash

###### USERS INPUT ############################################################

#Diffusion coefficient (m^2/s)
Diff=1e-8

#model='molDiff' ##Deff=Diff

#model='tortuosity'

model='PecletDependent' 
## Dipsersivity constant
## Pe<1
## Dx=Diff/tau*(1+betax**Pe^alpha1x)
## Dy=Diff/tau*(1+betay*Pe^alpha1y)
## Dz=Diff/tau*(1+betaz*Pe^alpha1z)
## Pe>1
## Dx=Diff/tau*(1+betax**Pe^alpha2x)
## Dy=Diff/tau*(1+betay*Pe^alpha2y)
## Dz=Diff/tau*(1+betaz*Pe^alpha2z)


#define the labels of the phases
phases=(1,2,3)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=(1.0,0.445167,0.0001)

##tortuosity iof micropores
tau=(1,2.636,0)

##micro pore size - Pe = UL/Diff/eps
Lpore=(3.69828e-7,3.69828e-07,0)

betax=(0,70,70)
alpha1x=(1.35,1.35,1.35)
alpha2x=(1.35,1.35,1.35)
betay=(0,4.6,4.6)
alpha1y=(1.0,1.0,1.0)
alpha2y=(0.68,0.68,0.68)
betaz=(0,4.6,4.6)
alpha1z=(1.0,1.0,1.0)
alpha2z=(0.68,0.68,0.68)

#### END OF USER INPUT #######################################################

rm -f constant/fvOptions

cp constant/transportProperties2 constant/transportProperties
sed -i "s/Diff/$Diff/g" constant/transportProperties

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f 0/eps ] && grep -q "frontAndBack" 0/eps; } || \
   { [ -f processor0/0/eps ] && grep -q "frontAndBack" processor0/0/eps; }; then
   dimension="2D"
else dimension="3D"

fi

echo "Create vector B"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    if [[ "${PLATFORM}" == "ARCHER2" ]]; then    
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createBCyclic.py $dimension $NPX $NPY $NPZ

	if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependentCyclic.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ
	elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro effective diffusion coefficient"
          srun python $GCFOAM_DIR/applications/utilities/pyTools/createDefftortCyclic.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
        fi

   else
        mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createBCyclic.py $dimension $NPX $NPY $NPZ

        if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependentCyclic.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ 
	elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro effective diffusion coefficient"
          mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDefftortCyclic.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
        fi
   fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createBCyclic.py $dimension $NPX $NPY $NPZ

    if [[ "$model" == "PecletDependent" ]];then
      echo "Calculate micro dispersivity"
      python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependentCyclic.py $dimension $Diff $micro_por $tau $Lpore $betax $alpha1x $alpha2x  $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ
    elif [[ "$model" == "tortuosity" ]];then
      echo "Calculate micro effective diffusion coefficient"
      python $GCFOAM_DIR/applications/utilities/pyTools/createDefftortCyclic.py $dimension $Diff $micro_por $tau $NPX $NPY $NPZ
    fi
fi

echo -e "Case initialised"


