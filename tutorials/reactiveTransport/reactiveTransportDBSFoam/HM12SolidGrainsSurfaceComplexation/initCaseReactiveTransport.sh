#!/bin/bash

###### USERS INPUT ############################################################

#Species in the domain
species=("Ca+2","Cl-","H+","OH-")

#Diffusion coefficients (m^2/s), one per species
Diff=(1e-9,1e-9,1e-9,1e-9)

#Initial concentration
C0=(9.9705234e-05,0.00019941047,1.0166822e-09,1.0296592e-05)

#inlet concntration
C_in=(0,0,9.974142e-08,1.0095413e-07,0)
#C_in=(0,1.9940888e-05,1.0034273e-08,1.0156911e-06,9.970394e-06)

#Surface master species and site density (kmol/m2)
Surface_masters=("Surf_a")
master_density=(2.4e-9)

#Surface species and initial concentrations (kmol/m2)
Surf_species=("Surf_a-","Surf_aCa+","Surf_aH")
Surf_C0=(8.7802555e-10,8.4844109e-10,6.7353336e-10)

model='molDiff' ##Deff=Diff

## model='tortuosity'
## D=Diff/tau

#model='PecletDependent' 
## Dipsersivity constant
## Pe<1
## Dx=Diff/tau*(1+betax*Pe^alpha1x)
## Dy=Diff/tau*(1+betay*Pe^alpha1y)
## Dz=Diff/tau*(1+betaz*Pe^alpha1z)
## Pe>1
## Dx=Diff/tau*(1+betax**Pe^alpha2x)
## Dy=Diff/tau*(1+betay*Pe^alpha2y)
## Dz=Diff/tau*(1+betaz*Pe^alpha2z)


#define the labels of the phases
phases=(0,1,2)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=(0.0001,0.4,1)

##tortuosity factor of micropores
tau=(2.5,2.5,1)

##micro pore size (m) - Pe = UL/Diff/eps
#Lpore=(5e-6,5e-6,5e-6)

#betax=(0.5,0.5,0)
#alpha1x=(1,1,1)
#alpha2x=(1,1,1)
#betay=(0.5,0.5,0)
#alpha1y=(1,1,1)
#alpha2y=(1,1,1)
#betaz=(0.5,0.5,0)
#alpha1z=(1,1,1)
#alpha2z=(1,1,1)

#### END OF USER INPUT


IFS=',' read -ra species_array <<< "${species[0]}"
IFS=',' read -ra Diff_array <<< "${Diff[0]}"
IFS=',' read -ra C0_array <<< "${C0[0]}"
IFS=',' read -ra C_in_array <<< "${C_in[0]}"
IFS=',' read -ra Surface_masters_array <<< "${Surface_masters[0]}"
IFS=',' read -ra master_density_array <<< "${master_density[0]}"
IFS=',' read -ra Surf_species_array <<< "${Surf_species[0]}"
IFS=',' read -ra Surf_C0_array <<< "${Surf_C0[0]}"

python $GCFOAM_DIR/applications/utilities/pyTools/createThermoPhysicalProperties.py $species $Diff $Surf_species $Surface_masters

NPX="$(tail -n 1 system/NPX)"
NPY="$(tail -n 1 system/NPY)"
NPZ="$(tail -n 1 system/NPZ)"

if { [ -f constant/polyMesh/boundary ] && grep -q "frontAndBack" constant/polyMesh/boundary; } || \
   { [ -f processor0/constant/polyMesh/boundary ] && grep -q "frontAndBack" processor0/constant/polyMesh/boundary; }; then
   dimension="2D"
else dimension="3D"

fi

echo -e "create concentration"
if [ -d "processor0" ]
then
    export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"

    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then
        for i in "${!species_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done
        for i in "${!Surface_masters_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            srun python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done

	if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
          for i in "${!species_array[@]}"; do
              srun python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
          done
	elif [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
          for i in "${!species_array[@]}"; do
              srun python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
          done
	fi
    else
        for i in "${!species_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done

        for i in "${!Surface_masters_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done
        if [[ "$model" == "PecletDependent" ]];then
          echo "Calculate micro dispersivity"
            for i in "${!species_array[@]}"; do
                mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
            done
        elif  [[ "$model" == "tortuosity" ]];then
          echo "Calculate micro dispersivity"
            for i in "${!species_array[@]}"; do
                mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
            done
	fi
    fi
else
        for i in "${!species_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createC.py $dimension $NPX $NPY $NPZ "${species_array[$i]}" "${C0_array[$i]}" "${C_in_array[$i]}"
        done
        for i in "${!Surface_masters_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surface_masters_array[$i]}" "${master_density_array[$i]}"
        done
        for i in "${!Surf_species_array[@]}"; do
            python $GCFOAM_DIR/applications/utilities/pyTools/createSurf.py $dimension $NPX $NPY $NPZ "${Surf_species_array[$i]}" "${Surf_C0_array[$i]}"
        done


    if [[ "$model" == "PecletDependent" ]];then
      echo "Calculate micro dispersivity"
      for i in "${!species_array[@]}"; do
          python $GCFOAM_DIR/applications/utilities/pyTools/createDeffPeDependent.py $dimension "${Diff_array[$i]}" $micro_por $tau $Lpore $betax $alpha1x $alpha2x $betay $alpha1y $alpha2y $betaz $alpha1z $alpha2z $NPX $NPY $NPZ "DT_${species_array[$i]}"
      done
    elif [[ "$model" == "tortuosity" ]];then
      echo "Calculate micro dispersivity"
      for i in "${!species_array[@]}"; do
          python $GCFOAM_DIR/applications/utilities/pyTools/createDefftort.py $dimension "${Diff_array[$i]}" $micro_por $tau $NPX $NPY $NPZ "D_${species_array[$i]}"
      done
    fi
fi



echo -e "Case initialised"
