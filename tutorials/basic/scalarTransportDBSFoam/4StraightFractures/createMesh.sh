#!/bin/bash

###### USERS INPUT ############################################################

#Define image name
Image_name="4StraightFractures"
#Image_name="Bentheimer400-5mum"
#Image_name="Est_3phase500cubed4micron"

#Image directory location ($PWD if current directory)
dir="$GCFOAM_IMG/raw"

#Choose image format
format='raw'

#Choose if the image is compressed or not
compressed='yes'

# Define image dimensions
x_dim=500
y_dim=500
z_dim=1

#segmentation='grayscale'
#Values of solid, pore, and minimum porosity value for the solid phase (note: if the image contains solid voxels, this CANNOT be 0)
#pores_value=255
#solid_value=0
#eps_min=0.0001

segmentation='phases'
#Values of solid and pore
pores_value=2
solid_value=0

#define the labels of the phases
phases=(0,1,2)

#define the porosity of each phase, note that the porosity of the solid phase CANNOT be 0, default to 0.0001
micro_por=(0.0001,0.4,1)

# Define cropping parameters
x_min=0
x_max=500
y_min=0
y_max=500
z_min=0
z_max=1

#padding for inlet/outlet
#This adds voxels in the image directly
padWidth=0

# define resolution (m)
res=0.0001

# number of cells of initial mesh
## They need to be a factor of the number of voxels in the padded image
n_x=250
n_y=250
n_z=1

#Mesh refinement level
localRefinement='false'
#number of layers around interface
layer=0
nlevel=0
refineStokes=0

#dimension="3D"
dimension="2D"
# flow direction 0, 1 or 2 
direction=0

# Number of processors in each direction
NPX=4
NPY=2
NPZ=1

#### END OF USER INPUT #######################################################

if [ $format != 'raw' ]
then
        echo "ERROR: only raw format is implemented for this solver"
        exit
fi

if [[ "${dimension}" == "2D" ]]
then
  if [ $direction == 2 ]
  then
     echo "ERROR: direction cannot be 2 in 2D simulations"
     exit
  elif [ $n_z -gt 1 ]
  then
     echo "ERROR: n_z= $n_z must be 1 in 2D simulations"
     exit
  elif [ $NPZ -gt 1 ]
  then
     echo "EROR: NPZ must be 1 in 2D simulations"
     exit
  elif [ $localRefinement == 'true' ]
  then
     echo "ERROR: no dynamicmeshRef for 2D simulations"
     exit
  fi
fi

#Insert dimensions in postProcessDict
x_1=0
y_1=0
z_1=0

xSize=$(echo "$x_max - $x_min" | bc)
ySize=$(echo "$y_max - $y_min" | bc)
zSize=$(echo "$z_max - $z_min" | bc)

x_2=$(expr $xSize*$res | bc)
y_2=$(expr $ySize*$res | bc)
z_2=$(expr $zSize*$res | bc)

cp system/postProcessDict1 system/postProcessDict
sed -i "s/x_1/$x_1/g" system/postProcessDict
sed -i "s/y_1/$y_1/g" system/postProcessDict
sed -i "s/z_1/$z_1/g" system/postProcessDict

sed -i "s/x_2/$x_2/g" system/postProcessDict
sed -i "s/y_2/$y_2/g" system/postProcessDict
sed -i "s/z_2/$z_2/g" system/postProcessDict

sed -i "s/flowdir/$direction/g" system/postProcessDict

mkdir -p constant/polyMesh

mkdir -p 0

filename=$Image_name\.$format

if [ $compressed == 'yes' ]
then
        filename=$Image_name\.$format\.tar.gz
fi

mkdir -p constant/triSurface
cp $dir\/$filename constant/triSurface/.
cd constant/triSurface
if [ $compressed == 'yes' ]
then
        tar -xf $filename
fi
cd ../..

NP=$(echo "$NPX*$NPY*$NPZ" | bc)

rm -f system/NP*

echo $NPX >> system/NPX
echo $NPY >> system/NPY
echo $NPZ >> system/NPZ

if [ $segmentation == 'grayscale' ]
then
	micro_por=0
	phases=0
elif [ $segmentation == 'phases' ]
then
    eps_min=0
fi

cyclic='no'

if [[ $NP > 1 ]]; then
    # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
    if [[ "${PLATFORM}" == "ARCHER2" ]]; then

        # error check that number of process files matches number of MPI tasks in parallel job
        if [[ ${NP} != ${SLURM_NTASKS} ]]
        then
           echo "ERROR: Number of MPI tasks does not equal number of processor files"
           exit
        fi

        # prepend spindle to srun command to pre-load python modules, otherwise comment-out sprindle line to simply use srun
        srun python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

    else
       mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

    fi
else
    python $GCFOAM_DIR/applications/utilities/pyTools/createMesh.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $NPX $NPY $NPZ

fi

if [ $localRefinement == 'true' ]
then
    echo -e "create topoSetDict"
    python $GCFOAM_DIR/applications/utilities/pyTools/createTopoSet.py $segmentation $eps_min $micro_por $refineStokes
    echo -e "refine mesh"

    cp system/fvSolution1 system/fvSolution
    sed -i "s/nSmooth/1/g" system/fvSolution
    sed -i "s/cSmooth/0.5/g" system/fvSolution

    for ((i=1; i < nlevel+1; i++))
    do
        if [[ $NP > 1 ]]; then
            # if PLATFORM is ARCHER2 then use srun, otherwise use serial version
            if [[ "${PLATFORM}" == "ARCHER2" ]]; then

                # error check that number of process files matches number of MPI tasks in parallel job
                if [[ ${NP} != ${SLURM_NTASKS} ]]
                then
                   echo "ERROR: Number of MPI tasks does not equal number of processor files"
                   exit
                fi

                for ((j=0; j < layer-1; j++))
                do
                  echo -e "smooth solid surface"
                  srun smoothSolidSurface -parallel > smoothSolidSurface.out
                done

                echo -e "refine fluid/solid interface"
                srun topoSet -parallel > topoSet.out
                srun refineHexMesh refinementRegion -overwrite -parallel > refineHexMesh.out

                echo -e "processMeshCellCenters"
                srun processMeshCellCenters -parallel > processMeshCellCenters.out

                echo "create refined eps"
                srun python $GCFOAM_DIR/applications/utilities/pyTools/createEpsRefinement.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $i $refineStokes $NPX $NPY $NPZ

            else
                for ((j=0; j < layer-1; j++))
                do
                    echo -e "smooth solid surface"
                    mpirun -np $NP smoothSolidSurface -parallel > smoothSolidSurface.out
                done

                echo -e "refine fluid/solid interface"
                mpirun -np $NP topoSet -parallel > topoSet.out
                mpirun -np $NP refineHexMesh refinementRegion -overwrite -parallel > refineHexMesh.out

                echo -e "processMeshCellCenters"
                mpirun -np $NP processMeshCellCenters -parallel > processMeshCellCenters.out

                echo "create refined eps"
                mpirun -np $NP python $GCFOAM_DIR/applications/utilities/pyTools/createEpsRefinement.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $i $refineStokes $NPX $NPY $NPZ

            fi

            rm -rf processor*/0/cellCenters
        else
            for ((j=0; j < layer-1; j++))
            do
                echo -e "smooth solid surface"
                smoothSolidSurface > smoothSolidSurface.out
            done

            echo -e "refine fluid/solid interface"
            topoSet > topoSet.out
            refineHexMesh refinementRegion -overwrite > refineHexMesh.out

            echo -e "processMeshCellCenters"
            processMeshCellCenters > processMeshCellCenters.out

            echo "create refined eps"
            python $GCFOAM_DIR/applications/utilities/pyTools/createEpsRefinement.py $x_dim $y_dim $z_dim $x_min $x_max $y_min $y_max $z_min $z_max $n_x $n_y $n_z $res $Image_name $padWidth $pores_value $solid_value $eps_min $dimension $direction $cyclic $segmentation $micro_por $phases $i $refineStokes $NPX $NPY $NPZ

            rm -rf 0/cellCenters
        fi
    done
fi

rm -rf constant/triSurface

echo -e "Mesh created. It is advised to check in paraview to confirm mesh and 0/eps are reasonable before running flow"
