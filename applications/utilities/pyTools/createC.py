import sys
from mpi4py import MPI
import numpy as np
import array
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dimension=sys.argv[1]
NPX=int(sys.argv[2])
NPY=int(sys.argv[3])
NPZ=int(sys.argv[4])

speciesName=sys.argv[5]
cinlet=float(sys.argv[6])



NP=NPX*NPY*NPZ

ipz = rank // (NPX * NPY)  # Get ipz by integer division
remainder = rank % (NPX * NPY)
ipy = remainder // NPX     # Get ipy
ipx = remainder % NPX      # Get ipx

                
###################################################################
###### p ##########################################################
###################################################################
 
if NP>1:
    output_path = "processor"+str(rank)+"/"
else:
    output_path = ""

data = [
'/*--------------------------------*- C++ -*----------------------------------*\\\\n',
'| =========                 |                                                 |\n',
'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
'|  \\\\    /   O peration     | Version:  2212                                  |\n',
'|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n',
'|    \\\\/     M anipulation  |                                                 |\n',
'\\*---------------------------------------------------------------------------*/\n',
'FoamFile\n',
"{\n",
"    version     2.0;\n",
"    format      ascii;\n",
"    arch        \"LSB;label=32;scalar=64\";\n",
"    class       volScalarField;\n",
"    location    \"0\";\n",
"    object      "+speciesName+";\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"dimensions      [0 -3 0 0 1 0 0];\n",
"\n",
f"internalField   uniform 0;\n",
'\n'
"boundaryField"'\n',
"{"'\n',
"    inlet"'\n',
"    {"'\n'
"        type            fixedValue;\n",
f"        value           uniform {str(cinlet)};\n",
"    }\n",
"    outlet\n"
"    {\n"
"        type           zeroGradient;\n",
"    }\n"
"    \"wall_.*\"""\n"
"    {\n",
"        type            zeroGradient;\n",
"    }\n"
]

if (dimension=="2D"):
   data.extend([
   "    frontAndBack\n"
   "    {\n",
   "        type            empty;\n",
   "    }\n"
   ])

##proc k to k-1
if (ipz>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])
            
##proc j to j-1                
if (ipy>0):
   data.extend([
   f"    procBoundary{rank}to{rank - NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc i to i-1                
if (ipx>0):
   data.extend([
   f"    procBoundary{rank}to{rank - 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc i to i+1                
if (ipx<NPX-1):
   data.extend([
   f"    procBoundary{rank}to{rank + 1}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc j to j+1                
if (ipy<NPY-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc k to k+1                
if (ipz<NPZ-1):
   data.extend([
   f"    procBoundary{rank}to{rank + NPX * NPY}\n",
   "    {\n",
   "        type            processor;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])
data.extend([
"}\n",
"\n",
"\n",
"// ************************************************************************* //"
])


#io_start_time = time.time()

# Write data to file
with open(output_path+'0/'+speciesName, 'w') as f:
  f.writelines(data)
  f.close()
