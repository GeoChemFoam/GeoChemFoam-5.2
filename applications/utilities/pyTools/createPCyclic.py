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
"    object      p;\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"dimensions      [0 2 -2 0 0 0 0];\n",
"\n",
"internalField   uniform 0;\n",
'\n'
"boundaryField"'\n',
"{"'\n',
"    left"'\n',
"    {"'\n'
"        type            cyclic;\n",
"    }\n"
"    right\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
"    bottom"'\n',
"    {"'\n'
"        type            cyclic;\n",
"    }\n"
"    top\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
]

if (dimension=="2D"):
   data.extend([
   "    frontAndBack\n"
   "    {\n",
   "        type            empty;\n",
   "    }\n"
   ])
else:
   data.extend([
   "    front"'\n',
   "    {"'\n'
   "        type            cyclic;\n",
   "    }\n"
   "    back\n"
   "    {\n",
   "        type            cyclic;\n",
   "    }\n"
   ])

##proc rank to rank-NPX*NPY*(NPZ-1) through back
if (ipz==NPZ-1) and (NPZ>2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           uniform 0;\n",
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
            

##proc rank to rank-NPX*NPY*(NPZ-1) through back
if (ipz==NPZ-1) and (NPZ==2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc rank to rank-NPX*(NPY-1) through top 
if (ipy==NPY-1) and (NPY>2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
   "    {\n",
   "        type            processorCyclic;\n",
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

##proc rank to rank-NPX*(NPY-1) through top 
if (ipy==NPY-1) and (NPY==2):
   data.extend([
   f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
   "    {\n",
   "        type            processorCyclic;\n",
   "        value           uniform 0;\n",
   "    }\n"
   ])

##proc rank to rank-(NPX-1) through right
if (ipx==NPX-1) and (NPX>2):
   data.extend([
   f"    procBoundary{rank}to{rank-(NPX-1)}throughright\n",
   "    {\n",
   "        type            processorCyclic;\n",
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

##proc rank to rank-(NPX-1) through right
if (ipx==NPX-1) and (NPX==2):
   data.extend([
   f"    procBoundary{rank}to{rank-(NPX-1)}throughright\n",
   "    {\n",
   "        type            processorCyclic;\n",
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

##proc rank to rank+(NPX-1) through left
if (ipx==0) and (NPX>1):
   data.extend([
   f"    procBoundary{rank}to{rank+(NPX-1)}throughleft\n",
   "    {\n",
   "        type            processorCyclic;\n",
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

##proc rank to rank+NPX*(NPY-1) through bottom
if (ipy==0) and (NPY>1):
   data.extend([
   f"    procBoundary{rank}to{rank+NPX*(NPY-1)}throughbottom\n",
   "    {\n",
   "        type            processorCyclic;\n",
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

##proc rank to rank+NPX*NPY*(NPZ-1) through front 
if (ipz==0) and (NPZ>1):
   data.extend([
   f"    procBoundary{rank}to{rank+NPX*NPY*(NPZ-1)}throughfront\n",
   "    {\n",
   "        type            processorCyclic;\n",
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
with open(output_path+'0/p', 'w') as f:
  f.writelines(data)
  f.close()
