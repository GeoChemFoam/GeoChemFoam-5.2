import numpy as np
import array
import os
import time
import sys
from mpi4py import MPI
import h5py


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dimension = sys.argv[1]
Diff = float(sys.argv[2])
micro_por = sys.argv[3]
tau = sys.argv[4]
nx1=int(sys.argv[5])
ny1=int(sys.argv[6])
nz1=int(sys.argv[7])
NPX = int(sys.argv[8])
NPY = int(sys.argv[9])
NPZ = int(sys.argv[10])



micro_por_array = [float(x) for x in  micro_por.split(',')]
tau_array = [float(x) for x in tau.split(',')]

NP=NPX*NPY*NPZ

ipz = rank // (NPX * NPY)  # Get ipz by integer division
remainder = rank % (NPX * NPY)
ipy = remainder // NPX     # Get ipy
ipx = remainder % NPX      # Get ipx

baseZ = nz1 // NPZ
remZ = nz1 % NPZ

if ipz < remZ:
   nzp = baseZ + 1
   startZ = ipz * nzp
else:
   nzp = baseZ
   startZ = remZ * (baseZ + 1) + (ipz - remZ) * baseZ


baseY = ny1 // NPY
remY = ny1 % NPY

if ipy < remY:
   nyp = baseY + 1
   startY = ipy * nyp
else:
   nyp = baseY
   startY = remY * (baseY + 1) + (ipy - remY) * baseY

baseX = nx1 // NPX
remX = nx1 % NPX

if ipx < remX:
   nxp = baseX + 1
   startX = ipx * nxp
else:
   nxp = baseX
   startX = remX * (baseX + 1) + (ipx - remX) * baseX

ncell = nxp*nyp*nzp
if NP>1:
  output_path = "processor"+str(rank)+"/"
else:
  output_path = ""


with h5py.File('Result.h5', 'r') as f:
    dsetEps = f['/eps'] 

    # Read the local chunk
    # Slicing format: [startX : startX + nxp, startY : startY + nyp, startZ : startZ + nzp]
    eps = dsetEps[startX:startX+nxp, startY:startY+nyp, startZ:startZ+nzp] 

    if (ipz>0):
      eps_z0 = dsetEps[startX:startX+nxp, startY:startY+nyp, startZ-1]
    if (ipy>0):
      eps_y0 = dsetEps[startX:startX+nxp, startY-1, startZ:startZ+nzp]
    if (ipx>0):
      eps_x0 = dsetEps[startX-1, startY:startY+nyp, startZ:startZ+nzp]
    if (ipx<NPX-1):
      eps_x1 = dsetEps[startX+nxp, startY:startY+nyp, startZ:startZ+nzp]
    if (ipy<NPY-1):
      eps_y1 = dsetEps[startX:startX+nxp, startY+nyp, startZ:startZ+nzp]
    if (ipz<NPZ-1):
      eps_z1 = dsetEps[startX:startX+nxp, startY:startY+nyp, startZ+nzp]

D = np.zeros((nxp,nyp,nzp),dtype="float64")
 
for k in range(0,nzp):
    for j in range(0,nyp):
        for i in range(0,nxp):
            epsVal=eps[i,j,k]
            lower = upper = None
            for n in range(len(micro_por_array) - 1):
                if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                    lower = micro_por_array[n]
                    upper = micro_por_array[n + 1]
                    break  # Exit once the bracket is found
            if lower is None or upper is None:
                raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
            delta=(epsVal-lower)/(upper-lower)
            tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
            D[i,j,k] = Diff*tinv

##get subarray for proc i-1
D_bx0 =np.zeros((nyp,nzp),dtype="float64")
##get subarray for proc i+1
D_bx1 =np.zeros((nyp,nzp),dtype="float64")
##get subarray for proc j-1
D_by0 =np.zeros((nxp,nzp),dtype="float64")
##get subarray for proc j+1
D_by1 =np.zeros((nxp,nzp),dtype="float64")
##get subarray for proc k-1
D_bz0 =np.zeros((nxp,nyp),dtype="float64")
##get subarray for proc k-1
D_bz1 =np.zeros((nxp,nyp),dtype="float64")

#proc k to k-1
if (ipz>0):
  for i in range (0,nxp):
      for j in range (0,nyp):
        epsVal=eps_z0[i,j]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_bz0[i,j] = Diff*tinv

#proc j to j-1
if (ipy>0):
  for i in range (0,nxp):
      for k in range (0,nzp):
        epsVal=eps_y0[i,k]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_by0[i,k] = Diff*tinv

#proc i to i-1
if (ipx>0):
  for j in range (0,nyp):
      for k in range (0,nzp):
        epsVal=eps_x0[j,k]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_bx0[j,k] = Diff*tinv

#proc i to i+1
if (ipx<NPX-1):

  for j in range (0,nyp):
      for k in range (0,nzp):
        epsVal=eps_x1[j,k]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_bx1[j,k] = Diff*tinv


#proc j to j+1
if (ipy<NPY-1):
  for i in range (0,nxp):
      for k in range (0,nzp):
        epsVal=eps_y1[i,k]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_by1[i,k] = Diff*tinv

#proc k to k+1
if (ipz<NPZ-1):
  for i in range (0,nxp):
      for j in range (0,nyp):
        epsVal=eps_z1[i,j]
        lower = upper = None
        for n in range(len(micro_por_array) - 1):
            if min(micro_por_array[n],micro_por_array[n+1]) <= epsVal  <= max(micro_por_array[n],micro_por_array[n + 1]):
                lower = micro_por_array[n]
                upper = micro_por_array[n + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{i},{j},{k}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[n]+delta*(tau_array[n+1]-tau_array[n]))
        D_bz1[i,j] = Diff*tinv

data = [
'/*--------------------------------*- C++ -*----------------------------------*\\\\n',
'| =========                 |                                                 |\n',
'| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n',
'|  \\\\    /   O peration     | Version:  2212                                  |\n',
'|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n',
'|    \\\\/     M anipulation  |                                                 |\n',
'\\*---------------------------------------------------------------------------*/\n',
"FoamFile\n",
"{\n",
"    version     2.0;\n",
"    format      ascii;\n",
"    arch        \"LSB;label=32;scalar=64\";\n",
"    class       volTensorField;\n",
"    location    \"0\";\n",
"    object      D;\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"dimensions      [0 2 -1 0 0 0 0];\n",
"\n",
"internalField   nonuniform List<tensor> \n",
f"{ncell}\n",
"(\n"
]

for k in range(0,nzp):
    for j in range(0,nyp):
        for i in range(0,nxp):
            sci_str = format(D[i,j,k], ".7e")
            mantissa, exponent = sci_str.split('e')
            mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
            D_string = mantissa + 'e' + exponent
            data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")

data.extend([
")\n",
";\n",
"\n",
"boundaryField\n",
"{\n"
"    inlet\n"
"    {\n",
"        type            zeroGradient;\n",
"    }\n"
"    outlet\n"
"    {\n",
"        type            zeroGradient;\n",
"    }\n"
"    \"wall_.*\"""\n"
"    {\n",
"        type            zeroGradient;\n",
"    }\n"
])

if (dimension=="2D"):
    data.extend([
    "    frontAndBack\n"
    "    {\n",
    "        type            empty;\n",
    "    }\n"
    ])

#proc k to k-1
if (ipz>0):
  data.extend([
  f"    procBoundary{rank}to{rank - NPX * NPY}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nxp * nyp}\n",
  "(\n"
  ])
  for j in range(nyp):
    for i in range(nxp):
       sci_str = format(D_bz0[i,j], ".7e")
       mantissa, exponent = sci_str.split('e')
       mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
       D_string = mantissa + 'e' + exponent
       data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")

  data.extend([
  ")\n",
  ";\n",
  "    }\n"
  ])

#proc j to j-1
if (ipy>0):
  data.extend([
  f"    procBoundary{rank}to{rank - NPX}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nxp * nzp}\n",
  "(\n"
  ])
  for k in range(nzp):
    for i in range(nxp):
       sci_str = format(D_by0[i,k], ".7e")
       mantissa, exponent = sci_str.split('e')
       mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
       D_string = mantissa + 'e' + exponent
       data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")

  data.extend([
  ")\n",
  ";\n",
  "    }\n"
  ])

#proc i to i-1
if (ipx>0):
  data.extend([
  f"    procBoundary{rank}to{rank - 1}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nyp * nzp}\n",
  "(\n"
  ])
  for k in range(nzp):
    for j in range(nyp):
        sci_str = format(D_bx0[j,k], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")

  data.extend([
  ")\n",
  ";\n",
  "    }\n"
  ])

##proc i to i+1
if (ipx<NPX-1):
  data.extend([
  f"    procBoundary{rank}to{rank + 1}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nyp * nzp}\n",
  "(\n"
  ])
  for k in range(nzp):
    for j in range(nyp):
        sci_str = format(D_bx1[j,k], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")

  data.extend([
  ")\n",
  ";\n",
  "    }\n"
  ])

#proc j to j+1
if (ipy<NPY-1):
  data.extend([
  f"    procBoundary{rank}to{rank + NPX}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nxp * nzp}\n",
  "(\n"
  ])
  for k in range(nzp):
    for i in range(nxp):
      sci_str = format(D_by1[i,k], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
  data.extend([
  ")\n",
  ";\n",
  "    }\n"
  ])

#proc k to k+1
if (ipz<NPZ-1):
  data.extend([
  f"    procBoundary{rank}to{rank + NPX * NPY}\n",
  "    {\n",
  "        type            processor;\n",
  "        value           nonuniform List<tensor> \n",
  f"{nxp * nyp}\n",
  "(\n"
  ])
  for j in range(nyp):
    for i in range(nxp):
       sci_str = format(D_bz1[i,j], ".7e")
       mantissa, exponent = sci_str.split('e')
       mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
       D_string = mantissa + 'e' + exponent
       data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
        
  data.extend([
  ")\n",
  ";\n",
  "    }\n",
  ])

data.extend([
"}\n",
"\n",
"\n",
"// ************************************************************************* //\n"
])

# Write data to file
with open(output_path+'0/D', 'w') as f:
  f.writelines(data)
  f.close()

