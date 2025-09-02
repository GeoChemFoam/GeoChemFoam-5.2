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
NPX = int(sys.argv[5])
NPY = int(sys.argv[6])
NPZ = int(sys.argv[7])



micro_por_array = [float(x) for x in  micro_por.split(',')]
tau_array = [float(x) for x in tau.split(',')]

NP=NPX*NPY*NPZ

ipz = rank // (NPX * NPY)  # Get ipz by integer division
remainder = rank % (NPX * NPY)
ipy = remainder // NPX     # Get ipy
ipx = remainder % NPX      # Get ipx

if NP>1:
  output_path = "processor"+str(rank)+"/"
else:
  output_path = ""

with open(output_path+"0/eps", 'r') as f:
    lines = f.readlines()

# Find line with 'internalField' and get the line with the number of entries
for i, line in enumerate(lines):
    if 'internalField' in line and 'nonuniform' in line:
        uniform=0
        n = int(lines[i + 1].strip())  # number of scalar values
        eps= np.empty(n, dtype=np.float64)
        start_idx = i + 3  # Skip the '(', go to the data
        for j in range(n):
            eps[j] = float(lines[start_idx + j].strip())
        i=start_idx+n
        break
    elif 'internalField' in line and 'uniform' in line:
        uniform=1
        n=1
        eps= np.empty(n, dtype=np.float64)
        parts = lines[i].replace(';', '').split()
        eps[0] = float(parts[-1].rstrip(';'))
        i +=1
        break

if (ipz>0):
    found_ipz0=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX*NPY}" in line:
            found_ipz0=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bz0=0
                    nbz0 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bz0 = np.empty(nbz0, dtype=np.float64)
                    for j in range(nbz0):
                        eps_bz0[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbz0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bz0=1
                    nbz0=1
                    eps_bz0 = np.empty(nbz0, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bz0[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipz0==1):
            break

if (ipy>0):
    found_ipy0=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX}" in line:
            found_ipy0=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_by0=0
                    nby0 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_by0 = np.empty(nby0, dtype=np.float64)
                    for j in range(nby0):
                        eps_by0[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nby0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_by0=1
                    nby0=1
                    eps_by0 = np.empty(nby0, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_by0[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipy0==1):
            break

if (ipx>0):
    found_ipx0=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-1}" in line:
            found_ipx0=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bx0=0
                    nbx0 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bx0 = np.empty(nbx0, dtype=np.float64)
                    for j in range(nbx0):
                        eps_bx0[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbx0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bx0=1
                    nbx0=1
                    eps_bx0 = np.empty(nbx0, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bx0[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipx0==1):
            break

if (ipx<NPX-1):
    found_ipx1=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+1}" in line:
            found_ipx1=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bx1=0
                    nbx1 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bx1 = np.empty(nbx1, dtype=np.float64)
                    for j in range(nbx1):
                        eps_bx1[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbx1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bx1=1
                    nbx1=1
                    eps_bx1 = np.empty(nbx1, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bx1[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipx1==1):
            break

if (ipy<NPY-1):
    found_ipy1=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+NPX}" in line:
            found_ipy1=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_by1=0
                    nby1 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_by1 = np.empty(nby1, dtype=np.float64)
                    for j in range(nby1):
                        eps_by1[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nby1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_by1=1
                    nby1=1
                    eps_by1 = np.empty(nby1, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_by1[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipy1==1):
            break

if (ipz<NPZ-1):
    found_ipz1=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+NPX*NPY}" in line:
            found_ipz1=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bz1=0
                    nbz1 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bz1 = np.empty(nbz1, dtype=np.float64)
                    for j in range(nbz1):
                        eps_bz1[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbz1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bz1=1
                    nbz1=1
                    eps_bz1 = np.empty(nbz1, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bz1[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipz1==1):
            break
    i += 1


D = np.zeros(n,dtype="float64")
 
for j in range(0,n):
    epsVal=eps[j]
    lower = upper = None
    for k in range(len(micro_por_array) - 1):
        if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
            lower = micro_por_array[k]
            upper = micro_por_array[k + 1]
            break  # Exit once the bracket is found
    delta=(epsVal-lower)/(upper-lower)
    tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
    D[j] = Diff*tinv

#proc k to k-1
if (ipz>0):
    D_bz0 =np.zeros(nbz0,dtype="float64")
    for j in range (0,nbz0):
      epsVal=eps_bz0[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
        if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
            lower = micro_por_array[k]
            upper = micro_por_array[k + 1]
            break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_bz0[j] = Diff*tinv

#proc j to j-1
if (ipy>0):
    D_by0 =np.zeros(nby0,dtype="float64")
    for j in range (0,nby0):
      epsVal=eps_by0[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
          if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
              lower = micro_por_array[k]
              upper = micro_por_array[k + 1]
              break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_by0[j] = Diff*tinv

#proc i to i-1
if (ipx>0):
    D_bx0 =np.zeros(nbx0,dtype="float64")
    for j in range (0,nbx0):
      epsVal=eps_bx0[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
          if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
              lower = micro_por_array[k]
              upper = micro_por_array[k + 1]
              break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_bx0[j] = Diff*tinv

#proc i to i+1
if (ipx<NPX-1):
    D_bx1 =np.zeros(nbx1,dtype="float64")
    for j in range (0,nbx1):
      epsVal=eps_bx1[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
          if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
              lower = micro_por_array[k]
              upper = micro_por_array[k + 1]
              break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_bx1[j] = Diff*tinv

#proc j to j+1
if (ipy<NPY-1):
    D_by1 =np.zeros(nby1,dtype="float64")
    for j in range (0,nby1):
      epsVal=eps_by1[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
          if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
              lower = micro_por_array[k]
              upper = micro_por_array[k + 1]
              break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_by1[j] = Diff*tinv

#proc k to k+1
if (ipz<NPZ-1):
    D_bz1 =np.zeros(nbz1,dtype="float64")
    for j in range (0,nbz1):
      epsVal=eps_bz1[j]
      lower = upper = None
      for k in range(len(micro_por_array) - 1):
          if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
              lower = micro_por_array[k]
              upper = micro_por_array[k + 1]
              break  # Exit once the bracket is found
      delta=(epsVal-lower)/(upper-lower)
      tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
      D_bz1[j] = Diff*tinv

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
]
if (uniform):
    data.append("internalField uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
else:
    data.extend([
    "internalField   nonuniform List<tensor> \n",
    f"{n}\n",
    "(\n"
    ])

    for j in range(0,n):
        sci_str = format(D[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
    data.extend([
    ")\n",
    ";\n",
    ])
data.extend([
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

##proc k to k-1
if (ipz>0):
    data.extend([
    f"    procBoundary{rank}to{rank - NPX*NPY}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_bz0):
      sci_str = format(D_bz0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbz0}\n",
      "(\n"
      ])
      for j in range(nbz0):
        sci_str = format(D_bz0[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")

#proc j to j-1
if (ipy>0):
    data.extend([
    f"    procBoundary{rank}to{rank - NPX}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_by0):
      sci_str = format(D_by0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nby0}\n",
      "(\n"
      ])
      for j in range(nby0):
        sci_str = format(D_by0[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")

#proc i to i-1
if (ipx>0):
    data.extend([
    f"    procBoundary{rank}to{rank - 1}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_bx0):
      sci_str = format(D_bx0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbx0}\n",
      "(\n"
      ])
      for j in range(nbx0):
        sci_str = format(D_bx0[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")

##proc i to i+1
if (ipx<NPX-1):
    data.extend([
    f"    procBoundary{rank}to{rank + 1}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_bx1):
      sci_str = format(D_bx1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbx1}\n",
      "(\n"
      ])
      for j in range(nbx1):
        sci_str = format(D_bx1[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")


#proc j to j+1
if (ipy<NPY-1):
    data.extend([
    f"    procBoundary{rank}to{rank + NPX}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_by1):
      sci_str = format(D_by1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nby1}\n",
      "(\n"
      ])
      for j in range(nby1):
        sci_str = format(D_by1[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")


#proc k to k+1
if (ipz<NPZ-1):
    data.extend([
    f"    procBoundary{rank}to{rank + NPX * NPY}\n",
    "    {\n",
    "        type            processor;\n",
    ])
    if (uniform_bz1):
      sci_str = format(D_bz1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      D_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({D_string} 0 0 0 {D_string} 0 0 0 {D_string});\n")
    else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbz1}\n",
      "(\n"
      ])
      for j in range(nbz1):
        sci_str = format(D_bz1[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        D_string = mantissa + 'e' + exponent
        data.append(f"({D_string} 0 0 0 {D_string} 0 0 0 {D_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
    data.append("    }\n")


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

