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
Lpore = sys.argv[5]
betax = sys.argv[6]
alpha1x = sys.argv[7]
alpha2x = sys.argv[8]
betay = sys.argv[9]
alpha1y = sys.argv[10]
alpha2y = sys.argv[11]
betaz = sys.argv[12]
alpha1z = sys.argv[13]
alpha2z = sys.argv[14]
NPX = int(sys.argv[15])
NPY = int(sys.argv[16])
NPZ = int(sys.argv[17])

micro_por_array = [float(x) for x in  micro_por.split(',')]
tau_array = [float(x) for x in tau.split(',')]
Lpore_array = [float(x) for x in Lpore.split(',')]
betax_array = [float(x) for x in betax.split(',')]
alpha1x_array = [float(x) for x in alpha1x.split(',')]
alpha2x_array = [float(x) for x in alpha2x.split(',')]
betay_array = [float(x) for x in betay.split(',')]
alpha1y_array = [float(x) for x in alpha1y.split(',')]
alpha2y_array = [float(x) for x in alpha2y.split(',')]
betaz_array = [float(x) for x in betaz.split(',')]
alpha1z_array = [float(x) for x in alpha1z.split(',')]
alpha2z_array = [float(x) for x in alpha2z.split(',')]

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

with open(output_path+"0/U", 'r') as f:
    lines = f.readlines()

# Find line with 'internalField' and get the line with the number of entries
for i, line in enumerate(lines):
    if 'internalField' in line and 'nonuniform' in line:
        n = int(lines[i + 1].strip())  # number of scalar values
        if (uniform==1):
            uniform=0
            epsVal=eps[0]
            eps=np.full(n,epsVal,dtype=np.float64)
        U= np.empty(n, dtype=np.float64)      
        start_idx = i + 3  # Skip the '(', go to the data
        for j in range(n):
            vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
            parts = vec_line.split()
            Ux, Uy, Uz = map(float, parts)
            U[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
        i=start_idx+n
        break
    elif 'internalField' in line and 'uniform' in line:
        parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
        Ux = float(parts[-3])
        Uy = float(parts[-2])
        Uz = float(parts[-1])
        magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
        U=np.full(n,magU,dtype=np.float64)
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
                    nbz0 = int(lines[i + 1].strip())
                    if (uniform_bz0==1):
                        uniform_bz0=0
                        epsVal=eps_bz0[0]
                        eps_bz0=np.full(nbz0,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_bz0 = np.empty(nbz0, dtype=np.float64)
                    for j in range(nbz0):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_bz0[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nbz0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_bz0=np.full(nbz0,magU,dtype=np.float64)
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
                    nby0 = int(lines[i + 1].strip())
                    if (uniform_by0==1):
                        uniform_by0=0
                        epsVal=eps_by0[0]
                        eps_by0=np.full(nby0,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_by0 = np.empty(nby0, dtype=np.float64)
                    for j in range(nby0):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_by0[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nby0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_by0=np.full(nby0,magU,dtype=np.float64)
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
                    nbx0 = int(lines[i + 1].strip())
                    if (uniform_bx0==1):
                        uniform_bx0=0
                        epsVal=eps_bx0[0]
                        eps_bx0=np.full(nbx0,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_bx0 = np.empty(nbx0, dtype=np.float64)
                    for j in range(nbx0):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_bx0[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nbx0  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_bx0=np.full(nbx0,magU,dtype=np.float64)
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
                    nbx1 = int(lines[i + 1].strip())
                    if (uniform_bx1==1):
                        uniform_bx1=0
                        epsVal=eps_bx1[0]
                        eps_bx1=np.full(nbx1,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_bx1 = np.empty(nbx1, dtype=np.float64)
                    for j in range(nbx1):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_bx1[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nbx1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_bx1=np.full(nbx1,magU,dtype=np.float64)
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
                    nby1 = int(lines[i + 1].strip())
                    if (uniform_by1==1):
                        uniform_by1=0
                        epsVal=eps_by1[0]
                        eps_by1=np.full(nby1,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_by1 = np.empty(nby1, dtype=np.float64)
                    for j in range(nby1):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_by1[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nby1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_by1=np.full(nby1,magU,dtype=np.float64)
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
                    nbz1 = int(lines[i + 1].strip())
                    if (uniform_bz1==1):
                        uniform_bz1=0
                        epsVal=eps_bz1[0]
                        eps_bz1=np.full(nbz1,epsVal,dtype=np.float64)
                    start_idx = i + 3
                    U_bz1 = np.empty(nbz1, dtype=np.float64)
                    for j in range(nbz1):
                        vec_line = lines[start_idx + j].strip().strip('()')  # remove surrounding ()
                        parts = vec_line.split()
                        Ux, Uy, Uz = map(float, parts)
                        U_bz1[j]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    i = start_idx + nbz1  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    parts = lines[i].replace(';', '').replace('(', '').replace(')', '').split()
                    Ux = float(parts[-3])
                    Uy = float(parts[-2])
                    Uz = float(parts[-1])
                    magU=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
                    U_bz1=np.full(nbz1,magU,dtype=np.float64)
                    break
                i += 1
        i +=1
        if (found_ipz1==1):
            break
    i += 1

Dx = np.zeros(n,dtype="float64")
Dy = np.zeros(n,dtype="float64")
Dz = np.zeros(n,dtype="float64")
for j in range(n):
    epsVal=eps[j]
    magU=U[j]
    lower = upper = None
    for k in range(len(micro_por_array) - 1):
        if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
            lower = micro_por_array[k]
            upper = micro_por_array[k + 1]
            break  # Exit once the bracket is found
    delta=(epsVal-lower)/(upper-lower)
    tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
    Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
    betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
    betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
    betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
    alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
    alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
    alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
    alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
    alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
    alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

    if (Pe<1):
        Dx[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
        Dy[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
        Dz[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
    else:
        Dx[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
        Dy[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
        Dz[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))
 
if (ipz>0):
    if (found_ipz0==0):
        print(f"rank is {rank}")
    Dx_bz0=np.empty(nbz0, dtype=np.float64)
    Dy_bz0=np.empty(nbz0, dtype=np.float64)
    Dz_bz0=np.empty(nbz0, dtype=np.float64)
    for j in range(nbz0):
        epsVal=eps_bz0[j]
        magU=U_bz0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_bz0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_bz0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_bz0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_bz0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_bz0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_bz0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

if (ipy>0):
    Dx_by0=np.empty(nby0, dtype=np.float64)
    Dy_by0=np.empty(nby0, dtype=np.float64)
    Dz_by0=np.empty(nby0, dtype=np.float64)
    for j in range(nby0):
        epsVal=eps_by0[j]
        magU=U_by0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_by0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_by0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_by0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_by0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_by0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_by0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

if (ipx>0):
    Dx_bx0=np.empty(nbx0, dtype=np.float64)
    Dy_bx0=np.empty(nbx0, dtype=np.float64)
    Dz_bx0=np.empty(nbx0, dtype=np.float64)
    for j in range(nbx0):
        epsVal=eps_bx0[j]
        magU=U_bx0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_bx0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_bx0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_bx0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_bx0[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_bx0[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_bx0[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

if (ipx<NPX-1):
    Dx_bx1=np.empty(nbx1, dtype=np.float64)
    Dy_bx1=np.empty(nbx1, dtype=np.float64)
    Dz_bx1=np.empty(nbx1, dtype=np.float64)
    for j in range(nbx1):
        epsVal=eps_bx1[j]
        magU=U_bx1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_bx1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_bx1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_bx1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_bx1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_bx1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_bx1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

if (ipy<NPY-1):
    if (found_ipy1==0):
        print(f"rank is {rank}")
    Dx_by1=np.empty(nby1, dtype=np.float64)
    Dy_by1=np.empty(nby1, dtype=np.float64)
    Dz_by1=np.empty(nby1, dtype=np.float64)
    for j in range(nby1):
        epsVal=eps_by1[j]
        magU=U_by1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_by1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_by1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_by1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_by1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_by1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_by1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

if (ipz<NPZ-1):
    Dx_bz1=np.empty(nbz1, dtype=np.float64)
    Dy_bz1=np.empty(nbz1, dtype=np.float64)
    Dz_bz1=np.empty(nbz1, dtype=np.float64)
    for j in range(nbz1):
        epsVal=eps_bz1[j]
        magU=U_bz1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        tinv=1.0/(tau_array[k]+delta*(tau_array[k+1]-tau_array[k]))
        Pe=magU/epsVal*(Lpore_array[k]+delta*(Lpore_array[k+1]-Lpore_array[k]))/Diff
        betaxVal=betax_array[k]+delta*(betax_array[k+1]-betax_array[k])
        betayVal=betay_array[k]+delta*(betay_array[k+1]-betay_array[k])
        betazVal=betaz_array[k]+delta*(betaz_array[k+1]-betaz_array[k])
        alpha1xVal=alpha1x_array[k]+delta*(alpha1x_array[k+1]-alpha1x_array[k])
        alpha1yVal=alpha1y_array[k]+delta*(alpha1y_array[k+1]-alpha1y_array[k])
        alpha1zVal=alpha1z_array[k]+delta*(alpha1z_array[k+1]-alpha1z_array[k])
        alpha2xVal=alpha2x_array[k]+delta*(alpha2x_array[k+1]-alpha2x_array[k])
        alpha2yVal=alpha2y_array[k]+delta*(alpha2y_array[k+1]-alpha2y_array[k])
        alpha2zVal=alpha2z_array[k]+delta*(alpha2z_array[k+1]-alpha2z_array[k])

        if (Pe<1):
            Dx_bz1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha1xVal))
            Dy_bz1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha1yVal))
            Dz_bz1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha1zVal))
        else:
            Dx_bz1[j] = Diff*tinv*(1+betaxVal*pow(Pe,alpha2xVal))
            Dy_bz1[j] = Diff*tinv*(1+betayVal*pow(Pe,alpha2yVal))
            Dz_bz1[j] = Diff*tinv*(1+betazVal*pow(Pe,alpha2zVal))

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
    sci_str = format(Dx[0], ".7e")
    mantissa, exponent = sci_str.split('e')
    mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
    Dx_string = mantissa + 'e' + exponent
    sci_str = format(Dy[0], ".7e")
    mantissa, exponent = sci_str.split('e')
    mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
    Dy_string = mantissa + 'e' + exponent
    sci_str = format(Dz[0], ".7e")
    mantissa, exponent = sci_str.split('e')
    mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
    Dz_string = mantissa + 'e' + exponent
    data.append(f"internalField uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
else:
    data.extend([
    "internalField   nonuniform List<tensor> \n",
    f"{n}\n",
    "(\n"
    ])

    for j in range(0,n):
        sci_str = format(Dx[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        Dx_string = mantissa + 'e' + exponent
        sci_str = format(Dy[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        Dy_string = mantissa + 'e' + exponent
        sci_str = format(Dz[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        Dz_string = mantissa + 'e' + exponent
        data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
    
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
      sci_str = format(Dx_bz0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_bz0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_bz0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbz0}\n",
      "(\n"
      ])
      for j in range(nbz0):
          sci_str = format(Dx_bz0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_bz0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_bz0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

##proc j to j-1
if (ipy>0):
  data.extend([
  f"    procBoundary{rank}to{rank - NPX}\n",
  "    {\n",
  "        type            processor;\n",
  ])
  if (uniform_by0):
      sci_str = format(Dx_by0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_by0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_by0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nby0}\n",
      "(\n"
      ])
      for j in range(nby0):
          sci_str = format(Dx_by0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_by0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_by0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")

      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

##proc i to i-1
if (ipx>0):
  data.extend([
  f"    procBoundary{rank}to{rank - 1}\n",
  "    {\n",
  "        type            processor;\n",
  ])
  if (uniform_bx0):
      sci_str = format(Dx_bx0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_bx0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_bx0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbx0}\n",
      "(\n"
      ])
      for j in range(nbx0):
          sci_str = format(Dx_bx0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_bx0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_bx0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
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
      sci_str = format(Dx_bx1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_bx1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_bx1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbx1}\n",
      "(\n"
      ])
      for j in range(nbx1):
          sci_str = format(Dx_bx1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_bx1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_bx1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

##proc j to j+1
if (ipy<NPY-1):
  data.extend([
  f"    procBoundary{rank}to{rank + NPX}\n",
  "    {\n",
  "        type            processor;\n",
  ])
  if (uniform_by1):
      sci_str = format(Dx_by1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_by1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_by1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nby1}\n",
      "(\n"
      ])
      for j in range(nby1):
          sci_str = format(Dx_by1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_by1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_by1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

#proc k to k+1
if (ipz<NPZ-1):
  data.extend([
  f"    procBoundary{rank}to{rank + NPX*NPY}\n",
  "    {\n",
  "        type            processor;\n",
  ])
  if (uniform_bz1):
      sci_str = format(Dx_bz1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dx_string = mantissa + 'e' + exponent
      sci_str = format(Dy_bz1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dy_string = mantissa + 'e' + exponent
      sci_str = format(Dz_bz1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Dz_string = mantissa + 'e' + exponent
      data.append(f"value uniform ({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string});\n")
  else:
      data.extend([
      "        value           nonuniform List<tensor> \n",
      f"{nbz1}\n",
      "(\n"
      ])
      for j in range(nbz1):
          sci_str = format(Dx_bz1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dx_string = mantissa + 'e' + exponent
          sci_str = format(Dy_bz1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dy_string = mantissa + 'e' + exponent
          sci_str = format(Dz_bz1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Dz_string = mantissa + 'e' + exponent
          data.append(f"({Dx_string} 0 0 0 {Dy_string} 0 0 0 {Dz_string})\n")
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

