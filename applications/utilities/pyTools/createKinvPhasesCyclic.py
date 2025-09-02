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
micro_por = sys.argv[2]
micro_k = sys.argv[3]
NPX = int(sys.argv[4])
NPY = int(sys.argv[5])
NPZ = int(sys.argv[6])



micro_por_array = [float(x) for x in  micro_por.split(',')]
micro_k_array = [float(x) for x in micro_k.split(',')]

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

if (ipz==NPZ-1) and (NPZ>2):
    found_ipz01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback" in line:
            found_ipz01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bz01=0
                    nbz01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bz01 = np.empty(nbz01, dtype=np.float64)
                    for j in range(nbz01):
                        eps_bz01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbz01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bz01=1
                    nbz01=1
                    eps_bz01 = np.empty(nbz01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bz01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipz01==1):
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

if (ipz==NPZ-1) and (NPZ==2):
    found_ipz01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback" in line:
            found_ipz01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bz01=0
                    nbz01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bz01 = np.empty(nbz01, dtype=np.float64)
                    for j in range(nbz01):
                        eps_bz01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbz01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bz01=1
                    nbz01=1
                    eps_bz01 = np.empty(nbz01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bz01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipz01==1):
            break

if (ipy==NPY-1) and (NPY>2):
    found_ipy01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop" in line:
            found_ipy01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_by01=0
                    nby01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_by01 = np.empty(nby01, dtype=np.float64)
                    for j in range(nby01):
                        eps_by01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nby01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_by01=1
                    nby01=1
                    eps_by01 = np.empty(nby01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_by01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipy01==1):
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

if (ipy==NPY-1) and (NPY==2):
    found_ipy01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop" in line:
            found_ipy01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_by01=0
                    nby01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_by01 = np.empty(nby01, dtype=np.float64)
                    for j in range(nby01):
                        eps_by01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nby01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_by01=1
                    nby01=1
                    eps_by01 = np.empty(nby01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_by01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipy01==1):
            break

if (ipx==NPX-1) and (NPX>2):
    found_ipx01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX+1}throughright" in line:
            found_ipx01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bx01=0
                    nbx01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bx01 = np.empty(nbx01, dtype=np.float64)
                    for j in range(nbx01):
                        eps_bx01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbx01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bx01=1
                    nbx01=1
                    eps_bx01 = np.empty(nbx01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bx01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipx01==1):
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

if (ipx==NPX-1) and (NPX==2):
    found_ipx01=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank-NPX+1}throughright" in line:
            found_ipx01=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bx01=0
                    nbx01 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bx01 = np.empty(nbx01, dtype=np.float64)
                    for j in range(nbx01):
                        eps_bx01[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbx01  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bx01=1
                    nbx01=1
                    eps_bx01 = np.empty(nbx01, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bx01[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipx01==1):
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


if (ipx==0) and (NPX>1):
    found_ipx10=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+NPX-1}throughleft" in line:
            found_ipx10=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bx10=0
                    nbx10 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bx10 = np.empty(nbx10, dtype=np.float64)
                    for j in range(nbx10):
                        eps_bx10[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbx10  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bx10=1
                    nbx10=1
                    eps_bx10 = np.empty(nbx10, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bx10[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipx10==1):
            break
    i += 1

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

if (ipy==0) and (NPY>1):
    found_ipy10=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+NPX*(NPY-1)}throughbottom" in line:
            found_ipy10=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_by10=0
                    nby10 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_by10 = np.empty(nby10, dtype=np.float64)
                    for j in range(nby10):
                        eps_by10[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nby10  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_by10=1
                    nby10=1
                    eps_by10 = np.empty(nby10, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_by10[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipy10==1):
            break
    i += 1

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

if (ipz==0) and (NPZ>1):
    found_ipz10=0
    while (i<len(lines)):
        line=lines[i].strip()
        if f"procBoundary{rank}to{rank+NPX*NPY*(NPZ-1)}throughfront" in line:
            found_ipz10=1
            # Search for internalField inside this patch
            while (i<len(lines)):
                if 'value' in lines[i] and 'nonuniform' in lines[i]:
                    uniform_bz10=0
                    nbz10 = int(lines[i + 1].strip())
                    start_idx = i + 3
                    eps_bz10 = np.empty(nbz10, dtype=np.float64)
                    for j in range(nbz10):
                        eps_bz10[j] = float(lines[start_idx + j].strip())
                    i = start_idx + nbz10  # skip ahead
                    break
                elif 'value' in lines[i] and 'uniform' in lines[i]:
                    uniform_bz10=1
                    nbz10=1
                    eps_bz10 = np.empty(nbz10, dtype=np.float64)
                    parts = lines[i].replace(';', '').split()
                    eps_bz10[0] = float(parts[-1].rstrip(';'))
                    break
                i += 1
        i +=1
        if (found_ipz10==1):
            break
    i += 1

Kinv = np.zeros(n,dtype="float64")
for j in range(n):
    epsVal=eps[j]
    lower = upper = None
    for k in range(len(micro_por_array) - 1):
        if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
            lower = micro_por_array[k]
            upper = micro_por_array[k + 1]
            break  # Exit once the bracket is found
    if lower is None or upper is None:
        raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
    delta=(epsVal-lower)/(upper-lower)
    kf1=1/micro_k_array[k]*np.power(lower+1e-13,3)/np.power(1-lower+1e-13,2)
    kf2=1/micro_k_array[k+1]*np.power(upper+1e-13,3)/np.power(1-upper+1e-13,2)
    kf=kf1+delta*(kf2-kf1)
    Kinv[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal+1e-13,3)

if (ipz==NPZ-1) and (NPZ>1):
    Kinv_bz01=np.empty(nbz01, dtype=np.float64)
    for j in range(nbz01):
        epsVal=eps_bz01[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bz01[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)
if (ipz>0):
    Kinv_bz0=np.empty(nbz0, dtype=np.float64)
    for j in range(nbz0):
        epsVal=eps_bz0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bz0[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipy==NPY-1) and (NPY>1):
    Kinv_by01=np.empty(nby01, dtype=np.float64)
    for j in range(nby01):
        epsVal=eps_by01[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_by01[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipy>0):
    Kinv_by0=np.empty(nby0, dtype=np.float64)
    for j in range(nby0):
        epsVal=eps_by0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_by0[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)


if (ipx==NPX-1) and (NPX>1):
    Kinv_bx01=np.empty(nbx01, dtype=np.float64)
    for j in range(nbx01):
        epsVal=eps_bx01[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bx01[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)


if (ipx>0):
    Kinv_bx0=np.empty(nbx0, dtype=np.float64)
    for j in range(nbx0):
        epsVal=eps_bx0[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array processor {rank}")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bx0[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipx<NPX-1):
    Kinv_bx1=np.empty(nbx1, dtype=np.float64)
    for j in range(nbx1):
        epsVal=eps_bx1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bx1[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipx==0) and (NPX>1):
    Kinv_bx10=np.empty(nbx10, dtype=np.float64)
    for j in range(nbx10):
        epsVal=eps_bx10[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bx10[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipy<NPY-1):
    Kinv_by1=np.empty(nby1, dtype=np.float64)
    for j in range(nby1):
        epsVal=eps_by1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_by1[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipy==0) and (NPY>1):
    Kinv_by10=np.empty(nby10, dtype=np.float64)
    for j in range(nby10):
        epsVal=eps_by10[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_by10[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipz<NPZ-1):
    Kinv_bz1=np.empty(nbz1, dtype=np.float64)
    for j in range(nbz1):
        epsVal=eps_bz1[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bz1[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

if (ipz==0) and (NPZ>1):
    Kinv_bz10=np.empty(nbz10, dtype=np.float64)
    for j in range(nbz10):
        epsVal=eps_bz10[j]
        lower = upper = None
        for k in range(len(micro_por_array) - 1):
            if min(micro_por_array[k],micro_por_array[k+1]) <= epsVal  <= max(micro_por_array[k],micro_por_array[k + 1]):
                lower = micro_por_array[k]
                upper = micro_por_array[k + 1]
                break  # Exit once the bracket is found
        if lower is None or upper is None:
            raise TypeError(f"eps[{j}] = {epsVal} is outside the range of micro_por_array")
        delta=(epsVal-lower)/(upper-lower)
        kf1=1/micro_k_array[k]*np.power(lower,3)/np.power(1-lower+1e-13,2)
        kf2=1/micro_k_array[k+1]*np.power(upper,3)/np.power(1-upper+1e-13,2)
        kf=kf1+delta*(kf2-kf1)
        Kinv_bz10[j]=kf*np.power(1-epsVal+1e-13,2)/np.power(epsVal,3)

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
"    class       volScalarField;\n",
"    location    \"0\";\n",
"    object      Kinv;\n",
"}\n",
"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n",
"\n",
"dimensions      [0 -2 0 0 0 0 0];\n",
"\n",
]
if (uniform):
    sci_str = format(Kinv[j], ".7e")
    mantissa, exponent = sci_str.split('e')
    mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
    Kinv_str = mantissa + 'e' + exponent
    data.append("internalField uniform {Kinv_str};\n")
else:
    data.extend([
    "internalField   nonuniform List<scalar> \n",
    f"{n}\n",
    "(\n"
    ])

    for j in range(0,n):
        sci_str = format(Kinv[j], ".7e")
        mantissa, exponent = sci_str.split('e')
        mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
        Kinv_str = mantissa + 'e' + exponent
        data.append(Kinv_str+"\n")
    
    data.extend([
    ")\n",
    ";\n",
    ])
data.extend([
"\n",
"boundaryField\n",
"{\n"
"    left\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
"    right\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
"    bottom\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
"    top\n"
"    {\n",
"        type            cyclic;\n",
"    }\n"
])

if (dimension=="2D"):
    data.extend([
    "    frontAndBack\n"
    "    {\n",
    "        type            empty;\n",
    "    }\n"
    ])

else:
    data.extend([
    "    front\n"
    "    {\n",
    "        type            cyclic;\n",
    "    }\n"
    "    back\n"
    "    {\n",
    "        type            cyclic;\n",
    "    }\n"
    ])

if (ipz==NPZ-1) and (NPZ>2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bz01):
      sci_str = format(Kinv_bz01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbz01}\n",
      "(\n"
      ])
      for j in range(nbz01):
          sci_str = format(Kinv_bz01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

##proc k to k-1
if (ipz>0):
  data.extend([
  f"    procBoundary{rank}to{rank - NPX*NPY}\n",
  "    {\n",
  "        type            processor;\n",
  ])
  if (uniform_bz0):
      sci_str = format(Kinv_bz0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbz0}\n",
      "(\n"
      ])
      for j in range(nbz0):
          sci_str = format(Kinv_bz0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipz==NPZ-1) and (NPZ==2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX*NPY*(NPZ-1)}throughback\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bz01):
      sci_str = format(Kinv_bz01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbz01}\n",
      "(\n"
      ])
      for j in range(nbz01):
          sci_str = format(Kinv_bz01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipy==NPY-1) and (NPY>2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_by01):
      sci_str = format(Kinv_by01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nby01}\n",
      "(\n"
      ])
      for j in range(nby01):
          sci_str = format(Kinv_by01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

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
      sci_str = format(Kinv_by0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nby0}\n",
      "(\n"
      ])
      for j in range(nby0):
          sci_str = format(Kinv_by0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipy==NPY-1) and (NPY==2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX*(NPY-1)}throughtop\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_by01):
      sci_str = format(Kinv_by01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nby01}\n",
      "(\n"
      ])
      for j in range(nby01):
          sci_str = format(Kinv_by01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipx==NPX-1) and (NPX>2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX+1}throughright\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bx01):
      sci_str = format(Kinv_bx01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbx01}\n",
      "(\n"
      ])
      for j in range(nbx01):
          sci_str = format(Kinv_bx01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

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
      sci_str = format(Kinv_bx0[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbx0}\n",
      "(\n"
      ])
      for j in range(nbx0):
          sci_str = format(Kinv_bx0[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipx==NPX-1) and (NPX==2):
  data.extend([
  f"    procBoundary{rank}to{rank-NPX+1}throughright\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bx01):
      sci_str = format(Kinv_bx01[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbx01}\n",
      "(\n"
      ])
      for j in range(nbx01):
          sci_str = format(Kinv_bx01[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

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
      sci_str = format(Kinv_bx1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbx1}\n",
      "(\n"
      ])
      for j in range(nbx1):
          sci_str = format(Kinv_bx1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipx==0) and (NPX>1):
  data.extend([
  f"    procBoundary{rank}to{rank+NPX-1}throughleft\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bx10):
      sci_str = format(Kinv_bx10[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbx10}\n",
      "(\n"
      ])
      for j in range(nbx10):
          sci_str = format(Kinv_bx10[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

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
      sci_str = format(Kinv_by1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nby1}\n",
      "(\n"
      ])
      for j in range(nby1):
          sci_str = format(Kinv_by1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipy==0) and (NPY>1):
  data.extend([
  f"    procBoundary{rank}to{rank+NPX*(NPY-1)}throughbottom\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_by10):
      sci_str = format(Kinv_by10[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nby10}\n",
      "(\n"
      ])
      for j in range(nby10):
          sci_str = format(Kinv_by10[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")

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
      sci_str = format(Kinv_bz1[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbz1}\n",
      "(\n"
      ])
      for j in range(nbz1):
          sci_str = format(Kinv_bz1[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

if (ipz==0) and (NPZ>1):
  data.extend([
  f"    procBoundary{rank}to{rank+NPX*NPY*(NPZ-1)}throughfront\n",
  "    {\n",
  "        type            processorCyclic;\n",
  ])
  if (uniform_bz10):
      sci_str = format(Kinv_bz10[0], ".7e")
      mantissa, exponent = sci_str.split('e')
      mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
      Kinv_str = mantissa + 'e' + exponent
      data.append(f"value uniform {Kinv_str};\n")
  else:
      data.extend([
      "        value           nonuniform List<scalar> \n",
      f"{nbz10}\n",
      "(\n"
      ])
      for j in range(nbz10):
          sci_str = format(Kinv_bz10[j], ".7e")
          mantissa, exponent = sci_str.split('e')
          mantissa = mantissa.rstrip('0').rstrip('.')  # Also remove trailing dot if all decimals are stripped
          Kinv_str = mantissa + 'e' + exponent
          data.append(Kinv_str+"\n")
      data.extend([
      ")\n",
      ";\n",
      ])
  data.append("    }\n")

data.extend([
"}\n",
"\n",
"\n",
"// ************************************************************************* //"
])

# Write data to file
with open(output_path+'0/Kinv', 'w') as f:
  f.writelines(data)
  f.close()

