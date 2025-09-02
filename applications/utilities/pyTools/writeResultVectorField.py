from mpi4py import MPI
import h5py
import numpy as np
import sys

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


nx1 = int(sys.argv[1])
ny1 = int(sys.argv[2])
nz1 = int(sys.argv[3])
NPX = int(sys.argv[4])
NPY = int(sys.argv[5])
NPZ = int(sys.argv[6])

NP=NPX*NPY*NPZ

if NP>1:
  output_path = "processor"+str(rank)+"/"
else:
  output_path = ""

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

ncell=nxp*nyp*nzp

# Each process creates its local U
eps_local = np.zeros((nxp, nyp, nzp), dtype='f8')

file = open(output_path+"0/eps","r")
Lines = file.readlines()
count=0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      k = count // (nxp * nyp)
      remainder = count % (nxp * nyp)
      j = remainder // nxp
      i = remainder % nxp
      eps_local[i,j,k]=ls
      count +=1
  if (ls=="("):
      wbool=1
file.close()

# Gather sizes and data to rank 0
all_eps = comm.gather(eps_local, root=0)
all_startX = comm.gather(startX, root=0)
all_startY = comm.gather(startY, root=0)
all_startZ = comm.gather(startZ, root=0)


# Each process creates its local U
U_local = np.zeros((nxp, nyp, nzp), dtype='f8')

file = open(output_path+"0/U","r")
Lines = file.readlines()
count=0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      k = count // (nxp * nyp)
      remainder = count % (nxp * nyp)
      j = remainder // nxp
      i = remainder % nxp
      Ux=float(ls.split("(")[1].split(")")[0].split()[0])
      Uy=float(ls.split("(")[1].split(")")[0].split()[1])
      Uz=float(ls.split("(")[1].split(")")[0].split()[2])
      U_local[i,j,k]=np.sqrt(Ux*Ux+Uy*Uy+Uz*Uz)
      count +=1
  if (ls=="("):
      wbool=1
file.close()

# Gather sizes and data to rank 0
all_U = comm.gather(U_local, root=0)
all_startX = comm.gather(startX, root=0)
all_nxp    = comm.gather(nxp, root=0)
all_startY = comm.gather(startY, root=0)
all_nyp    = comm.gather(nyp, root=0)
all_startZ = comm.gather(startZ, root=0)
all_nzp    = comm.gather(nzp, root=0)

if rank == 0:
    # Assemble global array
    eps_global = np.empty((nx1, ny1, nz1), dtype='f8')
    U_global = np.empty((nx1, ny1, nz1), dtype='f8')

    for r, block in enumerate(all_eps):
        # Use startX[r], startY[r], startZ[r] to place block in U_global
        eps_global[all_startX[r]:all_startX[r]+all_nxp[r], all_startY[r]:all_startY[r]+all_nyp[r], all_startZ[r]:all_startZ[r]+all_nzp[r]] = block

    for r, block in enumerate(all_U):
        # Use startX[r], startY[r], startZ[r] to place block in U_global
        U_global[all_startX[r]:all_startX[r]+all_nxp[r], all_startY[r]:all_startY[r]+all_nyp[r], all_startZ[r]:all_startZ[r]+all_nzp[r]] = block

    with h5py.File("Result.h5", "w") as f:
        f.create_dataset("eps", data=eps_global)
        f.create_dataset("U", data=U_global)

