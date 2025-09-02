import sys
import os
from mpi4py import MPI
import createprocessordir as createprocessordir
import createcellmesh as createcellmesh
import createpointmesh as createpointmesh
import createfacemesh as createfacemesh
import createfacemeshCyclic as createfacemeshCyclic
import createboundarymesh as createboundarymesh
import createboundarymeshCyclic as createboundarymeshCyclic
import createEpsGrayscale as createEpsGrayscale
import createEpsGrayscaleCyclic as createEpsGrayscaleCyclic
import createEpsPhases as createEpsPhases
import createEpsPhasesCyclic as createEpsPhasesCyclic
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

x_dim = int(sys.argv[1])
y_dim = int(sys.argv[2])
z_dim = int(sys.argv[3])
x_min = int(sys.argv[4])
x_max = int(sys.argv[5])
y_min = int(sys.argv[6])
y_max = int(sys.argv[7])
z_min = int(sys.argv[8])
z_max = int(sys.argv[9])
n_x = int(sys.argv[10])
n_y = int(sys.argv[11])
n_z = int(sys.argv[12])
res = float(sys.argv[13])
Image_name = sys.argv[14]
padWidth = int(sys.argv[15])
pores_value = float(sys.argv[16])
solid_value = float(sys.argv[17])
eps_min = float(sys.argv[18])
dimension = sys.argv[19]
direction = int(sys.argv[20])
cyclic = sys.argv[21]
segmentation = sys.argv[22]
micro_por=sys.argv[23]
phases=sys.argv[24]
NPX = int(sys.argv[25])
NPY = int(sys.argv[26])
NPZ = int(sys.argv[27])

#start_time = time.time()
NP=NPX*NPY*NPZ

if rank==0:
    print("create mesh")

if NP>1:

    output_path = "processor"+str(rank)+"/"
    #if rank == 0:
    #    print(f"create processor folders\n")

    createprocessordir.main(rank)


    #if rank == 0:
    #    print(f"create cells\n")

    createcellmesh.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, padWidth, direction, NPX, NPY, NPZ, rank)

else:

    output_path = ""

#if rank == 0:
    #print(f"create points\n")

createpointmesh.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, res, padWidth, direction, NPX, NPY, NPZ, rank, output_path)

#if rank == 0:
    #print(f"create faces\n")

if (cyclic=='yes'):
   createfacemeshCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, padWidth, direction, NPX, NPY, NPZ, rank, output_path) 
else: 
   createfacemesh.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, padWidth, direction, NPX, NPY, NPZ, rank, output_path)

#if rank == 0:
    #print(f"create boundary\n")

if (cyclic=='yes'):
    createboundarymeshCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, padWidth, dimension, direction, NPX, NPY, NPZ, rank, output_path)
else:
    createboundarymesh.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, padWidth, dimension, direction, NPX, NPY, NPZ, rank, output_path)


if rank == 0:
    print(f"create eps")

if (segmentation=='grayscale'):
  if (cyclic=='yes'):
    createEpsGrayscaleCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pores_value, solid_value, eps_min, dimension,direction, NPX, NPY, NPZ, rank, output_path)
  else:
    createEpsGrayscale.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pores_value, solid_value, eps_min, dimension,direction, NPX, NPY, NPZ, rank, output_path)
elif (segmentation=='phases'):
  if (cyclic=='yes'):
    createEpsPhasesCyclic.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pores_value, solid_value, micro_por, phases, dimension,direction, NPX, NPY, NPZ, rank, output_path)
  else:
    createEpsPhases.main(x_dim, y_dim, z_dim, x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, Image_name, padWidth, pores_value, solid_value, micro_por, phases, dimension,direction, NPX, NPY, NPZ, rank, output_path)
else:
  raise TypeError("only grayscale and phases segmentation accepted")


#end_time = time.time()
#elapsed_time = end_time - start_time

#if rank ==0:
    #print("Elapsed time createMesh.py:", elapsed_time, "seconds")

MPI.Finalize()

