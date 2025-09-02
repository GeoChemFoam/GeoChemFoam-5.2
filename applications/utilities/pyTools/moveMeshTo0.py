from mpi4py import MPI
import os
import shutil
import glob
import sys

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Process ID

proc_dir = f"processor{rank}"
constant_dir = os.path.join(proc_dir,'constant')
destination_dir = os.path.join(constant_dir, 'polyMesh')
# Find directories like processorX/[1-9]*
match_dirs = [d for d in glob.glob(os.path.join(proc_dir, '[1-9]*')) if os.path.isdir(d)]

for d in match_dirs:
  src_dir = os.path.join(d, 'polyMesh')
  if os.path.isdir(src_dir):
     for file_name in os.listdir(src_dir):
       src_file = os.path.join(src_dir, file_name)
       dst_file = os.path.join(destination_dir, file_name)

       if os.path.exists(src_file):
         # If destination file exists, remove it first (optional, or skip instead)
         if os.path.exists(dst_file):
             os.remove(dst_file)
         shutil.move(src_file, destination_dir)

  # After moving, remove the now-empty directory
  shutil.rmtree(src_dir)

# Synchronize all MPI processes
comm.Barrier()

