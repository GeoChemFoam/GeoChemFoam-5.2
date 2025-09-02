from mpi4py import MPI
import os
import shutil
import glob
import sys

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()  # Process ID

variable_names = sys.argv[1:]

proc_dir = f"processor{rank}"
zero_dir = os.path.join(proc_dir, "0")

# Find directories like processorX/[1-9]*
match_dirs = [d for d in glob.glob(os.path.join(proc_dir, '[1-9]*')) if os.path.isdir(d)]

for d in match_dirs:
    for variable_name in variable_names:
      src_file = os.path.join(d, variable_name)
      dst_file = os.path.join(zero_dir, variable_name)

      if os.path.exists(src_file):
          # If destination file exists, remove it first (optional, or skip instead)
          if os.path.exists(dst_file):
              os.remove(dst_file)
          shutil.move(src_file, zero_dir)

    # After moving, remove the now-empty directory
    shutil.rmtree(d)

# Synchronize all MPI processes
comm.Barrier()

