# HPC overview

## Important notes

* The `module purge` gets rid of any default compilers.  However, it takes out the queueing system and a few other
  things as collateral damage!!!
* The queueing system is restored later via the appropriate `module load`
* All of the below works once you attend to the code review on your PR.  I made the changes manually/locally, but you
  should confirm that they work against 0.2.1 for yourself!

## Compiling

The important point is that you **must not** use the HPC compiler, as it is "whatever came with the login node's version
of CentOS".  Rather, use the version with the conda env.

Setup:
```sh
# Remove all existing modules
module purge
# Load the lab conda module
module load krthornt/anaconda
# Source the env
source activate fwdpy11_0_2_0
```

Neovim is available, and its defaults are nicer for python editing:

```sh
module load krthornt/neovim
```

Compiling:

```sh
git clone https://github.com/DL42/fwdpy11_arg_example
cd fwdpy11_arg_example
git checkout ancestral_sample
# Not sure why the CPPFLAGS is needed, but it is for some reason:
CPPFLAGS="-I/data/apps/user_contributed_software/krthornt/anaconda3/envs/fwdpy11_0_2_0/lib/python3.6/site-packages/fwdpy11/headers -I/data/apps/user_contributed_software/krthornt/anaconda3/envs/fwdpy11_0_2_0/lib/python3.6/site-packages/fwdpy11/headers/fwdpp" CC=gcc CXX=g++ python3 setup.py build_ext -i
```

The following type of script should be in `/share/kevin2/$USER/project`, where `project` is something sensible, with
well-thought out subdirectories, etc..

```sh
#!/bin/bash

#$ -q krt2
#$ -pe openmp 128

# makes sure to write to output file in directory from which program was launched
cd $SGE_O_WORKDIR 

module load krthornt/anaconda
source activate fwdpy11_0_2_0

# You don't have permission to install to the lab 
# conda module, so you set PYTHONPATH to use your 
# module
PYTHONPATH=$HOME/src/fwdpy11_arg_example python script.py
```

The above script requests the very fast **krt2** queue and requests the entire node for each job.  `concurrent.futures`
will automagically attempt to fill all 128 cores.  The nodes have 512GB RAM.  Thus, if your jobs need > 3.5-4GB each,
set `max_workers` in `concurrent_futures` to something like 64. **This is important!!!**

```#$ -pe openmp 64-128
```
Then, pass $CORES to your python script, and use the value of CORES to set max_workers for `concurrent.futures`.

The contents of `script.py` are in the **same** folder and contain:

```py
import fwdpy11_arg_example

print(fwdpy11_arg_example.__file__)
```

The output I got was:

```
/data/users/krthornt/src/fwdpy11_arg_example/fwdpy11_arg_example/__init__.py
```

## Running jobs

At this point, you should log out/log in to HPC again, **or**:

```sh
# Bring the queuing module back, which you purged early on
module load GE
```

Submit your job:

```
qsub jobscript.sh
```

### Getting the work done:

* 1 job script = 1 param combo
* submit them all.  The concurrent.futures will manage the feeding of the 128 cores
* There are 14 128-core nodes.  USE THEM. The throughput per node is about 5x the dev server.  I imagine you'll have
  most results done by Friday.
* Rsync output back to dev server where there is a fully functional/up-to-date R env for plotting, etc.

### Useful test script

#!/bin/bash
​
​#$ -q krt2
​
​ls -lhrt $HOME > $HOME"_files.txt"
