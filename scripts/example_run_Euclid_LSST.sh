# @Author: lshuns
# @Date:   2021-09-27 09:50:44
# @Last Modified by:   lshuns
# @Last Modified time: 2025-11-06 12:52:29

# positional arguments:
#   taskIDs               Select a set of IDs for processes:
#                             0: generate an example configuration file
#                             1: simulate images
#                             2: swarp images
#                             3: detect objects
#                             4: measure photometry
#                             5: measure photo-z
#                             6_1: psf modelling
#                             6_2: measure shapes
#                             7: create a combined catalogue
#                             all: run all tasks in a sequence

python ../modules/Run.py 1 --runTag m000m000 --threads 7 --rng_seed 940120 --cosmic_shear 0 0 -c ./example_Euclid_LSST.ini --sep_running_log

