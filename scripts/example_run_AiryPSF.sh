# @Author: lshuns
# @Date:   2021-10-22 17:13:20
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-13 15:26:00

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

python ../modules/Run.py 3 --runTag m000m000 --threads 1 --rng_seed 940120 --cosmic_shear 0 0 -c ./example_AiryPSF_VIS.ini --sep_running_log

