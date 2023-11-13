# @Author: lshuns
# @Date:   2022-07-03 13:44:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-11-13 10:35:21

#### all nine bands 

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

# >>>> part 0
python ../../modules/Run.py 5 7 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./config_part0_i2_photoz.ini --sep_running_log 

python ../../modules/Run.py 5 7 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./config_part0_i2_photoz.ini --sep_running_log 

python ../../modules/Run.py 5 7 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part0_i2_photoz.ini --sep_running_log 

python ../../modules/Run.py 5 7 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part0_i2_photoz.ini --sep_running_log 

# # >>>> part 1
# python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 95\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part1_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9501\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part1_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9502\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part1_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9503\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part1_i2_photoz.ini --sep_running_log 

# # >>>> part 2
# python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 96\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part2_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9601\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part2_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9602\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part2_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9603\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part2_i2_photoz.ini --sep_running_log 

# # >>>> part 3
# python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 97\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part3_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9701\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part3_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9702\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part3_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9703\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part3_i2_photoz.ini --sep_running_log 

# # >>>> part 4
# python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 98\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part4_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9801\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part4_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9802\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part4_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9803\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part4_i2_photoz.ini --sep_running_log 

# # >>>> part 5
# python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 99\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part5_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9901\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part5_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9902\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part5_i2_photoz.ini --sep_running_log 

# python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9903\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part5_i2_photoz.ini --sep_running_log 


########### amsteldiep, 5 (5 parts)
# Elapsed:128:27:59.50,User=16787311.142,System=1460552.519,CPU=3945.6%.
########### ijmeer, 7 
# Elapsed:2:15:26.23,User=4369.145,System=4496.967,CPU=109.1%.