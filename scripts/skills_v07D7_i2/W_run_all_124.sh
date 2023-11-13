# @Author: lshuns
# @Date:   2022-07-03 13:44:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-11-13 10:34:20

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
python ../../modules/Run.py 1 2 4 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./config_part0_i2.ini --sep_running_log 

python ../../modules/Run.py 1 2 4 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./config_part0_i2.ini --sep_running_log 

python ../../modules/Run.py 1 2 4 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part0_i2.ini --sep_running_log 

python ../../modules/Run.py 1 2 4 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part0_i2.ini --sep_running_log 

# # >>>> part 1
# python ../../modules/Run.py 4 --runTag m283p283 --threads 48 --rng_seed 95\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part1_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283p283 --threads 48 --rng_seed 9501\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part1_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283m283 --threads 48 --rng_seed 9502\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part1_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag m283m283 --threads 48 --rng_seed 9503\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part1_i2.ini --sep_running_log 

# # >>>> part 2
# python ../../modules/Run.py 4 --runTag m283p283 --threads 48 --rng_seed 96\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part2_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283p283 --threads 48 --rng_seed 9601\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part2_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283m283 --threads 48 --rng_seed 9602\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part2_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag m283m283 --threads 48 --rng_seed 9603\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part2_i2.ini --sep_running_log 

# # >>>> part 3
# python ../../modules/Run.py 4 --runTag m283p283 --threads 48 --rng_seed 97\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part3_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283p283 --threads 48 --rng_seed 9701\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part3_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283m283 --threads 48 --rng_seed 9702\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part3_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag m283m283 --threads 48 --rng_seed 9703\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part3_i2.ini --sep_running_log 

# # >>>> part 4
# python ../../modules/Run.py 4 --runTag m283p283 --threads 48 --rng_seed 98\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part4_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283p283 --threads 48 --rng_seed 9801\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part4_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283m283 --threads 48 --rng_seed 9802\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part4_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag m283m283 --threads 48 --rng_seed 9803\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part4_i2.ini --sep_running_log 

# # >>>> part 5
# python ../../modules/Run.py 4 --runTag m283p283 --threads 48 --rng_seed 99\
#  --cosmic_shear -0.0283 0.0283 -c ./config_part5_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283p283 --threads 48 --rng_seed 9901\
#  --cosmic_shear 0.0283 0.0283 -c ./config_part5_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag p283m283 --threads 48 --rng_seed 9902\
#  --cosmic_shear 0.0283 -0.0283 -c ./config_part5_i2.ini --sep_running_log 

# python ../../modules/Run.py 4 --runTag m283m283 --threads 48 --rng_seed 9903\
#  --cosmic_shear -0.0283 -0.0283 -c ./config_part5_i2.ini --sep_running_log 



########### core 48, task 1
# Elapsed:30:36:49.58,User=2459521.339,System=85681.397,CPU=2760.2%.

########## core 48. task 2
# Elapsed:14:28:59.84,User=807223.723,System=97128.042,CPU=1734.4%.

########## core 48. task 4
# Elapsed:20:24:23.52,User=1993500.504,System=32956.879,CPU=2758.4%.