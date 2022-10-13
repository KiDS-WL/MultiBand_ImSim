# @Author: lshuns
# @Date:   2022-07-03 13:44:20
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-13 15:45:07

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
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part0.ini --sep_running_log 

# >>>> part 1
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 95\
 --cosmic_shear -0.0283 0.0283 -c ./config_part1.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9501\
 --cosmic_shear 0.0283 0.0283 -c ./config_part1.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9502\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part1.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9503\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part1.ini --sep_running_log 

# >>>> part 2
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 96\
 --cosmic_shear -0.0283 0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9601\
 --cosmic_shear 0.0283 0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9602\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9603\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part2.ini --sep_running_log 

# >>>> part 3
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 96_1\
 --cosmic_shear -0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 96_101\
 --cosmic_shear 0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 96_102\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 96_103\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 

# >>>> part 4
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 98\
 --cosmic_shear -0.0283 0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9801\
 --cosmic_shear 0.0283 0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9802\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9803\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part4.ini --sep_running_log 

# >>>> part 5
python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 99\
 --cosmic_shear -0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9901\
 --cosmic_shear 0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9902\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9903\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 

