# @Author: lshuns
# @Date:   2023-02-22 18:31:12
# @Last Modified by:   lshuns
# @Last Modified time: 2023-05-18 18:11:46

# # >>> sizeD
python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_sizeD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_sizeD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_sizeD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_sizeD.ini --sep_running_log 

########### ijmeer 50 cores
#### 1: Elapsed:40:14:29.04,User=6455533.494,System=238520.814,CPU=4620.7%.

#### 4 5: Elapsed:50:50:18.01,User=6034016.301,System=257174.909,CPU=3437.4%.