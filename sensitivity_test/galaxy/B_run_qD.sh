# @Author: lshuns
# @Date:   2023-02-01 09:29:06
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-13 11:16:48


# >>> qD
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_qD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_qD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_qD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_qD.ini --sep_running_log 

########### ijmeer 50 cores
##### 1
# Elapsed:40:26:51.28,User=6425680.694,System=237124.861,CPU=4575.7%.

########### amsteldiep 50 cores
##### 2 3: Elapsed:11:42:15.18,User=599505.669,System=81444.587,CPU=1616.1%.

##### 4
# Elapsed:6:19:55.02,User=653862.401,System=6664.552,CPU=2897.6%.
##### 6_1
# Elapsed:0:18.79,User=14.247,System=33.505,CPU=254.0%.