# @Author: lshuns
# @Date:   2023-02-13 15:04:20
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-31 10:54:53

# >>> nD
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_nD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_nD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_nD.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_nD.ini --sep_running_log 

########### amsteldiep, 50 cores
#### 1: Elapsed:43:16:25.86,User=7184630.817,System=190768.419,CPU=4734.3%.
#### 2 & 3: Elapsed:11:28:05.16,User=603567.381,System=80870.365,CPU=1657.8%.
#### 4 & 5: Elapsed:51:01:14.48,User=6055162.740,System=259486.722,CPU=3437.9%.