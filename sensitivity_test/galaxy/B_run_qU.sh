# @Author: lshuns
# @Date:   2023-02-04 11:34:56
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-13 12:29:04

# >>> qU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_qU.ini --sep_running_log 

########### amsteldiep 50 cores
#### 1 2 3: Elapsed:46:43:34.24,User=6418893.278,System=257755.726,CPU=3969.1%.
#### 4 5: Elapsed:51:13:42.81,User=6057974.650,System=260297.914,CPU=3425.9%.