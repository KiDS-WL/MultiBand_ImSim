# @Author: lshuns
# @Date:   2022-07-23 12:49:39
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-01 10:19:16

# >>> sizeU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_sizeU.ini --sep_running_log 

# >>> nU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_nU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_nU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_nU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_nU.ini --sep_running_log 

# >>> qU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_qU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_qU.ini --sep_running_log 

# #### amsteldiep 25 cores
#  Elapsed:51:20:52.46,User=3959454.070,System=90291.769,CPU=2190.7%.

# #### 1236_1 amsteldiep 50 cores
# Elapsed:31:30:23.48,User=4633108.049,System=124955.371,CPU=4194.9%.

# ##### 1236_1 ijmeer 50cores
# Elapsed:39:43:27.24,User=5288796.579,System=174602.163,CPU=3820.3%.