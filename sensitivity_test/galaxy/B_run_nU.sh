# @Author: lshuns
# @Date:   2023-02-13 15:07:16
# @Last Modified by:   lshuns
# @Last Modified time: 2023-03-31 12:11:01

# >>> nU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_nU.ini --sep_running_log

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_nU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_nU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_nU.ini --sep_running_log 

############# ijmeer 50 cores
#### 1: Elapsed:55:08:52.22,User=9054066.507,System=287975.713,CPU=4705.5%.

############# amsteldiep 50 cores
#### 2 3: Elapsed:11:26:50.22,User=603247.010,System=80206.583,CPU=1658.4%. 