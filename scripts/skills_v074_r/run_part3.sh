# @Author: lshuns
# @Date:   2021-12-15 12:01:14
# @Last Modified by:   lshuns
# @Last Modified time: 2021-12-15 13:43:51

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 50 --rng_seed 97 --cosmic_shear -0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 50 --rng_seed 9701 --cosmic_shear 0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 50 --rng_seed 9702 --cosmic_shear 0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 50 --rng_seed 9703 --cosmic_shear -0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 
