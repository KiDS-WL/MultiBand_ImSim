# @Author: lshuns
# @Date:   2021-12-21 17:01:41
# @Last Modified by:   lshuns
# @Last Modified time: 2022-04-21 15:23:27

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 54 --rng_seed 97\
 --cosmic_shear -0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 54 --rng_seed 9701\
 --cosmic_shear 0.0283 0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 54 --rng_seed 9702\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 54 --rng_seed 9703\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part3.ini --sep_running_log 
