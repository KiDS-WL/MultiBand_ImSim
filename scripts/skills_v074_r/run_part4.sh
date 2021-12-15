# @Author: lshuns
# @Date:   2021-12-15 12:03:15
# @Last Modified by:   lshuns
# @Last Modified time: 2021-12-15 13:44:25

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 50 --rng_seed 98 --cosmic_shear -0.0283 0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 50 --rng_seed 9801 --cosmic_shear 0.0283 0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 50 --rng_seed 9802 --cosmic_shear 0.0283 -0.0283 -c ./config_part4.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 50 --rng_seed 9803 --cosmic_shear -0.0283 -0.0283 -c ./config_part4.ini --sep_running_log 
