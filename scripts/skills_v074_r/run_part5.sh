# @Author: lshuns
# @Date:   2021-12-15 12:04:46
# @Last Modified by:   lshuns
# @Last Modified time: 2021-12-15 13:44:13

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 50 --rng_seed 99 --cosmic_shear -0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 50 --rng_seed 9901 --cosmic_shear 0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 50 --rng_seed 9902 --cosmic_shear 0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 50 --rng_seed 9903 --cosmic_shear -0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 
