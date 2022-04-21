# @Author: lshuns
# @Date:   2022-01-06 09:01:25
# @Last Modified by:   lshuns
# @Last Modified time: 2022-04-21 15:23:54

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 54 --rng_seed 99\
 --cosmic_shear -0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 54 --rng_seed 9901\
 --cosmic_shear 0.0283 0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 54 --rng_seed 9902\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 54 --rng_seed 9903\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part5.ini --sep_running_log 
