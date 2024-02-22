# @Author: lshuns
# @Date:   2023-02-27 14:01:02
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-11 14:14:48

# >>> sizeU
python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_sizeU.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_sizeU.ini --sep_running_log 