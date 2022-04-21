# @Author: lshuns
# @Date:   2021-12-21 16:59:47
# @Last Modified by:   lshuns
# @Last Modified time: 2022-04-21 15:23:16

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 63 --rng_seed 96\
 --cosmic_shear -0.0283 0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 63 --rng_seed 9601\
 --cosmic_shear 0.0283 0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 63 --rng_seed 9602\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part2.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 63 --rng_seed 9603\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part2.ini --sep_running_log 
