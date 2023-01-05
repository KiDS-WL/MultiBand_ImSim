# @Author: lshuns
# @Date:   2022-08-22 15:30:38
# @Last Modified by:   lshuns
# @Last Modified time: 2022-09-01 10:19:44

# >>> sizeU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_sizeU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_sizeU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_sizeU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_sizeU_2.ini --sep_running_log 

# >>> nU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_nU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_nU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_nU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_nU_2.ini --sep_running_log 

# >>> qU
python ../../modules/Run.py 7 --runTag m283p283 --threads 50 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./B_config_qU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283p283 --threads 50 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./B_config_qU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag p283m283 --threads 50 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./B_config_qU_2.ini --sep_running_log 

python ../../modules/Run.py 7 --runTag m283m283 --threads 50 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./B_config_qU_2.ini --sep_running_log 