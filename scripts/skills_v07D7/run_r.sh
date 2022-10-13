# @Author: lshuns
# @Date:   2022-06-08 11:28:45
# @Last Modified by:   lshuns
# @Last Modified time: 2022-10-13 15:45:23

### only r band

python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./config_part0_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./config_part0_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part0_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part0_r.ini --sep_running_log 


python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 95\
 --cosmic_shear -0.0283 0.0283 -c ./config_part1_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9501\
 --cosmic_shear 0.0283 0.0283 -c ./config_part1_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9502\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part1_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9503\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part1_r.ini --sep_running_log 


python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 96\
 --cosmic_shear -0.0283 0.0283 -c ./config_part2_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9601\
 --cosmic_shear 0.0283 0.0283 -c ./config_part2_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9602\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part2_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9603\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part2_r.ini --sep_running_log 


python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 97\
 --cosmic_shear -0.0283 0.0283 -c ./config_part3_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9701\
 --cosmic_shear 0.0283 0.0283 -c ./config_part3_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9702\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part3_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9703\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part3_r.ini --sep_running_log 


python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 98\
 --cosmic_shear -0.0283 0.0283 -c ./config_part4_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9801\
 --cosmic_shear 0.0283 0.0283 -c ./config_part4_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9802\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part4_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9803\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part4_r.ini --sep_running_log 


python ../../modules/Run.py 1 --runTag m283p283 --threads 48 --rng_seed 99\
 --cosmic_shear -0.0283 0.0283 -c ./config_part5_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283p283 --threads 48 --rng_seed 9901\
 --cosmic_shear 0.0283 0.0283 -c ./config_part5_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag p283m283 --threads 48 --rng_seed 9902\
 --cosmic_shear 0.0283 -0.0283 -c ./config_part5_r.ini --sep_running_log 

python ../../modules/Run.py 1 --runTag m283m283 --threads 48 --rng_seed 9903\
 --cosmic_shear -0.0283 -0.0283 -c ./config_part5_r.ini --sep_running_log 

