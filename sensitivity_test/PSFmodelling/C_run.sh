# @Author: lshuns
# @Date:   2022-07-08 13:58:14
# @Last Modified by:   lshuns
# @Last Modified time: 2023-04-10 12:32:24

python ../../modules/Run.py 7 --runTag m283p283 --threads 48 --rng_seed 94\
 --cosmic_shear -0.0283 0.0283 -c ./C_config.ini --sep_running_log

python ../../modules/Run.py 7 --runTag p283p283 --threads 48 --rng_seed 9401\
 --cosmic_shear 0.0283 0.0283 -c ./C_config.ini --sep_running_log

python ../../modules/Run.py 7 --runTag p283m283 --threads 48 --rng_seed 9402\
 --cosmic_shear 0.0283 -0.0283 -c ./C_config.ini --sep_running_log

python ../../modules/Run.py 7 --runTag m283m283 --threads 48 --rng_seed 9403\
 --cosmic_shear -0.0283 -0.0283 -c ./C_config.ini --sep_running_log


## amsteldiep, order 0 0
# finished in 17.627714981238046 h

# ## amsteldiep order 4 1
# finished in 17.891379081010818 h ======