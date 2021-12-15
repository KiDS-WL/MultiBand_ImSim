# @Author: lshuns
# @Date:   2021-12-15 09:53:50
# @Last Modified by:   lshuns
# @Last Modified time: 2021-12-15 13:42:33

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> --help info
# usage: Run.py [-h] [--runTag RUNTAG]
#               [--cosmic_shear COSMIC_SHEAR COSMIC_SHEAR]
#               [--shear_columns SHEAR_COLUMNS SHEAR_COLUMNS] [-c config_file]
#               [--threads number_threads] [--rng_seed RNG_SEED]
#               [--loglevel logging_level] [--sep_running_log] [--version]
#               taskIDs [taskIDs ...]
# MultiBand_ImSim v0.6.1: generate realistic multi-band sky images using galaxy & star mock catalogues.
# positional arguments:
#   taskIDs               Select a set of IDs for processes:
#                             0: generate an example configuration file
#                             1: simulate images
#                             2: swarp images
#                             3: detect objects
#                             4: measure photometry
#                             5: measure photo-z
#                             6_1: psf modelling
#                             6_2: measure shapes
#                             7: create a combined catalogue
#                             all: run all tasks in a sequence
# optional arguments:
#   -h, --help            show this help message and exit
#   --runTag RUNTAG       A tag to label the current run.
#   --cosmic_shear COSMIC_SHEAR COSMIC_SHEAR
#                         2D cosmic shear values: g1 g2. 
#                            (ignored if shear_columns is provided)
#   --shear_columns SHEAR_COLUMNS SHEAR_COLUMNS
#                         column names to the input shear. 
#                            (use cosmic_shear instead if need constant shear)
#   -c config_file, --config config_file
#                         A configuration file including all necessary setups.
#                            To generate an example file, use `python Run.py 0`.
#   --threads number_threads
#                         The maximum number of threads to use. 
#                            (default: 40 [half of the total CPU count])
#   --rng_seed RNG_SEED   base seed for all random generation.
#   --loglevel logging_level
#                         Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
#   --sep_running_log     If set, save log files for external code running info.
#   --version             The pipeline version.

# >>> running scripts: one line for one shear
python ../../modules/Run.py 1 2 3 7 --runTag m283p283 --threads 50 --rng_seed 94 --cosmic_shear -0.0283 0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283p283 --threads 50 --rng_seed 9401 --cosmic_shear 0.0283 0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag p283m283 --threads 50 --rng_seed 9402 --cosmic_shear 0.0283 -0.0283 -c ./config_part0.ini --sep_running_log 

python ../../modules/Run.py 1 2 3 7 --runTag m283m283 --threads 50 --rng_seed 9403 --cosmic_shear -0.0283 -0.0283 -c ./config_part0.ini --sep_running_log 
