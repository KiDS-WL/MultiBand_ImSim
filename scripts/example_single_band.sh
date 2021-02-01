# @Author: lshuns
# @Date:   2021-01-29 17:15:39
# @Last Modified by:   lshuns
# @Last Modified time: 2021-02-01 11:54:11

# usage: Run.py [-h] [--runTag RUNTAG] [-c config_file] [--threads number_threads]
#               [--loglevel logging_level] [--sep_running_log] [--version]
#               taskIDs [taskIDs ...]
# positional arguments:
#   taskIDs               Select a set of IDs for processes:
#                             0: generate an example configuration file
#                             1: simulate images
#                             2: swarp images
#                             3: detect objects
#                             4: measure photometry
#                             5: measure photo-z
#                             6: measure shapes
#                             7: create a combined catalogue
#                             all: run all tasks in a sequence
# optional arguments:
#   -h, --help            show this help message and exit
#   --runTag RUNTAG       A tag to label the current run.
#   -c config_file, --config config_file
#                         A configuration file including all necessary setups.
#                            To generate an example file, use `python Run.py 0`.
#   --threads number_threads
#                         The maximum number of threads to use. 
#                            (default: 12 [half of the total CPU count])
#   --loglevel logging_level
#                         Logging level: DEBUG, INFO, WARNING, ERROR, CRITICAL
#   --sep_running_log     If set, save log files for external code running info.
#   --version             The pipeline version.


#################### an example to simulate image and detect objects
### example1: KiDS-like swarp pattern
###      NOTE: 1 for simulation, 2 for SWarp, 3 for detection
###         complete available tasks see above, BUT some are not tested yet for v0.2
# python ../modules/Run.py 1 2 3 --runTag test_KiDS_like --threads 1 -c ./example_single_band.ini --sep_running_log
### example2: only simple stacked image without observation patterns
python ../modules/Run.py 1 3 --runTag test_simple --threads 1 -c ./example_single_band.ini --sep_running_log
