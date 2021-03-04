# @Author: lshuns
# @Date:   2021-01-29 17:15:39
# @Last modified by:   lshuns
# @Last modified time: 2021-03-04, 14:22:12
#
# usage: Run.py [-h] [--runTag RUNTAG] [--cosmic_shear COSMIC_SHEAR COSMIC_SHEAR]
#               [-c config_file] [--threads number_threads] [--loglevel logging_level]
#               [--sep_running_log] [--version]
#               taskIDs [taskIDs ...]
#
# MultiBand_ImSim v0.3: generate realistic multi-band sky images using galaxy & star mock catalogues.
#
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
#
# optional arguments:
#   -h, --help            show this help message and exit
#   --runTag RUNTAG       A tag to label the current run.
#   --cosmic_shear COSMIC_SHEAR COSMIC_SHEAR
#                         2D cosmic shear values: g1 g2.
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

#################### example scripts
### example1: only simple stacked image without observation patterns
###      NOTE: if you want a segmentation map from SEXtractor, set checkimage_type = SEGMENTATION in the configuration file (example_single_band.ini)
python ../modules/Run.py 1 3 7 --runTag test_simple --threads 1 --cosmic_shear 0.0 0.0 -c ./example_single_band.ini --sep_running_log

### example2: KiDS-like observation patterns
###     remember to adjust corresponding variables (survey = KiDS & image_type_list = sameExpo (or diffExpo)) in the configuration file (example_single_band.ini)
# python ../modules/Run.py 1 2 3 7 --runTag test_KiDS_like --threads 1 --cosmic_shear 0.0 0.0 -c ./example_single_band.ini --sep_running_log


### example3: KiDS-VIKING-like 9-band images and all the way to the shear catalogue
###     NOTE: lensfit is slow
# python ../modules/Run.py all --runTag test_9bands --threads 1 --cosmic_shear 0.0 0.0 -c ./example_nice_bands.ini --sep_running_log
