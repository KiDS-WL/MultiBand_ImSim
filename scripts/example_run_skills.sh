# @Author: lshuns
# @Date:   2021-06-25 18:42:24
# @Last Modified by:   lshuns
# @Last Modified time: 2021-11-09 13:17:43

# 0: generate an example configuration file\n\
# 1: simulate images\n\
# 2: swarp images\n\
# 3: detect objects\n\
# 4: measure photometry\n\
# 5: measure photo-z\n\
# 6: measure shapes\n\
# 7: create a combined catalogue\n\

python ../modules/Run.py 1 2 3 --runTag p283p283 --threads 48 --rng_seed 94012010 --cosmic_shear 0.0283 0.0283 -c ./example_skills_bandr.ini --sep_running_log

python ../modules/Run.py 1 2 3 --runTag m283m283 --threads 48 --rng_seed 94012030 --cosmic_shear -0.0283 -0.0283 -c ./example_skills_bandr.ini --sep_running_log
