# @Author: lshuns
# @Date:   2021-10-22 17:13:20
# @Last Modified by:   lshuns
# @Last Modified time: 2021-10-22 17:14:31

# 0: generate an example configuration file\n\
# 1: simulate images\n\
# 2: swarp images\n\
# 3: detect objects\n\
# 4: measure photometry\n\
# 5: measure photo-z\n\
# 6: measure shapes\n\
# 7: create a combined catalogue\n\

# ## All possible g1,g2 combinations
# g1Range = np.array([-0.0283,+0.0283,+0.0283,-0.0283])
# g2Range = np.array([+0.0283,+0.0283,-0.0283,-0.0283])

python ../modules/Run.py 1 --runTag m000m000 --threads 1 --rng_seed 940120 --cosmic_shear 0 0 -c ./example_AiryPSF_VIS.ini --sep_running_log

