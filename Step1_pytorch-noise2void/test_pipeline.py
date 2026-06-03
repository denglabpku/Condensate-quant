"""
This file will help to demonstrate pipeline for testing microscopy data using the Noise2Void algorithm.
"""
from noise2void.test_collection import testing_class

# %% Select file(s) to be processed (download if not present)

datasets_path = '/mnt/data1/WangBo-DataCenter/ImageData/20250211_OCT4_BRD4-Halo_liveSR/OCT4_BRD4_liveSR/'
denoise_model = 'OCT4_BRD4_liveSR_202508141924'  # A folder containing pth models to be tested
output_dir = '/mnt/data1/WangBo-DataCenter/ImageData/20250211_OCT4_BRD4-Halo_liveSR/OCT4_BRD4_liveSR/'
# %% First setup some parameters for testing
test_datasize = 100000                # the number of frames to be tested (test all frames if the number exceeds the total number of frames in a .tif file)
GPU = '0'                             # the index of GPU you will use for computation (e.g. '0', '0,1', '0,1,2')
patch_xy = 128                        # the width and height of 3D patches, 128
patch_t = 8                           # the time dimension of 3D patches, original 16
overlap_factor = 0.25                 # the overlap factor between two adjacent patches. 
                                      # Since the receptive field of 3D-Unet is ~90, seamless stitching requires an overlap (patch_xyt*overlap_factor）of at least 90 pixels.
masking_ratio = 0.2
neighbor_radius_xy = 3
neighbor_radius_t = 0
num_workers = 4                       # if you use Windows system, set this to 0.

# %% Setup some parameters for result visualization during testing period (optional)
visualize_images_per_epoch = False  # choose whether to display inference performance after each epoch

test_dict = {
    # dataset dependent parameters
    'patch_x': patch_xy,
    'patch_y': patch_xy,
    'patch_t': patch_t,
    'masking_ratio': masking_ratio,
    'neighbor_radius_xy': neighbor_radius_xy,
    'neighbor_radius_t': neighbor_radius_t,
    'overlap_factor':overlap_factor,
    'scale_factor': 1,                  # the factor for image intensity scaling
    'test_datasize': test_datasize,
    'datasets_path': datasets_path,
    'pth_dir': './pth',                 # pth file root path
    'denoise_model' : denoise_model,
    'output_dir' : output_dir,          # result file root path
    # network related parameters
    'fmap': 16,                         # the number of feature maps
    'GPU': GPU,
    'num_workers': num_workers,
    'visualize_images_per_epoch': visualize_images_per_epoch
}
# %%% Testing preparation
# first we create a testing class object with the specified parameters
tc = testing_class(test_dict)
# start the testing process
tc.run()
