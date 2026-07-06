# Condensate-quant

Code collection for detecting, segmenting, and quantifying biomolecular condensates in live-cell super-resolution microscopy images. The repository includes workflows for multi-channel chromatic alignment, Noise2Void self-supervised denoising, Richardson-Lucy deconvolution, HMRF-based condensate segmentation and visualization.

## Repository Structure

```text
Condensate-quant/
  README.md
  Function/                              # MATLAB utility functions
    bfmatlab/                            # Bio-Formats MATLAB tools
    CondensateAnalysis/                  # Condensate segmentation and analysis functions
    plotFunction/                        # Plotting utilities
  onnx/                                  # Pre-trained Noise2Void ONNX models
  PSF/                                   # PSF files for different channels and imaging modes
  Step0_channel_alignment/               # Multi-channel chromatic alignment
  Step1_pytorch-noise2void/              # PyTorch Noise2Void training and testing
  Step2_deconvolution-segmentation/      # Denoising, deconvolution, segmentation, and quantification
```

## Workflow Overview

The analysis workflow contains three main steps:

1. `Step0_channel_alignment/`
   - Estimate geometric transformation matrices between fluorescence channels using bead calibration data.
   - Apply chromatic correction to ND2 or batch image datasets through a MATLAB App.
2. `Step1_pytorch-noise2void/`
   - Train a Noise2Void-style self-supervised denoising model using PyTorch.
   - Save PyTorch checkpoints and export ONNX models for MATLAB inference.
3. `Step2_deconvolution-segmentation/`
   - Read OME-TIFF time-series image stacks.
   - Apply sliding-window denoising using a pre-trained ONNX Noise2Void model.
   - Perform Richardson-Lucy deconvolution using experimentally measured PSFs.
   - Generate nucleus masks and identify condensate regions using HMRF segmentation.
   - Export condensate masks, centers, boundaries, interfaces, and time-resolved morphometric measurements.

## Requirements

### MATLAB

A recent MATLAB version is recommended, together with the following toolboxes:

- Image Processing Toolbox
- Deep Learning Toolbox
- Statistics and Machine Learning Toolbox
- Computer Vision Toolbox, for some registration and geometric transformation functions

Bio-Formats MATLAB tools are included in:

```text
Function/bfmatlab/
```

Before running MATLAB scripts, add the repository and all subfolders to the MATLAB path:

```matlab
addpath(genpath('C:\path\to\Condensate-quant'));
```

### Python

Noise2Void training and testing scripts are located in `Step1_pytorch-noise2void/`. A conda environment is recommended:

```bash
cd Step1_pytorch-noise2void
conda create -n n2v python=3.10 -y
conda activate n2v
pip install -r requirements.txt
```

On Windows, set the following option in both `train_pipeline.py` and `test_pipeline.py`:

```python
num_workers = 0
```

This avoids common multiprocessing issues during data loading.

## Data

Example and raw data information is listed in `DataSource.txt` files under the corresponding analysis folders, for example:

```text
Step2_deconvolution-segmentation/DataSource.txt
```

The current data source points to the Zenodo DOI:

```text
10.5281/zenodo.18548898
```

`Step2_deconvolution-segmentation/` contain example inputs and generated outputs, which can be used to inspect expected file naming and output organization.

## Step 0: Multi-Channel Chromatic Alignment

Directory:

```text
Step0_channel_alignment/
```

Main files:

- `script_get_transformation_matrix.m`: computes channel-to-channel transformation matrices from fluorescent bead localization results.
- `ImageRegistrator.mlapp`: MATLAB App for applying transformation matrices to biological image data.
- `tform_20260123.mat`: example or precomputed transformation matrix file.

### Input

The registration script expects bead localization CSV files exported from Fiji TrackMate. Each calibration dataset should usually contain:

```text
*_405channel.csv
*_488channel.csv
*_560channel.csv
*_640channel.csv
```

### Usage

1. Modify `filename_list` in `script_get_transformation_matrix.m` to point to your bead calibration datasets.
2. Run the script to generate 2D and 3D geometric transformation matrices.
3. Open `ImageRegistrator.mlapp` and select the input images, output directory, and transformation matrix file.
4. Choose the imaging mode, either Live-SR or Confocal, and export chromatically corrected images.

## Step 1: Noise2Void Denoising Model Training and Testing

Directory:

```text
Step1_pytorch-noise2void/
```

Main files:

- `train_pipeline.py`: trains the Noise2Void model.
- `test_pipeline.py`: runs denoising with a trained model.
- `noise2void/`: core model, data processing, training, and inference code.
- `pth/`: PyTorch checkpoint output directory.
- `onnx/`: ONNX model output directory.
- `requirements.txt`: Python dependencies.

### Training

Modify the dataset path and training parameters in `train_pipeline.py`:

```python
datasets_path = 'your/training/tif/folder'
GPU = '0' # the index of GPU used for computation (e.g. '0', '0,1', '0,1,2')
patch_xy = 128 # the width and height of 3D patches
patch_t = 8 # the time dimension of 3D patches
pth_dir = './pth' # pth file and visualization result file path
```

Then run:

```bash
cd Step1_pytorch-noise2void
python train_pipeline.py
```

Training outputs are saved in `pth/`, and exported ONNX models are saved in `onnx/`.

### Testing

Modify the dataset path, model folder name, and output directory in `test_pipeline.py`:

```python
datasets_path = 'your/input/tif/folder'
denoise_model = 'your_model_folder_name'
output_dir = 'your/output/folder'
```

Then run:

```bash
python test_pipeline.py
```

## Step 2: Denoising, Deconvolution, Condensate Segmentation, and Quantification

Directory:

```text
Step2_deconvolution-segmentation/
```

Main script:

```text
denoise_deconv_condensate_2D_time_series_HMRF_boundary.m
```

### Main Inputs

- OME-TIFF time-series images (`.tif`)
- Noise2Void ONNX model, for example:

```text
onnx/N2V_3D_OCT4_BRD4_liveSR_125_8_E10_xy3z0.onnx
```

- Experimentally measured PSF, for example:

```text
PSF/psf_SR_channel561_2D.tif
```

- Key parameters:
  - `pixelSize`: physical pixel size in nm.
  - `resize_factor`: sub-pixel interpolation factor.
  - `min_radius`: minimum condensate radius threshold in nm.
  - `min_partition_coefficient`: minimum enrichment threshold.
  - `nclust` and `seg_point`: HMRF cluster number and condensate class threshold.

### Parameters to Modify Before Running

Open `denoise_deconv_condensate_2D_time_series_HMRF_boundary.m` and update the following paths and parameters for your data:

```matlab
onnx_path = '../onnx/N2V_3D_OCT4_BRD4_liveSR_125_8_E10_xy3z0.onnx';
psf_SR_560 = TIFFreader('../PSF/psf_SR_channel561_2D.tif', 'double');
filepath_list = {'path/to/your/tif/folder/'};
pixelSize = 95; % Physical pixel size in nanometers (nm)
resize_factor = 10; % Sub-pixel interpolation factor for morphological precision
min_radius = 50; % minimal condensate radius (nm)
min_partition_coefficient = 1.25;   % minimal condensate enrichment level
seg_point = 5;% Threshold for condensate class, above this will be considered as condensates
```

Then run the script in MATLAB.

### Outputs

For each input `.tif` file, the script creates an output folder with the same base name and generates:

```text
*-denoised.tif          # Noise2Void denoising result
*-denoised-deconv.tif   # Denoised and deconvolved result
*-CDmask.tif            # Condensate mask with boundary visualization
*-CDcenter.tif          # Condensate centers
*-CDinterface.tif       # Condensate interfaces
*-CDboundary.tif        # Condensate boundaries
*-HMRFseg.png           # HMRF segmentation preview for the first frame
*.mat                   # Image data, masks, condensate structures, and quantification results
```

The `condensate_result` structure contains:

```text
CDcenter
CDinterface
CDboundary
CDmask
labels
area
```

The second half of the main MATLAB script includes statistical analysis and visualization of condensate number, area, and radius over time. If preprocessing results are already available as `.mat` files, the statistics section can be run independently without repeating denoising, deconvolution, and segmentation.

## Notes

- Several scripts contain example relative or absolute paths. Update them to match your local data paths before running.
- On Windows, set `num_workers = 0` in the Python scripts.
- If GPU inference is unavailable in MATLAB, set:

```matlab
is_GPU_avaliable = false;
```

- Choose a PSF file that matches both the imaging mode and fluorescence channel, such as `SR/noSR` and `405/488/561/642`.
- HMRF segmentation can be sensitive to `nclust`, `seg_point`, `min_radius`, and `min_partition_coefficient`. It is recommended to inspect the output images on a small subset of data before batch processing.

## Citation

For questions, contact:

```text
wdeng@pku.edu.cn
```
