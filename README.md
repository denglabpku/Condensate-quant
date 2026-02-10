# Analysis scripts for condensate analysis and burst modeling in biological imaging data
## Overview

This repository contains scrips associated with the paper "**Coordinated dynamics of condensates and enhancer-promoter looping revealed by live-cell imaging**". It aims to support high-throughput quantitative analysis of super-resolution images and time series data, focusing on dynamic processes such as transcriptional bursting and condensate dynamics.

## Fig1. condensate analysis

This MATLAB pipeline processes time-lapse fluorescence microscopy data to quantitatively analyze the morphology and dynamics of OCT4 and BRD4 biomolecular condensates at the single-cell level. The workflow integrates deep-learning–based denoising, optical deconvolution, nucleus segmentation, and HMRF-based condensate identification to extract robust morphometric features for downstream statistical analysis.

#### **Inputs**

- Multi-channel OME-TIFF image stacks (`.tif`)
  - Time-lapse, multi-z, multi-color images
  - Channels include OCT4 / BRD4 (condensate-forming proteins)
- Pre-trained Noise2Void (N2V) deep-learning model (`.onnx`)
  - Trained on live-SR OCT4/BRD4 datasets
- Experimentally measured PSF image (`.tif`)
  - Corresponding to the 561 nm detection channel
- User-defined acquisition parameters
  - Physical pixel size (nm)
  - Sub-pixel interpolation factor for morphometric precision

#### **Main Processing Steps**

#### **1. Data loading and organization**

- Import OME-TIFF images using Bio-Formats.
- Reconstruct the full 5D image stack organized as:
  - `(Y, X, Z, Channel, Time)`
- Automatically iterate over multiple image files and experiments.

#### **2. Deep-learning–based image denoising (Noise2Void)**

- Apply a pre-trained Noise2Void neural network to suppress shot noise while preserving fine condensate structures.
- Perform sliding-window patch inference with overlap to avoid edge artifacts.
- Normalize intensity per patch and restore the original intensity scale after inference.
- GPU acceleration is optionally enabled to improve performance.

#### **3. Richardson–Lucy deconvolution**

- Perform Richardson–Lucy (RL) deconvolution on denoised images using an experimentally measured PSF.
- Padding is applied to minimize boundary artifacts.
- Deconvolution enhances spatial resolution and sharpens condensate boundaries.

#### **4. Nucleus segmentation**

- Generate a global nucleus mask from denoised images using Gaussian filtering and intensity thresholding.
- The nucleus mask restricts downstream condensate detection to nuclear regions.

#### **5. HMRF-based condensate segmentation**

- Apply Hidden Markov Random Field (HMRF) segmentation to determine adaptive, frame-specific intensity thresholds.
- Cluster pixel intensities into multiple classes and identify high-intensity condensate states.
- Use HMRF-derived thresholds to binarize condensate regions within the nucleus.

#### **6. Condensate topology extraction**

- Perform sub-pixel upsampling to improve morphological precision.
- For each frame, extract:
  - Condensate center
  - Condensate interface
  - Condensate boundary
- Label individual condensates and compute their area distributions.

#### **7. Morphometric quantification**

- Calculate condensate physical properties over time, including:
  - Area (nm²)
  - Effective radius (derived from area)
- Aggregate statistics across frames, cells, and experiments.

#### **8. Statistical analysis and visualization**

- Compute ensemble-averaged condensate metrics across datasets.
- Visualize temporal dynamics using mean ± standard deviation curves.

#### **Outputs**

- MATLAB `.mat` files containing:
  - Denoised and deconvolved image data
  - Nucleus masks
  - Condensate segmentation results
  - Condensate area distributions per frame
  - Time-resolved morphometric measurements

## Fig2-4. Multi-channel image analysis

This MATLAB pipeline processes multi-channel 3D fluorescence microscopy data to quantitatively analyze the spatial relationship between DNA/RNA foci and biomolecular condensates (OCT4 and BRD4) at the single-cell level. The workflow integrates deep-learning–based denoising, optical deconvolution, nucleus-restricted condensate segmentation, and sub-pixel distance measurements to extract robust DNA/RNA–condensate spatial metrics.

#### **Inputs**

- Multi-channel OME-TIFF image stacks (`.tif`)
  - Time-lapse, multi-z, multi-color images
  - Channels:
    - Promoter (405 nm)
    - RNA (488 nm)
    - SCR/OCT4 (561 nm)
    - OCT4/BRD4 (640 nm)
- Pre-trained Noise2Void (N2V) deep-learning model (`.onnx`)
  - Optimized for live-SR 2D time-lapse imaging
- Experimentally measured PSF images (`.tif`)
  - Channel-specific PSFs for 405/488 / 561 / 640 nm
- User-defined parameters
  - Physical pixel size (nm)
  - ROI size around DNA/RNA foci
  - Sub-pixel interpolation factor for geometric precision

#### **Main Processing Steps**

#### **1. Data loading and organization**

- Import OME-TIFF stacks using Bio-Formats.
- Reconstruct full image data with dimensions:
  - `(Y, X, Z, Channel, Time)`
- Automatically iterate over multiple datasets and experiments.

#### **2. Deep-learning–based denoising (Noise2Void)**

- Apply a pre-trained 2D Noise2Void neural network to all channels independently.
- Perform sliding-window inference with overlapping patches to suppress edge artifacts.
- Normalize intensity per patch and restore original intensity scaling after inference.
- Optional GPU acceleration is used for efficient processing.

#### **3. Multi-channel Richardson–Lucy deconvolution**

- Apply Richardson–Lucy deconvolution separately to each fluorescence channel using the corresponding experimental PSF.
- Padding is introduced to minimize boundary artifacts.
- Both deconvolved and re-convolved images are retained for quality control and export.

#### **4. Nucleus segmentation**

- Generate a global nucleus mask using the OCT4 channel.
- The nucleus mask is used to constrain downstream foci detection and condensate segmentation to nuclear regions.

#### **5. DNA/RNA foci detection and tracking**

- Detect DNA/RNA foci in 3D using Laplacian-of-Gaussian (LoG) filtering.
- Identify the optimal focal plane based on 3D spot quality metrics.
- Refine 2D DNA/RNA foci positions frame-by-frame within the nucleus.
- Compute RNA intensity normalized to local nuclear background.
- Extract time-resolved DNA/RNA-centered ROIs for visualization and downstream analysis.

#### **6. Condensate segmentation (OCT4 and BRD4)**

- Perform Hidden Markov Random Field (HMRF) segmentation to adaptively classify pixel intensities within the nucleus.
- Identify high-intensity condensate states using channel-specific cluster thresholds.
- For each frame and condensate channel, extract:
  - Condensate center
  - Condensate interface
  - Condensate boundary
  - Condensate mask and labeled regions
- Segmentation is performed at sub-pixel resolution via image upsampling.

#### **7. Sub-pixel ROI extraction**

- Define a fixed-size ROI centered on the DNA/RNA locus.
- Rescale ROIs by a user-defined interpolation factor to enable precise spatial measurements.
- Extract DNA, RNA, OCT4, and BRD4 ROIs consistently across all frames.

#### **8. DNA/RNA–condensate spatial measurements**

For each time point and condensate channel:

- Compute DNA/RNA distance to:
  - Condensate boundary
  - Condensate center
  - Condensate centroid
- Estimate condensate geometric properties:
  - Physical radius
  - Equivalent radius (area-based)
  - Boundary-to-center and boundary-to-centroid distances
- Convert all distances to physical units (nm).

#### **9. Visualization and temporal analysis**

- Generate overlay plots showing:
  - DNA/RNA foci positions
  - Condensate masks and boundaries
- Plot time-resolved RNA intensity and DNA/RNA foci–condensate distance trajectories.
- Compare DNA/RNA foci proximity to OCT4 versus BRD4 condensates within the same nucleus.

#### **Outputs**

- Exported image files (`.tif`)
  - Denoised and deconvolved multi-channel stacks
  - foci-centered ROIs
  - Condensate center, interface, boundary, and mask images
- MATLAB `.mat` files containing:
  - DNA/RNA foci positions, intensities, and background measurements
  - Condensate segmentation results and geometric descriptors
  - Time-resolved DNA/RNA foci–condensate distance measurements
- Publication-ready figures (`.png`)
  - Condensate segmentation visualization
  - DNA/RNA foci–condensate spatial relationship over time

## Fig5. burst modeling

### 1. RNA detection and tracking

This MATLAB script processes time-lapse RNA imaging data to detect and track transcriptional spots over time.

- **Inputs**:
  - A `.tif` image stack with RNA signal
  - An Ilastik-generated probability map (`.h5` file) for segmentation
  - An ROI file (`-ROI.txt`) defining cell regions
- **Main steps**:
  1. Load RNA images and probability maps
  2. Apply thresholding to identify potential RNA spots
  3. Detect the brightest spot in each frame within defined ROIs
  4. Track spot positions over time, fill missing data by interpolation
  5. Measure RNA signal intensity and correct for local background
- **Outputs**:
  - `spots`: all detected spot coordinates and intensities per frame
  - `tracks`: tracked spot trajectories and background-corrected RNA intensity traces

These processed data are saved as `.mat` files, ready for downstream burst modeling and analysis.

### 2. HMM modeling and burst analysis

This MATLAB script analyzes RNA bursting dynamics by applying a Hidden Markov Model (HMM) and fitting burst duration distributions.

- **Inputs**:
  - `.mat` files containing RNA tracks and intensity traces (from Step 1) under different experimental conditions
- **Main steps**:

1. **Merge data**:
    Load and combine RNA tracks from different treatment groups into a single dataset.
2. **Visualize RNA traces**:
    Generate heatmaps showing RNA intensity over time across all cells.
3. **Hidden Markov Model (HMM) analysis**:
   - Use a two-state (ON/OFF) Gaussian HMM to segment transcriptional states.
   - Estimate model parameters (state transition probabilities, means, variances) using the EM algorithm.
   - Infer the most probable state sequence for each cell over time (using the Viterbi algorithm).
4. **Burst duration fitting**:
   - Collect ON and OFF state durations across all cells.
   - Fit one-, two-, and three-component exponential models to the complementary cumulative distribution (1-CDF) of burst durations.

- **Outputs**:
  - Updated RNA track data including inferred states and burst statistics
  - Plots of heatmaps, single-cell traces, and fitted burst duration distributions
  - Estimated parameters of the HMM and exponential fits

This step produces a quantitative description of transcriptional bursting dynamics, ready for statistical comparison between experimental groups.

## Citation

If you use this code, please cite the original paper.

For questions, please open an Issue on GitHub or contact the author: wangbo@stu.pku.edu.cn.

