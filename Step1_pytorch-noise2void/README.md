# pytorch-noise2void

A lightweight PyTorch implementation of Noise2Void-style self-supervised image denoising. This folder contains training and testing pipelines, example datasets, pretrained weights, and export helpers.

### Contents

- `datasets/` — example / expected dataset layout
- `noise2void/` — core Noise2Void model and training code
- `onnx/` — exported ONNX models
- `pth/` — saved PyTorch checkpoints
- `requirements.txt` — Python dependencies
- `train_pipeline.py` — training script
- `test_pipeline.py` — inference / evaluation script

### Quick start

Create an environment and install dependencies use `conda`:

```bash
conda create -n n2v python=3.10 -y
conda activate n2v
pip install -r requirements.txt
```

### Training

Prepare your dataset under `datasets/` with a simple folder structure, e.g.:

```
datasets/
  train/
    img_000.tif
    img_001.tif
  val/
    img_000.tif
```

Specify the dataset path and then run the training pipeline:

```bash
python train_pipeline.py
```

Checkpoints will be written to `pth/`  and `onnx/`.

### Inference / Evaluation

Run denoising / evaluation with:

```bash
python test_pipeline.py
```

### Dataset expectations

- Single-channel (grayscale) images supported, common microscopy formats (TIFF, ND2).
- Training uses patches sampled from noisy images (self-supervised). Check `noise2void/` for exact preprocessing and augmentation details.

### Citation

If you use this code, please cite the original paper.

For questions, please open an Issue on GitHub or contact the author: wangbo@stu.pku.edu.cn.
