# Cosmic-Ray-Rejection

# Cosmic Ray Rejection in Astronomical Imaging

## Overview
This project implements methods for identifying and removing cosmic ray artifacts in astronomical imaging and spectroscopic data.
Cosmic rays introduces high-intensity pixels that can significantly impact scientific analysis, making their detection and removal a critical step in data reduction pipelines.

The project focuses on developing robust algorithms to detect and mask these artifacts in both multi-frame and single-frame datasets.

## Objectives
- Detect and flag cosmic ray–affected pixels in imaging and spectroscopic data
- Construct bias-corrected and cleaned science frames
- Develop reusable routines for bad pixel masking
- Compare multi-frame and single-frame cosmic ray rejection strategies

## Methods

### 1. Bias Correction and Image Combination (Keck/LRIS Blue)
- Constructed a median bias frame from calibration images
- Subtracted bias from science frames
- Combined multiple exposures to enhance signal-to-noise ratio
- Identified cosmic ray artifacts via outlier detection

### 2. Cosmic Ray Rejection in Multi-Frame Data (LRIS Red)
- Applied statistical comparison across multiple frames
- Flagged pixels with anomalously high counts relative to neighboring frames
- Generated bad pixel masks
- Produced cleaned, combined spectral images

### 3. Cosmic Ray Detection in Single Frames
- Developed methods to detect cosmic rays without multiple exposures
- Used local pixel statistics and spatial filtering
- Identified sharp, high-intensity outliers characteristic of cosmic ray hits
- Generated bad pixel masks for individual images

## Tools & Technologies
- Python
- NumPy
- Matplotlib
- Astropy (FITS file handling)

## Results
- Successful identification of cosmic ray–affected pixels
- Improved image quality after artifact removal
- Robust bad pixel masks for both multi-frame and single-frame datasets
- Demonstrated effectiveness of statistical and spatial filtering techniques

## Data
The dataset consists of Keck/LRIS spectroscopic observations. Due to file size constraints, raw FITS files are not included in this repository.
