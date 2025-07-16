# Analysis scripts for condensate analysis and burst modeling in biological imaging data
## Overview

This repository contains scrips associated with the paper "Pulsatile fusion of enhancer and promoter condensates drives transcriptional bursting". It aims to support high-throughput quantitative analysis of super-resolution images and time series data, focusing on dynamic processes such as transcriptional bursting and condensate dynamics.

## Parts of burst modeling

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

## Parts of condensate analysis

This MATLAB script processes 4-channel time-lapse fluorescence images at the single-cell level, extracting foci positions, intensities, and condensate boundaries for downstream analysis.

- **Inputs**:
  - Multi-channel `.tif` image stacks (DNA, RNA, OCT4, BRD4)
  - TrackMate `.txt` files containing single-cell foci coordinates
- **Main steps**:

1. **Load data**:
   - Read 4-channel image series and organize by frame, z-slice, and channel.
   - Extract max-intensity-foci Z-layer for each channel.
2. **Nucleus segmentation**:
   - Use Gaussian filtering and intensity thresholding to create a nucleus mask from the OCT4 channel.
3. **Single-cell foci detection (DNA and RNA channels)**:
   - Refine foci positions frame by frame.
   - Compute background and normalized intensity around each foci.
   - Crop ROIs centered on the foci for further analysis.
4. **Condensate analysis (OCT4 and BRD4 channels)**:
   - Apply rolling ball background subtraction and Gaussian filtering.
   - Perform HMRF segmentation to identify condensate regions and determine intensity thresholds.
   - Extract condensate center, interface, and boundary.
5. **Spatial measurements**:
   - Compute distances of DNA and RNA foci to the condensate boundary across all frames.
   - Calculate the distance between DNA and RNA foci within the same nucleus.

- **Outputs**:
  - Processed ROI images and segmentation masks in `.tif` format.
  - MATLAB `.mat` file containing:
    - Foci positions, intensities, and background
    - Condensate segmentation results and thresholds
    - Distance measurements between foci and condensates

This script transforms raw time-lapse image data into quantitative single-cell and subnuclear features for later statistical analysis.