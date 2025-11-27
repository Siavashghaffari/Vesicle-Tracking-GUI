# Vesicle Transcytosis Analysis Tool

A MATLAB-based graphical user interface (GUI) application for automated detection, tracking, and quantification of vesicle transcytosis events in endothelial cells using Total Internal Reflection Fluorescence (TIRF) microscopy.

## Overview

This tool provides comprehensive analysis of endothelial transcytosis by processing time-lapse TIRF microscopy image sequences. The software automatically identifies vesicle trafficking events, tracks individual vesicles over time, and quantifies transcytosis based on fluorescence intensity dynamics. The methodology implemented in this tool has been published and validated for studying cargo transport across endothelial barriers.

## Citation

If you use this software in your research, please cite the following protocol paper where the methodology is described:

**Ghaffari, S., Jang, E. & Lee, W. L. (2022).** Quantifying endothelial transcytosis with total internal reflection fluorescence microscopy (TIRF). In *Fluorescent Microscopy* (pp. 115-124). New York, NY: Springer US.

### Related Publications

This tool has been utilized in the following research:

1. **Ghaffari, S., Naderi Nabi, F., Sugiyama, M. G., & Lee, W. L. (2018).** Estrogen inhibits LDL (low-density lipoprotein) transcytosis by human coronary artery endothelial cells via GPER (G-protein–coupled estrogen receptor) and SR-BI (scavenger receptor class B type 1). *Arteriosclerosis, Thrombosis, and Vascular Biology*, 38(10), 2283-2294.

2. **Ghaffari, S., Jang, E., Naderinabi, F., Sanwal, R., Khosraviani, N., Wang, C., ... & Lee, W. L. (2021).** Endothelial HMGB1 is a critical regulator of LDL transcytosis via an SREBP2–SR-BI axis. *Arteriosclerosis, Thrombosis, and Vascular Biology*, 41(1), 200-216.

## Features

- **Automated Vesicle Detection**: Adaptive thresholding algorithm for identifying vesicle particles
- **Particle Tracking**: Multi-frame trajectory reconstruction with gap-closing capabilities
- **Transcytosis Event Identification**: Intensity-based detection of exocytosis events
- **Quantitative Analysis**:
  - Mean squared displacement (MSD) calculations
  - Diffusion coefficient determination
  - Trajectory characterization (length, velocity, displacement)
  - Transcytosis fraction quantification
- **Quality Filtering**: Size and circularity-based filtering of detected events
- **Batch Processing**: Analyze multiple image sequences in a single run
- **Data Export**: Comprehensive output in MATLAB format and tab-delimited text files
- **Visualization**: Automated generation of analysis overlays and diagnostic images

## System Requirements

### Software Dependencies
- MATLAB (R2015a or later recommended)
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox (optional, for enhanced analysis)

### Hardware Requirements
- Minimum 8 GB RAM (16 GB recommended for large datasets)
- Multi-core processor recommended for improved performance

## Installation

1. Clone or download this repository to your local machine:
   ```bash
   git clone https://github.com/yourusername/Vesicle-Tracking-GUI.git
   ```

2. Add the tool directory to your MATLAB path:
   ```matlab
   addpath(genpath('/path/to/Vesicle-Tracking-GUI'))
   ```

3. Verify installation by launching the GUI:
   ```matlab
   mahsa  % Launches the mahsa.m GUI
   ```

## Usage

### Quick Start

1. **Launch the GUI**:
   ```matlab
   mahsa  % Launches the mahsa.m GUI
   ```

2. **Select Input Directory**: Click the folder selection button and navigate to the directory containing your TIFF image sequences.

3. **Configure Analysis Parameters** (optional): Adjust detection and tracking parameters via the GUI interface.

4. **Select Output Directory**: Choose where analysis results will be saved.

5. **Run Analysis**: Click the "Run" button to begin processing.

6. **Review Results**: Analysis outputs will be saved in automatically generated folders with comprehensive tracking data and visualizations.

### Input Data Format

- **File Type**: TIFF image stacks (.tiff, .tif)
- **Image Requirements**:
  - Time-lapse sequences from TIRF microscopy
  - Single-channel fluorescence images
  - 8-bit or 16-bit grayscale
  - Sequential numbering recommended

### Analysis Pipeline

The tool executes the following workflow:

1. **Image Loading**: Import TIFF stacks and normalize intensity ranges
2. **Preprocessing**: Gaussian filtering and local background correction
3. **Vesicle Detection**: Adaptive thresholding based on local intensity
4. **Size and Shape Filtering**: Remove artifacts based on area and circularity
5. **Particle Tracking**: Link detected particles across frames using nearest-neighbor algorithm
6. **Intensity Analysis**: Measure fluorescence intensity dynamics for each track
7. **Transcytosis Classification**: Identify exocytosis events based on intensity slope changes
8. **Quantification**: Calculate transcytosis fraction and track statistics
9. **Data Export**: Generate analysis files and diagnostic visualizations

### Key Analysis Parameters

The following parameters can be adjusted for optimal detection (defaults shown):

- **minSize**: Minimum vesicle diameter (pixels) - *default: 3*
- **maxSize**: Maximum vesicle diameter (pixels) - *default: 9*
- **circFilter**: Circularity threshold (0-1, where 1 = perfect circle) - *default: 0.2*
- **cutOff**: Intensity threshold multiplier over local background - *default: 1.11*
- **psfFWHM**: Point spread function full-width half-maximum (pixels) - *default: 0.5*
- **Shift**: Maximum inter-frame displacement (pixels) - *default: 8*
- **mem**: Maximum number of frames a particle can disappear before track termination - *default: 2*
- **good**: Minimum track length (frames) to be included in analysis - *default: 5*
- **slopeCutoff**: Standard deviations required for exocytosis detection - *default: 2.5*

To modify these parameters, edit the default values in `transcytosisIntensity.m` (lines 54-77) or pass custom `InputCond` and `TrackCond` structures to the function.

## Output Files

For each analyzed image sequence, the tool generates:

### Data Files
- **transcytosis_analysis.mat**: Complete MATLAB structure containing all analysis results
- **ExoData.mat**: Aggregated data from multiple experiments
- **FinalData.txt**: Tab-delimited summary table with:
  - Filename
  - Number of transcytosis events
  - Total number of tracks
  - Transcytosis fraction

### Image Outputs
- **Pre-filter threshold.tif**: Binary masks before filtering
- **Filtered Threshold.tif**: Binary masks after size/shape filtering
- Additional diagnostic images (optional)

### Data Structure

The `outputData` structure contains:

```matlab
outputData
├── inputInfo           % Analysis parameters and file information
├── trackCond           % Tracking parameters
├── sites               % Detected vesicle properties per frame
├── sTable              % Particle position table
├── tracks              % Linked particle trajectories
├── trackData           % Individual track analysis
│   ├── trackLength     % Duration of each track (seconds)
│   ├── pos             % X-Y coordinates
│   ├── vel             % Instantaneous velocities
│   ├── DcoEff          % Diffusion coefficient
│   ├── Gamma           % MSD curve exponent
│   └── RawPos          % Intensity measurements
├── populationStats     % Aggregate statistics
├── possibleExo         % Candidate transcytosis events
├── confirmedExo        % Validated transcytosis events
└── graphData           % Visualization data
```

## Methodology

### Vesicle Detection

The detection algorithm employs adaptive thresholding to identify fluorescent vesicles:

1. Gaussian filtering based on estimated point spread function
2. Local background calculation using moving average filter
3. Ratio-based thresholding (signal/background > cutoff)
4. Morphological filtering based on size and circularity

### Particle Tracking

Trajectories are constructed using a nearest-neighbor algorithm with the following capabilities:

- Frame-to-frame particle linking based on maximum displacement criterion
- Gap-closing to handle temporary particle disappearance
- Minimum track length filtering to eliminate spurious detections
- Unique ID assignment for each trajectory

The tracking implementation is based on the algorithm described in:
> Crocker, J. C., & Grier, D. G. (1996). Methods of digital video microscopy for colloidal studies. *Journal of Colloid and Interface Science*, 179(1), 298-310.

### Transcytosis Event Detection

Transcytosis events are identified by analyzing fluorescence intensity dynamics:

1. **Intensity Measurement**: Vesicle intensity quantified at each time point
2. **Slope Calculation**: Frame-to-frame intensity changes computed
3. **Baseline Characterization**: Mean and standard deviation of intensity slope
4. **Event Detection**: Significant negative slope change indicating cargo release
5. **Validation**: Confirmation based on track endpoint and mobility criteria

An event is classified as transcytosis when:
- Final intensity slope < (mean slope - 2.5 × SD)
- Track ends before final frame (indicating vesicle disappearance)
- Gamma (MSD exponent) < 0.87 (excluding highly immobile particles)

### Diffusion Analysis

For each trajectory, the tool calculates:

- **Mean Squared Displacement (MSD)**: `MSD(τ) = ⟨[r(t+τ) - r(t)]²⟩`
- **Diffusion Coefficient**: From the linear portion of MSD vs. time
- **Anomalous Diffusion Exponent (Gamma)**: From fit to `MSD = 4Dτ^γ`
  - γ = 1: Normal diffusion
  - γ < 1: Subdiffusion (confined/hindered motion)
  - γ > 1: Superdiffusion (directed motion)

## Troubleshooting

### Common Issues

**Issue**: No particles detected
- **Solution**: Reduce `cutOff` parameter or verify image quality and contrast

**Issue**: Too many false positives
- **Solution**: Increase `minSize`, `circFilter`, or `cutOff` parameters

**Issue**: Tracks frequently broken
- **Solution**: Increase `Shift` (max displacement) or `mem` (memory frames)

**Issue**: No transcytosis events detected
- **Solution**: Verify that vesicles undergo exocytosis in your dataset, or adjust `slopeCutoff`

**Issue**: Out of memory errors
- **Solution**: Process smaller image stacks or increase available RAM

## Advanced Usage

### Programmatic Access

For batch processing or custom workflows, call the core function directly:

```matlab
% Define analysis parameters
InputCond.minSize = 3;
InputCond.maxSize = 9;
InputCond.circFilter = 0.2;
InputCond.cutOff = 1.11;
InputCond.psfFWHM = 0.5;
InputCond.Shift = 8;
InputCond.timeVect = [0 0.15];  % [start_time interval_in_seconds]
InputCond.scale = 0.151498;      % µm per pixel
InputCond.saveGraph = 1;
InputCond.saveVid = 1;
InputCond.saveDiag = 0;
InputCond.finalSteps = 2;
InputCond.slopeCutoff = 2.5;

% Define tracking parameters
TrackCond.mem = 2;
TrackCond.dim = 2;
TrackCond.good = 5;
TrackCond.quiet = 1;

% Run analysis
outputData = transcytosisIntensity('input_file.tif', InputCond, TrackCond);
```

### Custom Analysis Workflows

The modular design allows integration into custom pipelines:

```matlab
% Load image stack
imageStack = LoadStack('mydata.tif');

% Process with custom parameters
results = transcytosisIntensity('mydata.tif', InputCond, TrackCond);

% Aggregate results from multiple files
aggregatedData = colDataLee();

% Access specific metrics
transcytosisFraction = aggregatedData.populationData(:,3);
```

## Contributing

Contributions to improve the tool are welcome. Please consider:

- Bug reports and feature requests via GitHub Issues
- Code contributions via Pull Requests
- Documentation improvements
- Example datasets and analysis protocols



## Authors

This work was developed by **Siavash Ghaffari**. For any questions, feedback, or additional information, please feel free to reach out. Your input is highly valued and will help improve and refine this pipeline further.



## Version History

- **v1.0** (May 2017): Initial release with GUI interface and core analysis pipeline
- Current version includes batch processing, enhanced filtering, and comprehensive data export

---

**Note**: This tool is designed specifically for TIRF microscopy data. Adaptation to other imaging modalities may require parameter optimization and validation.
