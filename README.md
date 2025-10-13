# Pipeline for generation of Patient-specific Voxel Phantom
A Python pipeline for the analysis of segmented medical images (`.nrrd` labelmaps). This pipeline performs statistical analysis, generates 2D slice visualizations, and converts the data into patient-specific voxel phantoms (`.vox` files) compatible with the PENELOPE/penEasy Monte Carlo (MC) simulation framework.

## Introduction
In medical physics and radiation dosimetry, accurate modeling of patient anatomy is crucial for computational radiation dosimetry employing MC simulations. This pipeline addresses the need for a streamlined, automated tool to convert segmented 3D medical images (in NRRD format) into patient-specific computational phantoms.

It processes a labelmap file where each integer value corresponds to a specific organ or tissue, alongside an organ list and an optional CT scan. The primary output is a `.vox` file, ready for use in advanced dosimetry simulations with PENELOPE/penEasy, along with a comprehensive report and slice-by-slice visualizations.

This tool is designed for medical physicists, researchers, and students working with computational dosimetry and medical imaging.

## Description
The key features of the pipeline are:
* **NRRD Processing**: Loads segmented phantom labelmaps and optional CT data from `.nrrd` files.
* **Automatic Cropping**: Intelligently removes empty space (air) around the phantom to optimize phantom size and simulation speed.
* **CT Data Integration**: Converts Hounsfield Units (HU) from a CT scan into mass density (g/cm³) using a user-provided calibration curve.
* **Statistical Analysis**: Calculates key statistics for each organ, including volume, mass, mean/median density, Voxel Bounding Box (VBB), and Center of Mass (COM).
* **Data Integrity Checks**: Includes several safeguards:
    * Corrects negative densities that may arise from CT conversion.
    * Replaces densities of organs with large deviations from ICRP reference values, preventing unrealistic data.
    * Ensures consistency by setting air densities to ICRP reference values.
* **Phantom Generation**: Creates two types of voxel phantoms:
    1.  `phantom_avgden.vox`: A phantom where all voxels of an organ are assigned its average density.
    2.  `phantom_ctden.vox`: A more detailed phantom where each voxel has a specific density derived from the CT scan (if provided).
* **Visualization**: Generates and saves 2D slice images of the phantom along all three anatomical planes (Axial, Sagittal, Coronal) for visual verification.
* **Comprehensive Reporting**: Outputs a detailed `report.out` file containing all calculated statistics in well-formatted tables.


This script orchestrates the entire process, including:
1. Loading a segmented phantom labelmap from an NRRD file.
2. Loading a corresponding organlist from a CSV file.
3. Cropping the phantom to remove surrounding empty space.
4. Optionally loading a CT NRRD file, converting its HU values to density, and calculating organ-specific statistics.
5. Performing data integrity checks and creating phantom arrays based on average organ densities or per-voxel CT densities.
6. Generating and saving 2D slice visualizations of the phantoms.
7. Calculating Voxel Bounding Box (VBB) and Center of Mass (COM) for all organs.
8. Optionally creating .vox files formatted for Monte Carlo simulations.
9. Generating a comprehensive 'report.out' file with all calculated data tables.

## Dependencies/Prerequisites

To run this script, you'll need Python 3.x, as well as the following libraries, which are included with the anaconda distribution, except for pynrrd.
- Numpy
- Pandas
- Matplotlib
- Tabulate
- Time
- [pynrrd](https://github.com/mhe/pynrrd)

You can also install the dependencies using pip, via the terminal command:
Example:

```
pip install numpy
```

And use analogous commands for the other libraries.

It's possible to install pynrrd using pip or anaconda, via the following terminal commands:

```
pip install pynrrd
```

or 

```
conda install conda-forge::pynrrd
```


## Execution

To get started, place the script and your input files in the same directory. Then, simply run the script from your terminal:

```
python nrrd_to_3d_pipeline.py
```

The script will then guide you through a highly user-friendly and configurable through the command-line interface, consisting of a series of interactive prompts. Most parameters, like file names and phantom dimensions, are provided during runtime.

There are two output streams: terminal and file 'report.out'. There is info that is outputted to both streams. However, while the file 'report.out' is more for reporting results/calculations performed in the phantom, the terminal is more for reporting the progress of execution, if there are issues with data integrity, etc.

During execution of the program, the execution time is prompted various times throughout the program, so the user can understand which steps of processing take more or less time.

When the user is prompted for a filename (or path if the file is not in the same durectory as the script), a default filename is provided. The default filename will only be used if the user leaves the input prompt blank.

## Example


## Disclaimer

This software is provided "as is" for academic, scientific, and research purposes only. It is not a medical device, nor intended for clinical use, and has not been approved by any medical regulatory authority (e.g., EMA, FDA). The authors and contributors provide this software "as is", without any warranty of any kind, either expressed or implied. This includes, but is not limited to, the implied warranties of merchantability and fitness for a particular purpose. 

The users are solely responsible for data privacy management when using this program.
- Users must ensure they have the legal and ethical authority to use and process any patient data.
- Users must de-identify and anonymize all personal health information (PHI) in accordance with all applicable data protection laws and regulations (e.g., GDPR in Europe, HIPAA in the United States) before using it with this software.
- The authors of this software take no responsibility for data breaches or misuse of personal data by users of this script.

By using this software, you acknowledge that you have read, understood, and agree to be bound by the terms of this disclaimer.

## License
Copyright (C) 2025 Jorge Cebola Borbinha
This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License Version 3 (AGPLv3) as published by the Free Software Foundation.
See the GNU Affero General Public License for more details, available at: [AGPL v3 LICENSE](https://github.com/jorge-borbinha/Patient-specific_voxel_phantom_generation/blob/main/LICENSE.md)


