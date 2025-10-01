# Pipeline for generation of Patient-specific Voxel Phantom
Pipeline for patient-specific voxel phantom file generation


## Contextualization


## Description


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


## Usage

To get started, place the script and your input files in the same directory. Then, simply run the script from your terminal:

```
python nrrd_to_3d_pipeline.py
```

The script will then guide you through a command-line interface, consisting of a series of interactive prompts, to gather all the necessary phantom information and file names.


## Execution

The script is highly user-friendly and configurable through the command-line interface. Most parameters, like file names and phantom dimensions, are provided during runtime

## Example


