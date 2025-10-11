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
See the GNU Affero General Public License for more details, available at: [AGPL v3](https://www.gnu.org/licenses/agpl-3.0.html)


