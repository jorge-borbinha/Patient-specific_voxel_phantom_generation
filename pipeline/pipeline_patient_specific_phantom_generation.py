#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Python pipeline for analysis of NRRD segmented phantom labelmap files: 
# performs statistical analysis, visualization and converts to .vox files 
# readable by Monte Carlo software.
#
# @author:  Jorge Cebola Borbinha
# @github:  jorge-borbinha
# @website: jorge-borbinha.github.io
#
# Copyright (C) 2025 Jorge Cebola Borbinha
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Mozilla Public License 2.0
#
# For more information, see https://www.mozilla.org/en-US/MPL/2.0/
#
# This software is for educational and research purposes only. The user is 
# solely responsible for ensuring data privacy and must de-identify all patient
# data before use. The authors provide this software "as is" without any
# warranty. By using this software, you agree to these terms.
#

import numpy as np
import pandas as pd
from tabulate import tabulate
import matplotlib.pylab as plt
from scipy import ndimage
import nrrd
import time
import os
import sys

#### #### #### #### #### #### #### #### #### #### #### #### ####
#### CONFIGURATION #### #### #### #### #### #### #### #### ####

DENS_MAX_REL_DIFF = 75 # Max relative difference (75%) before using ICRP density.
DEFAULT_NRRD_FILE = 'phantom_labelmap.nrrd'
DEFAULT_ORGANLIST_FILE = 'organlist.csv'
DEFAULT_CT_FILE = 'ct.nrrd'
DEFAULT_VOX_AVG_FILENAME = 'phantom_avgden.vox'
DEFAULT_VOX_CT_FILENAME = 'phantom_ctden.vox'

#### #### #### #### #### #### #### #### #### #### #### #### ####
#### DEFINE FUNCTIONS #### #### #### #### #### #### #### #### #### 

def report_end_time(start_time):
    """Calculates and prints the total execution time of the script."""
    # get the end time
    end_time = time.time()
    # get the execution time in minutes
    elapsed_time = (end_time - start_time)/60
    print(f"\n> Execution time: {elapsed_time:2.1f} minutes.\n")
    pass

def prompt_yes_no(question):
    """Prompts the user with a yes/no question and returns their choice."""
    while True:
        choice = input(f">> {question} (y/n): ").lower()
        if choice in ['y', 'n']:
            return choice
        print("Invalid input. Please enter 'y' or 'n'.")

def format_header(title_text, width=9, border_char=">>>> "):
    """
    Creates a formatted, decorative header with a prompted title, to help
    making the output more readable.

    Args:
        title_text (str): Text to be displayed in the header.
        width (int): Total width of the header. Times to repeat the border_char on top and bottom lines.
        border_char (str): Character to use width times for the top and bottom borders.

    Returns:
        str: A formatted, multi-line string representing the header.
    """
    # Create the top and bottom border lines by repeating the character
    border_line = border_char * width

    # Combine the border and title into the final header string
    return f"{border_line}\n{border_char}{border_char} {title_text}\n{border_line}\n"

def get_file_path_from_user(prompt, default_path=None):
    """
    Prompts the user for a file path. Suggests a default if provided. If the 
    program does not provide a default value, the prompt is returned to the user.

    Args:
        prompt (str): Question to ask the user.
        default_path (str, optional): Default file path.

    Returns:
        str: File path provided by the user.
    """
    if default_path:
        prompt_with_default = f">> {prompt} [default: {default_path}]: "
        path = input(prompt_with_default)
        if not path:
            return default_path
        return path
    else:
        return input(f"{prompt}: ")

def _correct_ct_densities(phantom_ct, phantom, organ_df):
    """
    Corrects the CT density map by setting air densities to ICRP reference
    and ensuring no negative densities exist.

    Args:
        phantom_ct (np.ndarray): The 3D array of densities calculated from CT HU values.
        phantom (np.ndarray): The 3D array of organ IDs.
        organ_df (pd.DataFrame): DataFrame with organ information.

    Returns:
        np.ndarray: The corrected density map.
    """
    print("\nCorrecting densities: Setting air to ICRP values and handling negative densities...\n")
    # Set Air Density in phantom_ct to ICRP reference value for consistency
    # Organ ID 0 is air outside the body/phantom, the max Organ ID is air inside body
    air_outside_phantom_mask = (phantom == organ_df["Organ_ID"].iloc[0])
    phantom_ct[air_outside_phantom_mask] = organ_df["Density_ICRP"].iloc[0]
    air_inside_phantom_mask = (phantom == organ_df["Organ_ID"].iloc[-1])
    phantom_ct[air_inside_phantom_mask] = organ_df["Density_ICRP"].iloc[-1]
    print("For consistency, the ICRP reference air density was set for air inside and ouside phantom.")
    
    # Do not allow negative densities
    if np.any(phantom_ct < 0):
        mask = phantom_ct < 0
        print(f"There are negative densities in CT file, comprising a total of {mask.sum()} voxels.")
        print(f"The average of negative desities is {phantom_ct[mask].mean():.4f} g/cm^3.")
        print(f"The standard deviation of negative desities is {phantom_ct[mask].std():.4f} g/cm^3.")
        print("For the voxels containing negative densities, the ICRP reference air density will be passed on.\n")
        phantom_ct[mask] = organ_df["Density_ICRP"].iloc[0]
        
    return phantom_ct

def calculate_ct_statistics(phantom_labelmap, ct_density_map, organ_df, voxel_volume_mm3):
    """
    Calculates statistics for each organ/structure directly from the CT data.

    Args:
        phantom_labelmap (np.ndarray): 3D array of organ IDs.
        ct_density_map (np.ndarray): 3D array of densities g/cm^3.
        organ_df (pd.DataFrame): dataframe with organ information.
        voxel_volume_mm3 (float): voxel volume mm^3.

    Returns:
        ct_df(pd.DataFrame): dataframe with detailed statistics for each organ.
    """
    print("Calculating statistics from CT data for each organ and storing to database...\n")
    stats_list = []
    for organ_id in organ_df['Organ_ID']:
        # Create a mask to isolate voxels for the current organ
        mask = (phantom_labelmap == organ_id)
        
        # If the organ is not present in the phantom, fill with zeros
        if not np.any(mask):
            stats = {
                "Organ_ID": organ_id, "Organ": organ_df.loc[organ_id, "Organ"],
                "Den_Mean": 0, "Den_Median": 0, "Den_Standard_deviation": 0
            }
            stats_list.append(stats)
            continue

        # Get the HU and Density values for the current organ
        density_values = ct_density_map[mask]

        # Calculate statistics
        stats = {
            "Organ_ID": organ_id,
            "Organ": organ_df.loc[organ_id, "Organ"],
            "Density_Mean": np.mean(density_values),
            "Density_Median": np.median(density_values),
            "Density_Standard_deviation": np.std(density_values)
        }
        stats_list.append(stats)
    
    # Create the final DataFrame
    ct_df = pd.DataFrame(stats_list)
    return ct_df

def auto_crop_phantom(ph_data, padding=2):
    """
    Automatically calculates the bounding box for non-air structures and crops
    the phantom in x, y and z dimensions with a specified padding (i.e. voxel 
    layer tolerance beyond the phantom body). Voxels with Organ_ID = 0 represent 
    air outside phantom (and possibly other structures) and are cropped away.

    Args:
        ph_data (np.ndarray): 3D phantom data array retreived from nrrd file
        padding (int): voxel layer tolerance to add around the phantom bounding 
                       box on each side (x, y, z)

    Returns:
        tuple: A tuple containing:
            - np.ndarray: cropped phantom
            - tuple: The crop indices (crop_x_min, crop_x_max, crop_y_min, 
                     crop_y_max, crop_z_min, crop_z_max), which can be reused 
                     to crop other arrays (e.g. CT data)
    """
    # Air outside body (and other structures) is assumed to have Organ_ID = 0
    # Find the coordinates of all non-air voxels, i.e. the voxels that were segmented in the nrrd labelmap file
    non_air_out_mask = (ph_data != 0)
    bin_phantom = non_air_out_mask.astype(int)

    # Edge case where the nrrd file is not segmented, i.e. all Organ_ID = 0
    if non_air_out_mask.size == 0:
        print("Warning: No segmented structures found in the phantom. No cropping will be performed.\n")
        # Return the original phantom and full-dimension slices
        crop_x_min, crop_x_max = 0, ph_data.shape[0]
        crop_y_min, crop_y_max = 0, ph_data.shape[1]
        return ph_data, (crop_x_min, crop_x_max, crop_y_min, crop_y_max)

    # Determine the min and max indices for x and y axes from non-air phantom
    phantom_location = ndimage.find_objects(bin_phantom)
    x_axis, y_axis, z_axis = phantom_location[0]
    
    x_min, x_max = x_axis.start, x_axis.stop
    y_min, y_max = y_axis.start, y_axis.stop 
    z_min, z_max = z_axis.start, z_axis.stop 

    # Apply padding and ensure indices remain within the array's bounds
    # use max(0, ...) to prevent negative indices and min(shape, ...) to prevent indices from exceeding the array size
    # For the min boundary, subtract the padding
    crop_x_min = max(x_min - padding, 0)
    crop_y_min = max(y_min - padding, 0)
    crop_z_min = max(z_min - padding, 0)
    
    # For the max boundary, add the padding, +1 to max values -> NumPy's Inclusive-start vs. Exclusive-stop slicing rule
    crop_x_max = min(x_max + padding + 1, bin_phantom.shape[0])
    crop_y_max = min(y_max + padding + 1, bin_phantom.shape[1])
    crop_z_max = min(z_max + padding + 1, bin_phantom.shape[2])
    
    # Effectively crop the phantom data
    cropped_phantom = ph_data[crop_x_min:crop_x_max, crop_y_min:crop_y_max, crop_z_min:crop_z_max]

    # The slice info is returned to be reused for the CT data
    crop_indices = (crop_x_min, crop_x_max, crop_y_min, crop_y_max, crop_z_min, crop_z_max)
    
    return cropped_phantom, crop_indices

def visualize_phantom_slice(phantom_data, voxel_size, slice_axis, slice_index, cmap='viridis', save_path=None):
    """
    Visualizes 2D slices of a 3D phantom data array with correct aspect ratio 
    in the along the x, y and z axis (these correspond to the sagittal, coronal
    and axial planes, respectively). 3D visualization of the phantom (even using
    plotly) is quite computationally demanding, due to the very high number of 
    voxels in a voxel phantom.

    Args:
        phantom_data (np.ndarray): The 3D array to be visualized (e.g., phantom_ct).
        voxel_size (np.ndarray): The physical dimensions of the voxels (from the header).
        slice_axis (str): The axis to slice along ('x', 'y', or 'z').
        slice_index (int): The index of the slice to be displayed.
        cmap (str, optional): The colormap to use for the plot. Defaults to 'viridis'.
    """
    # Check if the chosen index is valid for the specified axis
    if slice_axis == 'x' and slice_index >= phantom_data.shape[0]:
        print(f"Error: Column index {slice_index} is out of bounds for the x-axis (max index is {phantom_data.shape[0] - 1}).")
        return
    elif slice_axis == 'y' and slice_index >= phantom_data.shape[1]:
        print(f"Error: Row index {slice_index} is out of bounds for the y-axis (max index is {phantom_data.shape[1] - 1}).")
        return
    elif slice_axis == 'z' and slice_index >= phantom_data.shape[2]:
        print(f"Error: Slice index {slice_index} is out of bounds for the z-axis (max index is {phantom_data.shape[2] - 1}).")
        return

    # Select the appropriate 2D slice and calculate the aspect ratio
    if slice_axis == 'x':
        slice_2d = phantom_data[slice_index, :, :]
        plt_xlabel = 'Y-axis (no. voxels)'
        plt_ylabel = 'Z-axis (no. voxels)'
        title = f'YZ-Plane (Sagittal Plane) Column at X-index {slice_index}'
        aspect_ratio = voxel_size[2] / voxel_size[1]
    elif slice_axis == 'y':
        slice_2d = phantom_data[:, slice_index, :]
        plt_xlabel = 'X-axis (no. voxels)'
        plt_ylabel = 'Z-axis (no. voxels)'
        title = f'XZ-Plane (Coronal Plane) Row at Y-index {slice_index}'
        aspect_ratio = voxel_size[2] / voxel_size[0]
    elif slice_axis == 'z':
        slice_2d = phantom_data[:, :, slice_index]
        plt_xlabel = 'X-axis (no. voxels)'
        plt_ylabel = 'Y-axis (no. voxels)'
        title = f'XY-Plane (Axial Plane) at Z-index {slice_index}'
        aspect_ratio = voxel_size[1] / voxel_size[0]
    else:
        print("Error: Invalid slice_axis. Please choose 'x', 'y', or 'z'.")
        return

    # Dynamically calculate figure size to match the physical aspect ratio of phantom slice
    # to prevent excessive white space around the plot -> because phantom size can vary
    # significantly
    # Determine physical height and width of the slice to be plotted (in mm)
    ## NumPy's Data Order: When we slice the 3D phantom for an XY-plane, you get a 2D array 
    ## with the shape (no_x_voxels, no_y_voxels). The X-dimension is the first axis 
    ## (the "columns") and the Y-dimension is the second axis (the "rows"). 
    ## This is in voxel phantom notation, numpy notation exchanges columns and rows.
    ## imshow's Plotting Order: By default, imshow plots the array's first axis along the 
    ## vertical Y-axis and the second axis along the horizontal X-axis.
    ## If we plotted the slice directly, imshow would place the data's X-axis vertically and 
    ## its Y-axis horizontally, resulting in an image that is rotated 90 degrees from the 
    ## expected anatomical view.
    
    if slice_axis == 'x':  # YZ-Plane: plot is Z vs Y, slice_2d.T has shape (Nz, Ny)
        physical_height = slice_2d.T.shape[0] * voxel_size[2]
        physical_width = slice_2d.T.shape[1] * voxel_size[1]
    elif slice_axis == 'y':  # XZ-Plane: plot is Z vs X, slice_2d.T has shape (Nz, Nx)
        physical_height = slice_2d.T.shape[0] * voxel_size[2]
        physical_width = slice_2d.T.shape[1] * voxel_size[0]
    else:  # XY-Plane: plot is Y vs X, slice_2d.T has shape (Ny, Nx)
        physical_height = slice_2d.T.shape[0] * voxel_size[1]
        physical_width = slice_2d.T.shape[1] * voxel_size[0]

    # Define base size for figure longest dimension in inches for readability
    base_size = 8.0

    # Calculate the figure width and height to maintain the physical aspect ratio
    if physical_width >= physical_height:
        fig_width = base_size
        fig_height = base_size * (physical_height / physical_width)
    else:
        fig_height = base_size
        fig_width = base_size * (physical_width / physical_height)

    # Add some extra width to the figure to include colorbar
    fig_width += 2.0
    
    # Calculate the global min and max of the entire phantom to ensure a consistent colormap across all plots
    vmin = np.nanmin(phantom_data)
    vmax = np.nanmax(phantom_data)

    # Create the plot with the new, dynamically calculated figure size
    plt.figure(figsize=(fig_width, fig_height))
    
    # Use imshow to display the 2D slice
    # Pass the calculated aspect ratio to the 'aspect' argument
    im = plt.imshow(slice_2d.T, cmap=cmap, origin='lower', aspect=aspect_ratio, vmin=vmin, vmax=vmax)
    
    # Add a color bar to show the density scale
    cbar = plt.colorbar(im)
    cbar.set_label('Density (g/cm³)')
    
    # Set plot labels and title
    plt.xlabel(plt_xlabel)
    plt.ylabel(plt_ylabel)
    plt.title(title)
    
    # Save figure if path is provided and then close it to free memory
    if save_path:
        try:
            # Use tight_layout to ensure all plot elements fit well before saving
            plt.tight_layout()
            plt.savefig(save_path, dpi=150) # Save with a decent resolution
        except Exception as e:
            print(f"Error saving file to {save_path}: {e}")

    # Close plot figure -> essential for batch saving to prevent memory leaks
    plt.close()

def save_all_phantom_slices(phantoms, voxel_size):
    """
    Manages the batch saving of all phantom slices based on user input.

    Args:
        phantoms (dict): A dictionary where keys are phantom names (str) and
                         values are the 3D phantom data (np.ndarray).
        voxel_size (np.ndarray): The physical dimensions of the voxels.
    """
    choice = prompt_yes_no("Do you want to save the plots from XY-,XZ- and YZ-planes to a folder?")

    if choice == 'y':
        save_location = input(">> Enter the full path to the directory where you want to save the plots: ")
        
        # Check if the provided path is a valid directory. Exit if not
        if not os.path.isdir(save_location):
            print(f"\nError: The path '{save_location}' is not a valid directory. The plots will not be saved.")
            # We return here instead of sys.exit to allow the script to continue
            return

        # Create the target folder called "phantom_plots" inside the user-provided path
        base_plots_folder = os.path.join(save_location, "phantom_plots")
        os.makedirs(base_plots_folder, exist_ok=True) # exist_ok=True prevents an error if the folder already exists
        print(f"\nPlots will be saved in: {base_plots_folder} ")
        print("(this process may take some time)\n")
        
        # Define the subdirectories to be created
        phantom_names = ['phantom_avgden', 'phantom_ct']
        plane_names = ['XY', 'XZ', 'YZ']
        
        # Construct the paths and create all the necessary subdirectories in advance
        for phantom_name in phantom_names:
            for plane_name in plane_names:
                # Construct the full path for each subdirectory: e.g., /path/phantom_plots/phantom_avgden/XY
                specific_folder_path = os.path.join(base_plots_folder, phantom_name, plane_name)
                os.makedirs(specific_folder_path, exist_ok=True)

        # Now the directory structure is in place, start saving files
        # Loop through each phantom: phantom_avgden, phantom_ct
        for name, data in phantoms.items():
            # If a phantom is None (e.g., phantom_ct when no file was given), skip it.
            if data is None:
                print(f"\n--- Skipping phantom: {name} (not generated) ---")
                continue

            print(f"\n--- Processing phantom: {name} ---")
            
            # Save all XY-plane slices (iterating along z-axis)
            total_slices = data.shape[2]
            print(f"Saving {total_slices} XY-plane slices...")
            for i in range(total_slices):
                filename = f"{name}_XY_z={i}.jpg"
                full_path = os.path.join(base_plots_folder, name, 'XY', filename)
                visualize_phantom_slice(data, voxel_size, 'z', i, save_path=full_path)

            # Save all XZ-plane slices (iterating along y-axis)
            total_slices = data.shape[1]
            print(f"Saving {total_slices} XZ-plane slices...")
            for i in range(total_slices):
                filename = f"{name}_XZ_y={i}.jpg"
                full_path = os.path.join(base_plots_folder, name, 'XZ', filename)
                visualize_phantom_slice(data, voxel_size, 'y', i, save_path=full_path)

            # Save all YZ-plane slices (iterating along x-axis)
            total_slices = data.shape[0]
            print(f"Saving {total_slices} YZ-plane slices...")
            for i in range(total_slices):
                filename = f"{name}_YZ_x={i}.jpg"
                full_path = os.path.join(base_plots_folder, name, 'YZ', filename)
                visualize_phantom_slice(data, voxel_size, 'x', i, save_path=full_path)
        print("\nAll slices have been saved successfully!\n")
    else:
        print("The program will continue without saving slice images.\n")

def create_vox_file (phantom_vox_file_name, phantom, voxel_size, material_array, density_array):
    """
    Creates a .vox voxel phantom file in the format required by the PENELOPE/penEasy 
    Monte Carlo simulation framework. First, a fixed header is written. Then, the 
    material and density data for each voxel are appended.

    Args:
        phantom_vox_file_name (str): The name of the output .vox file.
        phantom (np.ndarray): The 3D phantom array (used for shape info).
        voxel_size (np.ndarray): The physical dimensions of the voxels in mm.
        material_array (np.ndarray): 1D array of material IDs for each voxel.
        density_array (np.ndarray): 1D array of densities (g/cm^3) for each voxel.
    """
    
    # Using 'with open' for safer file handling. This ensures the file is
    # always closed properly, even if errors occur.
    # Create the .vox file
    with open(phantom_vox_file_name, 'w') as f:
        # VOX file header
        f.write('[SECTION VOXELS HEADER v.2008-04-13]\n')
        f.write(f'  {np.shape(phantom)[0]} {np.shape(phantom)[1]} {np.shape(phantom)[2]} \n')
        f.write(f'  {voxel_size[0]*0.1:.6f} {voxel_size[1]*0.1:.6f} {voxel_size[2]*0.1:.6f} \n')
        f.write(' 1\n')
        f.write(' 2\n')
        f.write(' 0\n')
        f.write('[END OF VXH SECTION]\n')

    # Combine material and density arrays to use np.savetxt(), which is very efficient.
    combined_array_mat_den = np.column_stack((material_array, density_array))
    with open(phantom_vox_file_name,'a') as f1:
        np.savetxt(f1, combined_array_mat_den, fmt = '%3d %7.4f')

def calculate_vbb_com(organ_df, phantom, phantom_ct, phantom_avgden, voxel_size):
    """
    Calculates the Voxel Bounding Box (VBB) and Center of Mass (COM) for each
    organ.
    
    Args:
        organ_df (pd.DataFrame): dataframe with organ information.
        phantom (np.ndarray): 3D array of organ IDs.
        phantom_ct (np.ndarray or None): 3D array of CT densities.
        phantom_avgden (np.ndarray): 3D array of average organ densities.
        voxel_size (np.ndarray): physical dimensions of voxels.
    
    Returns:
        organ_df(pd.DataFrame): updated organ_df with VBB and COM data.
    """
    print("Calculating voxel bounding boxes and center of mass for all organs...\n")
    
    # Initialize lists to store the VBB data
    vbb_data = {
        'x_min_vox': np.nan, 'x_max_vox': np.nan, 'y_min_vox': np.nan, 'y_max_vox': np.nan, 'z_min_vox': np.nan, 'z_max_vox': np.nan
    }

    # find_objects returns a list of tuples(min,max) of the min and max values for each organ in phantom array, for the X,Y and Z arrays
    # index i represents organ i+1. Organ_ID=0 is not needed because air outside phantom is the background
    vbb_slices = ndimage.find_objects(phantom)
    
    # Add the new columns to the DataFrame
    for key, value in vbb_data.items():
        organ_df[key] = value
    
    # Fill organ_df will results
    for i,slices in enumerate(vbb_slices):
        organ_id = i+1
        if slices is not None and organ_id in organ_df.index:
            organ_df.loc[organ_id, 'x_min_vox'] = slices[0].start
            organ_df.loc[organ_id, 'x_max_vox'] = slices[0].stop - 1
            organ_df.loc[organ_id, 'y_min_vox'] = slices[1].start
            organ_df.loc[organ_id, 'y_max_vox'] = slices[1].stop - 1
            organ_df.loc[organ_id, 'z_min_vox'] = slices[2].start
            organ_df.loc[organ_id, 'z_max_vox'] = slices[2].stop - 1

    # Calculate x/y/z size and values in cm
    organ_df['x_size_vox'] = organ_df['x_max_vox'] - organ_df['x_min_vox'] + 1
    organ_df['y_size_vox'] = organ_df['y_max_vox'] - organ_df['y_min_vox'] + 1
    organ_df['z_size_vox'] = organ_df['z_max_vox'] - organ_df['z_min_vox'] + 1

    organ_df['x_min_cm'] = organ_df['x_min_vox'] * voxel_size[0] * 0.1
    organ_df['x_max_cm'] = organ_df['x_max_vox'] * voxel_size[0] * 0.1
    organ_df['x_size_cm'] = organ_df['x_size_vox'] * voxel_size[0] * 0.1
    organ_df['y_min_cm'] = organ_df['y_min_vox'] * voxel_size[1] * 0.1
    organ_df['y_max_cm'] = organ_df['y_max_vox'] * voxel_size[1] * 0.1
    organ_df['y_size_cm'] = organ_df['y_size_vox'] * voxel_size[1] * 0.1 # Corrected from x_size_vox
    organ_df['z_min_cm'] = organ_df['z_min_vox'] * voxel_size[2] * 0.1
    organ_df['z_max_cm'] = organ_df['z_max_vox'] * voxel_size[2] * 0.1
    organ_df['z_size_cm'] = organ_df['z_size_vox'] * voxel_size[2] * 0.1 # Corrected from x_size_vox
    
    
    #### Center of Mass (com) Calculation ####

    # Create a mask of the organs that were actually segmented and are therefore actually present in the phantom
    # Get a list of organ IDs that are actually present in the phantom
    organs_in_phantom_mask = organ_df["Voxel_number"] > 0
    organ_ids_present = organ_df[organs_in_phantom_mask]['Organ_ID'].tolist()

    # Calculate com using phantom_avgden
    # Since each organ has unform density across all voxels, it's the same as the geometric center (centroid)
    print("Calculating center of mass for 'phantom_avgden'...\n")
    organ_df['com_x_vox_avg'] = np.nan
    organ_df['com_y_vox_avg'] = np.nan
    organ_df['com_z_vox_avg'] = np.nan

    centers_of_mass_avg = ndimage.center_of_mass(phantom_avgden, labels=phantom, index=organ_ids_present)
    
    # calculate organ_df with the com results from phantom_avgden
    for organ_id, com_coords in zip(organ_ids_present, centers_of_mass_avg):
        organ_df.loc[organ_id, 'com_x_vox_avg'] = com_coords[0]
        organ_df.loc[organ_id, 'com_y_vox_avg'] = com_coords[1]
        organ_df.loc[organ_id, 'com_z_vox_avg'] = com_coords[2]
        
    # Calculate com coordinates in cm for phantom_avgden
    organ_df['com_x_cm_avg'] = organ_df['com_x_vox_avg'] * voxel_size[0] * 0.1
    organ_df['com_y_cm_avg'] = organ_df['com_y_vox_avg'] * voxel_size[1] * 0.1
    organ_df['com_z_cm_avg'] = organ_df['com_z_vox_avg'] * voxel_size[2] * 0.1

    # calculate organ_df with the com results using phantom_ct, if it exists
    if phantom_ct is not None:
        print("Calculating center of mass for 'phantom_ct'...\n")
        organ_df['com_x_vox_ct'] = np.nan
        organ_df['com_y_vox_ct'] = np.nan
        organ_df['com_z_vox_ct'] = np.nan

        centers_of_mass_ct = ndimage.center_of_mass(phantom_ct, labels=phantom, index=organ_ids_present)
    
        # Populate the DataFrame with the com results from phantom_ct
        for organ_id, com_coords in zip(organ_ids_present, centers_of_mass_ct):
            organ_df.loc[organ_id, 'com_x_vox_ct'] = com_coords[0]
            organ_df.loc[organ_id, 'com_y_vox_ct'] = com_coords[1]
            organ_df.loc[organ_id, 'com_z_vox_ct'] = com_coords[2]
            
        # Calculate com coordinates in cm for phantom_ct
        organ_df['com_x_cm_ct'] = organ_df['com_x_vox_ct'] * voxel_size[0] * 0.1
        organ_df['com_y_cm_ct'] = organ_df['com_y_vox_ct'] * voxel_size[1] * 0.1
        organ_df['com_z_cm_ct'] = organ_df['com_z_vox_ct'] * voxel_size[2] * 0.1
    else:
        pass


    return organ_df

def create_simulation_files(phantom, voxel_size, phantom_avgden, phantom_mat_id, phantom_ct):
    """Asks the user if they want to create .vox files and creates them if requested.
    
    Args:
        phantom (np.ndarray): The 3D phantom array (used for shape info).
        voxel_size (np.ndarray): The physical dimensions of the voxels in mm.
        phantom_avgden (np.ndarray): 3D array containing average density in the voxels, calculated on average per each organ in the phantom.
        phantom_mat_id (np.ndarray): 3D array containing Material_ID for each voxel.
        phantom_ct (np.ndarray): 3D array containing exact density (g/cm^3) for each voxel, calculated from CT image.
    
    """
    choice = prompt_yes_no("Do you want to create .vox files for simulation?")

    if choice == 'y':
        # Create .VOX files, one with avg densities and another with per-voxel densities (if available)
        print('Creating the .VOX files...\n')
        phantom_vox_file_names=(DEFAULT_VOX_AVG_FILENAME,DEFAULT_VOX_CT_FILENAME)
        
        # Flatten arrays in 'F' (Fortran) order to match simulation software expectations
        phantom_avgden_flat = phantom_avgden.flatten(order='F')
        phantom_mat_id_flat = phantom_mat_id.flatten(order='F')
        
        # Create phantom with average organ densities
        create_vox_file(phantom_vox_file_names[0], phantom, voxel_size, phantom_mat_id_flat, phantom_avgden_flat)
        print(f"{phantom_vox_file_names[0]} created: represents the average density for all voxels comprising the organ.\n")

        
        # Only create the per-voxel density phantom if the CT data was provided
        if phantom_ct is not None:
            phantom_ct_flat = phantom_ct.flatten(order='F')
            # Create phantom with singular voxel densities
            create_vox_file(phantom_vox_file_names[1], phantom, voxel_size, phantom_mat_id_flat, phantom_ct_flat)
            print(f"{phantom_vox_file_names[1]} created: represents the exact density for each voxel from the CT image.")
        else:
            print(f"'{phantom_vox_file_names[1]}' was not created because no CT data was provided.")
    else:
        print("The .vox files will not be created.\n")


def main():
    """ 
    Main function to run the patient-specific phantom generation pipeline.
    
    This script orchestrates the entire process, including:
    1. Loading a segmented phantom labelmap from an NRRD file.
    2. Loading a corresponding organlist from a CSV file.
    3. Cropping the phantom to remove surrounding empty space.
    4. Optionally loading a CT NRRD file, converting its HU values to density,
       and calculating organ-specific statistics.
    5. Performing data integrity checks and creating phantom arrays based on
       average organ densities or per-voxel CT densities.
    6. Generating and saving 2D slice visualizations of the phantoms.
    7. Calculating Voxel Bounding Box (VBB) and Center of Mass (COM) for all organs.
    8. Optionally creating .vox files formatted for Monte Carlo simulations.
    9. Generating a comprehensive 'report.out' file with all calculated data tables.
    
    Returns:
        dict: A dictionary containing key variables/dfs for inspection in an interactive environment.
    """
    
    #### #### GET START TIME #### ####
    start_time = time.time()
    #### #### #### #### #### #### #### 

    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### LOAD PHANTOM DATA AND HEADER FROM NRRD FILE TO 3D-ARRAY
    
    print("\n")
    print(format_header("Welcome! ", width=14,border_char=">>>> "))
    print("This program is a pipeline to convert NRRD labelmap to a 3D voxel phantom, for visualization and simulation using Monte Carlo software.\n")
    
    
    ## get nrrd phantom file name from user, "phantom.nrrd" is the default
    nrrd_file = get_file_path_from_user("Enter the name of nrrd file containing the segmented phantom (labelmap):", DEFAULT_NRRD_FILE)
    
    print(f"\nLoading {nrrd_file}...\n")
    try:
        ph_data, header = nrrd.read(nrrd_file, index_order='F')
    except FileNotFoundError:
        print(f"Error: The file '{nrrd_file}' was not found. Verify the path provided.")
        sys.exit()

    voxel_size = np.sum(header["space directions"],axis=0)
    voxel_volume = voxel_size[0]*voxel_size[1]*voxel_size[2]
    
    print(">> Indices of the loaded nrrd phantom data:")
    print(f">>>> X-axis:  -> min:  0   -> max: {ph_data.shape[0]:3d} ")
    print(f">>>> Y-axis:  -> min:  0   -> max: {ph_data.shape[1]:3d} ")
    print(f">>>> Z-axis:  -> min:  0   -> max: {ph_data.shape[2]:3d} \n\n")
    
    report_content = []
    report_content.append(format_header("Welcome!", width=14, border_char = ">>>> "))
    report_content.append(" \n\n\n")
    report_content.append("This program is a pipeline to convert NRRD labelmap to a 3D voxel phantom, for visualization and simulation using Monte Carlo software. \n\n")
    report_content.append(f">> Number of voxels in {nrrd_file}: {ph_data.size} \n")
    report_content.append(f">> Number of missing values in {nrrd_file}: {np.count_nonzero(np.isnan(ph_data))} \n \n")   
    report_content.append(">> Indices of the segmented phantom data:\n")
    report_content.append(f">>>> X-axis:  -> min:  0   -> max: {ph_data.shape[0]:3d} \n")
    report_content.append(f">>>> Y-axis:  -> min:  0   -> max: {ph_data.shape[1]:3d} \n")
    report_content.append(f">>>> Z-axis:  -> min:  0   -> max: {ph_data.shape[2]:3d} \n\n")

    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### LOAD ORGANLIST FILE AND ORGAN DATABASE
    
    # get organlist file from user, "organlist.csv" is default
    organlist_file = get_file_path_from_user("Enter the name of the organlist file", DEFAULT_ORGANLIST_FILE)
    try:
        organ_df = pd.read_csv(organlist_file, sep=None, engine='python')
        # Validate organlist file, chack if required columns are present
        required_cols = ['Organ_ID', 'Organ', 'Material_ID', 'Material', 'Density_ICRP']
        missing_cols = [col for col in required_cols if col not in organ_df.columns]
        # List the columns missing from the organlist file
        if missing_cols:
            print(f"Error: The organlist file '{organlist_file}' is missing required columns: {', '.join(missing_cols)}")
            sys.exit()
    except FileNotFoundError:
        print(f"Error: The file '{organlist_file}' was not found. Verify the path provided.")
        sys.exit()

    organ_number = len(organ_df)

    # Count voxels per organ
    print("\nCropping phantom data to delete dead space and counting voxels per organ... \n")
    voxel_number = np.bincount(ph_data.ravel(), minlength=organ_number) # order='F', here for bincount doesn't matter
    organ_df["Voxel_number"] = voxel_number # Assign voxel_number to dataframe column

    # Cut the phantom to delete dead space (e.g. air, machine equipment) around the patient
    print("Automatically calculating cropping boundaries for the phantom...\n")
    phantom, crop_indices = auto_crop_phantom(ph_data, padding=2)
    # phantom contains organ IDs and is cropped. According to convention, air outside body = 0 and air inside body = organ_id.max()
    crop_x_min, crop_x_max, crop_y_min, crop_y_max, crop_z_min, crop_z_max = crop_indices
    print(f"Number of voxels in cropped phantom: {phantom.size}\n")
    
    report_content.append("The phantom data was automatically cropped to delete dead space (e.g. air, CT equipment).\n\n")
    report_content.append(f">> Number of voxels in cropped phantom: {phantom.size} \n\n")
    report_content.append(">> Cropping indices for the segmented phantom data:\n")
    print(">> Indices of the cropped phantom:")
    # Check if the phantom was cropeed in the X-axis and perform analogous verification for the other axis
    if crop_x_min == 0 and crop_x_max == ph_data.shape[0]: 
        print(">>>> X-axis: Phantom was not cropped in the X-axis. ")
        report_content.append(">>>> X-axis: Phantom was not cropped in the X-axis. \n")
    else: 
        print(f">>>> X-axis:  -> min: {crop_x_min:3d}  -> max: {crop_x_max - 1:3d} ")
        report_content.append(f">>>> X-axis: -> min: {crop_x_min} -> max: {crop_x_max - 1} \n")
    if crop_y_min == 0 and crop_y_max == ph_data.shape[1]:
        print(">>>> Y-axis: Phantom was not cropped in the Y-axis. ")
        report_content.append(">>>> Y-axis: Phantom was not cropped in the Y-axis. \n")
    else: 
        print(f">>>> Y-axis:  -> min: {crop_y_min:3d}  -> max: {crop_y_max - 1:3d} ")
        report_content.append(f">>>> Y-axis:  -> min: {crop_y_min}  -> max: {crop_y_max - 1} \n")
    if crop_z_min == 0 and crop_z_max == ph_data.shape[2]:
        print(">>>> Z-axis: Phantom was not cropped in the Z-axis. \n")
        report_content.append(">>>> Z-axis: Phantom was not cropped in the Z-axis. \n\n")
    else: 
        print(f">>>> Z-axis:  -> min: {crop_z_min:3d}  -> max: {crop_z_max - 1:3d} \n\n")
        report_content.append(f">>>> Z-axis:  -> min: {crop_z_min}  -> max: {crop_z_max - 1} \n\n")

    # Add phantom dimensions to report content
    report_content.append(">> Indices of the new cropped phantom:\n")
    report_content.append(f">>>> X-axis:  -> min:  0   -> max: {phantom.shape[0]:3d} \n")
    report_content.append(f">>>> Y-axis:  -> min:  0   -> max: {phantom.shape[1]:3d} \n")
    report_content.append(f">>>> Z-axis:  -> min:  0   -> max: {phantom.shape[2]:3d} \n\n")
    report_content.append(f">> Voxel size in x: {voxel_size[0] * 0.1:.5f} cm \n") 
    report_content.append(f">> Voxel size in y: {voxel_size[1] * 0.1:.5f} cm \n") 
    report_content.append(f">> Voxel size in z: {voxel_size[2] * 0.1:.5f} cm \n") 
    report_content.append(f">> Voxel volume   : {voxel_volume / 1000:0.5f} cm^3 \n\n")

    phantom_size_x = phantom.shape[0] * voxel_size[0]
    phantom_size_y = phantom.shape[1] * voxel_size[1]
    phantom_size_z = phantom.shape[2] * voxel_size[2]
    
    report_content.append(f">> Phantom size in x: {phantom_size_x * 0.1:.3f} cm \n") 
    report_content.append(f">> Phantom size in y: {phantom_size_y * 0.1:.3f} cm \n") 
    report_content.append(f">> Phantom size in z: {phantom_size_z * 0.1:.3f} cm \n") 
    report_content.append(f">> Phantom volume: {(phantom_size_x*phantom_size_y*phantom_size_z)/1000:.3f} cm^3 \n\n") 
    
    # Get number of voxels per organ on phantom after cutting (verify if cutting did not alter phantom)
    voxel_number_phantom = np.bincount(phantom.ravel(), minlength=organ_number)

    # Calculate the relative difference between the voxel organ number in data and phantom arrays
    # Use np.divide to handle division by zero safely
    voxel_number_quotient = np.divide(voxel_number, voxel_number_phantom, out=np.zeros_like(voxel_number, dtype=float), where=voxel_number_phantom!=0)
    
    # calculate volume in each organ and add to organ_df
    organ_df["Volume"] = organ_df["Voxel_number"] * voxel_volume / 1000.0  # organ Volume is in cm^3, voxel volume is in mm^3

    if np.allclose(voxel_number_quotient[1:40], 1):
        report_content.append(f"No differences were found between the number of voxels for the segmented structures in the nrrd\
phantom loaded and in the cropped phantom, except for the air outside phantom, which is {(voxel_number_quotient[0]*100):1.2f}% \
times higher in the labelmap.\n\n\n")
    
    print(f"Finished loading and cropping {nrrd_file}.\n")
    

    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### LOAD CT FILE
    
    # Initialize CT nrrd file, ct dataframe and calculate statistics.
    
    # Ask if user has CT file that corresponds to the segmented nrrd file
    phantom_ct = None # Initialize phantom_ct to None
    choice = prompt_yes_no("Do you have the CT nrrd file matching the segmented nrrd file?")

    if choice == 'y':
        # get ct file name from user, default is "ct_pre.nrrd"
        ct_nrrd_file = get_file_path_from_user("Enter the name of the CT nrrd file", DEFAULT_CT_FILE)
        print(f"\nLoading {ct_nrrd_file}... \n")
        
        try:
            ct_data, ct_header = nrrd.read(ct_nrrd_file, index_order='F')
        except FileNotFoundError:
            print(f"Error: The file '{ct_nrrd_file}' was not found. Verify the file provided.")
            sys.exit()
        
        # Cut ct_data to match the phantom dimensions
        ct_hu = ct_data[crop_x_min:crop_x_max, crop_y_min:crop_y_max, crop_z_min:crop_z_max]
        
        print(f"Loaded the {ct_nrrd_file} and cropped it to the phantom size.\n")
        report_content.append(f"Loaded the {ct_nrrd_file} and cropped it to the phantom size.\n\n")
        
        # --- Ask user for the CT HU calibration curve parameters ---
        print("Please provide the CT Hounsfield Units (HU) calibration curve parameters for the formula: HU = a * Density - b \n")
        while True:
            try:
                # Prompt user for the slope (a) and intercept (b)
                
                #hu_slope_a = float(input(">> Enter the slope 'a': "))
                #hu_intercept_b = float(input(">> Enter the intercept 'b': "))
                
                hu_slope_a = 902.3
                hu_intercept_b = 906.38
        
                # The slope 'a' cannot be zero
                if hu_slope_a == 0:
                    print("Error: The slope 'a' cannot be zero. Please enter a valid number.")
                    continue # Ask again
        
                break # Exit the loop if inputs are valid
        
            except ValueError:
                print("Invalid input. Please enter valid numbers for the slope and intercept.")
        # ----------------------------------------------------------------

        # Convert from HU to density using the calibration curve: HU = hu_slope_a * density - hu_intercept_b
        # phantom_ct now has densities
        phantom_ct = (ct_hu + hu_intercept_b) / hu_slope_a
        
        report_content.append(f"Used the calibration curve to convert Hounsfield Unit (HU) to density: HU = {hu_slope_a} * density - {hu_intercept_b} \n\n")
        
        # Correct air densities and handle negative values
        # phantom_ct now has normalized density values
        phantom_ct = _correct_ct_densities(phantom_ct, phantom, organ_df)

        # Calculate statistics directly from the CT data
        ct_df = calculate_ct_statistics(phantom, phantom_ct, organ_df, voxel_volume)
        
        # organ_df needs the newly calculated average densities
        organ_df["DensityAvg_CT"] = ct_df["Density_Mean"]

        #calculate organ mass from voxel-specific densities
        organ_df["Mass_CTDens(g)"] = (ct_df["Density_Mean"] * organ_df["Volume"])
        report_content.append("Organ mass was calculated for each Organ ID according to CT Density.\n\n")
        print("\nOrgan mass was calculated for each Organ ID according to CT Density.\n")
        print(f"Finished processing {ct_nrrd_file}.\n")
        print("The 'phantom_ct' was created. 'phantom_ct' contains the exact density for each voxel, obtained from the CT image and calibration curve.\n")
        report_content.append("The 'phantom_ct' was created. 'phantom_ct' contains the exact density for each voxel, obtained from the CT image and calibration curve.\n")

    else: # This block runs if choice is 'n'
        print("No CT nrrd file provided. The program will use ICRP reference densities instead.\n")
        # Fallback to using ICRP densities if no CT data is available
        organ_df["DensityAvg_CT"] = organ_df["Density_ICRP"]

    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### DATA INTEGRITY VERIFICATIONS, CALCULATIONS AND CREATION OF DATABASES
    
    print("Performing data integrity verifications and calculations...\n")

    # The rest of the script now uses organ_df["DensityAvg_CT"], which either uses data from the CT file
    # or from the ICRP reference values.
    # Check if there are negative density values -> safeguard
    if (organ_df['DensityAvg_CT'] < 0).values.any():
        for l in organ_df.index:
            if (organ_df.loc[l,"DensityAvg_CT"] <= 0):
                print(f"Organ ID {l} for {organ_df.loc[l,'Organ']} has negative average density of {organ_df.loc[l,'DensityAvg_CT']:.4f} for {organ_df.loc[l,'Organ']:.4f}.")
                print("This value is out of range. The ICRP density will be used as a reference value.\n")
                organ_df.loc[l,"DensityAvg_CT"] = organ_df.loc[l,"Density_ICRP"]
    
    if phantom_ct is not None:
        # If relative difference between DensityAvg_CT and Density_ICRP is too big (e.g.>75%), use ICRP reference density
        # Preserve original data for DensityAvg_CT before any corrections
        organ_df['DensityAvg_CT_old'] = organ_df['DensityAvg_CT']
        
        # Calculate relative difference between DensityAvg_CT and ICRP_Density.
        # do not multiply by 100, keep as a quotient because tabulate will transform the value to a percentage posteriorly 
        organ_df["RD_Dens_ICRP_CT"] = ( (organ_df["Density_ICRP"] - organ_df["DensityAvg_CT"]) / organ_df["Density_ICRP"] )
        
        # Use abs() to catch large deviations in both positive and negative directions.
        large_diff_mask = abs(organ_df["RD_Dens_ICRP_CT"]) > (DENS_MAX_REL_DIFF / 100.0)
        number_organs_corrected = large_diff_mask.sum()
        
        if number_organs_corrected > 0:
            print(f"Warning: Found {number_organs_corrected} organ(s) with a density deviation > 90%.")
            print(f"The organ(s) are: {organ_df.loc[large_diff_mask,'Organ']}")
            print("Their densities will be replaced with ICRP reference values.\n")
            # For the rows where the difference is too large, update DensityAvg_CT with the ICRP value
            organ_df.loc[large_diff_mask, 'DensityAvg_CT'] = organ_df.loc[large_diff_mask, 'Density_ICRP']
            
        # Calculate organ mass from organ DensityAvg_CT or Density_ICRP, if CT image was not provided
        organ_df["Mass_DensCT"] = organ_df["DensityAvg_CT"] * organ_df["Voxel_number"] * (voxel_volume/1000)
    else:
        organ_df["Mass"] = organ_df["Density_ICRP"] * organ_df["Voxel_number"] * (voxel_volume/1000)

    # Create avg density phantom and material phantom for posterior cleaner visualizations and file writing
    phantom_avgden = phantom.astype(np.float32) #initialize phantom_avgden
    phantom_mat_id = phantom.astype(np.float32)
    for l in organ_df.index: 
        phantom_avgden[phantom==l] = organ_df.loc[l,'DensityAvg_CT']
        phantom_mat_id[phantom==l] = organ_df.loc[l,'Material_ID']
    
    if phantom_ct is not None:
        print("The 'phantom_avgden' was created. For each organ, the density was calculated for all voxels comprising the organ. Then, \
the calculated average density was assigned to all voxels in the organ.\n\n")
        report_content.append("The 'phantom_avgden' was created. For each organ, the density was calculated for all voxels comprising \
the organ. Then, the calculated average density was assigned to all voxels in the organ.\n\n\n\n")
    else:
        print("The 'phantom_avgden' was created. For each organ, the ICRP reference density value was assigned to all voxels comprising \
the organ. Therefore, all organs have the ICRP reference density value as the average density. \n\n")
        report_content.append("The 'phantom_avgden' was created. For each organ, the ICRP reference density value was assigned to all \
voxels comprising the organ. Therefore, all organs have the ICRP reference density value as the average density. ICRP reference values \
were employed because no CT file or calibration curve were provided. \n\n\n\n")

    # Create the Material Database
    # The material database will be different if a ct file was not given
    print("Creating material database...\n")
    
    # if the is no column called "Mass_DensCT", no ct file was added
    if "Mass_DensCT" not in organ_df.columns:
        organ_df["Mass_DensCT"] = organ_df["Mass"]
        
    # Use groupby().agg() to perform multiple aggruptations directly
    mat_df = organ_df.groupby('Material_ID').agg(
        Material=('Material', 'first'),  # Get the material name (it's the same for all Organs composed by the same material)
        Voxel_number=('Voxel_number', 'sum'),  # Sum the voxels for each Organ_ID composed by the same material
        Total_mass=('Mass_DensCT', 'sum'),  # Calculate sum of mass of each organ
        Total_volume=('Volume','sum'),  # Calculate the sum of the volume of each organ
        Density_ICRP=('Density_ICRP', 'first')  # Get the ICRP density (it's the same for all Organs composed by the same material)
    ).reset_index()
    
    # Calculate the wighted average density for each material
    mat_df["Avg_Mat_Density"] = np.divide(mat_df["Total_mass"], mat_df["Total_volume"], out=np.zeros_like(mat_df['Total_mass']), 
                                        where=mat_df['Total_volume']!=0)
    
    
    if phantom_ct is not None:
        # define RD_Den as the relative diff between the ICRP 110 den and the CT calculated den
        # do not multiply by 100, keep as a quotient because tabulate will transform the value to a percentage posteriorly 
        mat_df['RD_Dens_ICRP_CT'] = ( (mat_df['Density_ICRP'] - mat_df['Avg_Mat_Density']) / mat_df['Density_ICRP'] )
    

    report_end_time(start_time)
    
    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### CREATE PHANTOM VISUALIZATION AND SAVE TO FILE

    # Create dictionary of existing phantoms to process
    phantoms_to_process = {
        'phantom_avgden': phantom_avgden,
        'phantom_ct': phantom_ct  # This will be None if the user answered 'n'
    }
    
    # Run function to save all slices in batch.
    save_all_phantom_slices(phantoms_to_process, voxel_size)
    
    report_end_time(start_time)
    
    
    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### CREATE VBB AND COM DATABASE AND WRITE TABLE TO REPORT



    # Add Voxels Bounding Box (VBB) to organ_db
    organ_df = calculate_vbb_com(organ_df, phantom, phantom_ct, phantom_avgden, voxel_size)

    # use tabulate to print organ_df
    organ_df_table_cols = ['Organ_ID','Organ','Material_ID','Material','Voxel_number', 'Volume','Density_ICRP']
    organ_table_fmt = ['d','','d','','d', '0.4f','0.4f']
    
    if phantom_ct is not None:
        # insert phantom_ct columns, i.e. density, RD and mass, to the list
        organ_df_table_cols.extend(['DensityAvg_CT', 'RD_Dens_ICRP_CT', 'Mass_DensCT'])
        organ_table_fmt.extend(['0.4f','2.1%','0.4f'])
    else:
        # insert Mass, calculated from Density_ICRP to the list
        organ_df_table_cols.append('Mass')
        organ_table_fmt.append('0.4f')
    
    organ_table_fmt = list(organ_table_fmt) # convert to list to insert in tabulate
    
    organ_table = tabulate(organ_df[organ_df_table_cols], tablefmt = 'simple', showindex = "never",
                           headers = organ_df_table_cols, 
                           floatfmt = organ_table_fmt, 
                           numalign ='right', stralign ='center')

    # Add organ_ID table to report content
    report_content.append(format_header("ORGAN DATABASE", width=20, border_char = ">>>> "))
    report_content.append('\n\n')
    report_content.append(organ_table)
    report_content.append("\n\n\n Volume: Calculated volume for each organ/structure (cm^3).\n")
    report_content.append(" Density_ICRP: Reference density value for each organ/structure according to the ICRP (g/cm^3).\n")
    report_content.append(" DensityAvg_CT: Average density value for each organ/structure calculated from the CT image and calibration curve (if provided), otherwise equal to Density_ICRP (g/cm^3).\n")
    if phantom_ct is not None:
        report_content.append(" RD_Dens_ICRP_CT: Relative difference between Density_ICRP and DensityAvg_CT.\n")
    report_content.append(" Mass_DensCT: Organ/structure mass calculated using DensityAvg_CT (g).\n")
    report_content.append('\n\n\n')

    # use tabulate to print mat_df
    mat_df_table_cols = ['Material_ID','Material','Voxel_number', 'Total_mass', 'Total_volume', 'Density_ICRP', 'Avg_Mat_Density']
    mat_table_fmt = ['d','','d','0.4f','0.4f','0.4f','0.4f']
    
    if phantom_ct is not None:
        mat_df_table_cols.extend(['RD_Dens_ICRP_CT'])
        mat_table_fmt.append('2.1%')
    else:
        pass
    
    mat_table_fmt = list(mat_table_fmt) # convert to list to insert in tabulate    

    mat_table = tabulate(mat_df[mat_df_table_cols], tablefmt = 'simple', showindex = "never",
                         headers = mat_df_table_cols,
                         floatfmt = mat_table_fmt, 
                         numalign = 'right', stralign = 'center')

    # Add mat_ID table to report content
    report_content.append("\n\n\n")
    report_content.append(format_header("MATERIAL DATABASE", width=20, border_char = ">>>> "))
    report_content.append('\n\n')
    report_content.append(mat_table)
    report_content.append("\n\n\n Avg_Mat_Density: Average density value calculated for material comprised from various organs/structures (g/cm^3).\n")
    if phantom_ct is not None:
        report_content.append(" RD_Dens_ICRP_CT: Relative difference between Density_ICRP and DensityAvg_CT.\n")
    report_content.append('\n\n\n')
    
    if phantom_ct is not None:
        # use tabulate to print ct_df
        ct_table = tabulate(ct_df, tablefmt = 'simple', showindex = "never",
                            headers = ct_df.columns.tolist(),
                            floatfmt = ['d', '', '0.4f', '0.4f', '0.4f'], 
                            numalign = 'right', stralign = 'center')
    
        # Add mat_ID table to report content
        report_content.append("\n\n\n")
        report_content.append(format_header("CT DATABASE", width=20, border_char = ">>>> "))
        report_content.append('\n\n')
        report_content.append(ct_table)
        report_content.append("\n\n\n Note: Density values are reported in g/cm^3.\n")
        report_content.append('\n\n\n')

    vbb_cols = ['Organ_ID','Organ','x_min_vox','x_max_vox','x_size_vox','x_min_cm','x_max_cm','x_size_cm',
                'y_min_vox','y_max_vox','y_size_vox','y_min_cm','y_max_cm','y_size_cm',
                'z_min_vox','z_max_vox','z_size_vox','z_min_cm','z_max_cm','z_size_cm']

    # use tabulate to print VBB from organ_df
    vbb_table = tabulate(organ_df[vbb_cols].fillna(0), tablefmt='simple', showindex="never", # Use fillna for pretty printing
                        headers = organ_df[vbb_cols].columns.tolist(),
                        floatfmt= ['d', '', '.0f','.0f','.0f','0.4f','0.4f','0.4f',
                        '.0f','.0f','.0f','0.4f','0.4f','0.4f',
                        '.0f','.0f','.0f','0.4f','0.4f','0.4f'], 
                        numalign='right', stralign='center')

    # Add VBB table to report content
    report_content.append("\n\n\n")
    report_content.append(format_header("VOXELS BOUNDING BOX (VBB) DATABASE", width=20, border_char = ">>>> "))
    report_content.append('\n\n')
    report_content.append(vbb_table)
    report_content.append("\n\n\n Note: To understand the phantom system of coordinates, visualize the phantom along the XY-,XZ- and YZ-planes. \n")
    report_content.append('\n\n\n')
    
    com_avgden_cols = ['Organ_ID','Organ','com_x_vox_avg','com_y_vox_avg','com_z_vox_avg',
                       'com_x_cm_avg','com_y_cm_avg','com_z_cm_avg']

    # use tabulate to print COM from organ_df, as calculated from phantom_avgden
    com_avg_table = tabulate(organ_df[com_avgden_cols].fillna(0), tablefmt='simple', showindex="never",
                        headers = ('Organ_ID','Organ','com_x_vox','com_y_vox','com_z_vox', 'com_x_cm','com_y_cm','com_z_cm'),
                        floatfmt= ['d', '', '.0f','.0f','.0f','0.4f','0.4f','0.4f'], 
                        numalign='right', stralign='center')

    # Add COM table to report content
    report_content.append("\n\n\n")
    report_content.append(format_header("CENTER OF MASS (COM) DATABASE - PHANTOM AVGDEN", width=20, border_char = ">>>> "))
    report_content.append('\n\n')
    report_content.append(com_avg_table)
    report_content.append('\n\n\n')

    
    if phantom_ct is not None:
        com_ct_cols = ['Organ_ID','Organ','com_x_vox_ct','com_y_vox_ct','com_z_vox_ct', 
                       'com_x_cm_ct','com_y_cm_ct','com_z_cm_ct']
        
        # use tabulate to print COM from organ_df, as calculated from phantom_ct
        com_ct_table = tabulate(organ_df[com_ct_cols].fillna(0), tablefmt='simple', showindex="never",
                            headers = ('Organ_ID','Organ','com_x_vox','com_y_vox','com_z_vox', 'com_x_cm','com_y_cm','com_z_cm'),
                            floatfmt= ['d', '', '.0f','.0f','.0f','0.4f','0.4f','0.4f'], 
                            numalign='right', stralign='center')
    
        # Add COM table to report content
        report_content.append("\n\n\n")
        report_content.append(format_header("CENTER OF MASS (COM) DATABASE - PHANTOM CT", width=20, border_char = ">>>> "))
        report_content.append('\n\n')
        report_content.append(com_ct_table)
        report_content.append('\n\n\n')
    


    
    report_end_time(start_time)
    
    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### CREATE SIMULATION FILE
    
    # Call function to handle .vox file creation.
    create_simulation_files(phantom, voxel_size, phantom_avgden, phantom_mat_id, phantom_ct)
    
    #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### FINAL ADJUSTMENTS
    
    # Pack strategic variables into a dictionary to return for inspection in Spyder's variable explorer
    # Use .get('variable_name') for variables that may not exist (i.e. CT-related ones)
    results = {'organ_df': organ_df,
                    'mat_df': mat_df,
                    'ct_df': locals().get('ct_df', None),
                     'ph_data': ph_data,
                    'phantom': phantom,
                    'phantom_avgden': phantom_avgden,
                    'phantom_ct': phantom_ct,
                    'phantom_mat_id': phantom_mat_id,
                    'voxel_size': voxel_size,
                    'voxel_volume': voxel_volume,
                    'crop_indices': crop_indices,
                    'header': header,
                    'hu_slope_a': locals().get('hu_slope_a', None),
                    'hu_intercept_b': locals().get('hu_intercept_b', None),
    }
    
    report_content.append('\n\n\n\n\n\n')
    report_content.append(format_header("Thank your for using this pipeline! Come back soon!", width=14, border_char = ">>>> "))
    report_content.append('\n')
    
    # open "report.out" as f1 and write the entire report file from report_content list
    with open('report.out', 'w') as f1:
        f1.writelines(report_content)
    
    report_end_time(start_time)
    
    print('Report on segmented phantom and CT nrrd files written to report.out.\n\n\n')

    print('\n')
    print(format_header('Thank you for using this pipeline! Come Back soon!', width=14, border_char = ">>>> "))

    return results


if __name__ == "__main__":
    results_dict = main()
# unpack dictionary to the global namespace
if results_dict:
    globals().update(results_dict)
