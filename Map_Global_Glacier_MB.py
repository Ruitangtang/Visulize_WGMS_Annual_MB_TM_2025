# Script Name: Map_Global_Glacier_MB.py
# Author: Ruitang Yang, tangruiyang123@gmail.com
# Description: This script maps global glacier mass balance (Grid, unit is m w.e. a-1 or Gt a-1).
# Date Created: 2025-08-17
# Version: 1.0

# Import necessary libraries
import argparse
import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
from cmcrameri import cm  #colormap from crameri
# Import cmcrameri (install : # Using conda :conda install cmcrameri; or using pip: pip install cmcrameri)


#%% Define functions for mapping glacier mass balance

# Function to read the nc file (read the attributes and variables)
def read_nc_file(file_path):
    """
    Reads a netCDF file and prints its attributes and variables.

    --- parameters ---
    file_path: str
        The path to the netCDF file to read.

    --- returns ---
    dataset: netCDF4.Dataset
        The opened netCDF dataset.
    variable_names: list
        A list of variable names present in the dataset.
    """
    try:
        dataset = nc.Dataset(file_path, 'r')
        print(f"Successfully opened netCDF file: {file_path}")
        
        # Print the global attributes
        print("Global attributes:")
        for attr in dataset.ncattrs():
            print(f"{attr}: {getattr(dataset, attr)}")

        # Print the variable names and their details,
        print("\nVariables:")
        variable_names = []
        for var_name in dataset.variables:
            var = dataset.variables[var_name]
            print(f"{var_name}: Dimensions: {var.dimensions}, Shape: {var.shape}, Data Type: {var.dtype}, Variable type: {type(var)}, Units: {getattr(var, 'units', 'N/A')}")
            variable_names.append(var_name)

        # Example of reading a variable (replace 'variable_name' with your variable)
        if 'variable_name' in dataset.variables:
            data = dataset.variables['variable_name'][:]
            print("\nData for 'variable_name':", data)

        return dataset, variable_names

    except Exception as e:
        print(f"Error reading netCDF file: {e}")
        return None, None


# Function to extract variable data
def extract_variable_data(dataset, variable_name):
    """
    Extracts data for a specific variable from the dataset.

    --- parameters ---
    dataset: netCDF4.Dataset
        The netCDF dataset to extract data from.
    variable_name: str
        The name of the variable to extract.

    --- returns ---
    var_data: netCDF4.Variable
        The extracted variable data.
    """
    try:
        var_data = dataset.variables[variable_name]
        return var_data
    except Exception as e:
        print(f"Error extracting {variable_name} data: {e}")
        return None


# Function to calculate IQR
def calculate_iqr(data):
    """
    Calculate the IQR and adjusted limits for the provided data, ignoring NaNs.
    --- Parameters ---
    data: np.ndarray
        The input data array for which to calculate the IQR.
    --- Returns ---
    q1: float
        The first quartile (25th percentile) of the data.
    q3: float
        The third quartile (75th percentile) of the data.
    iqr: float
        The interquartile range (IQR) of the data.
    """
    q1 = np.nanpercentile(data, 25)  # First quartile (25th percentile)
    q3 = np.nanpercentile(data, 75)  # Third quartile (75th percentile)
    iqr = q3 - q1  # Interquartile range
    return q1, q3, iqr


# Function to create global grid map GIF
def create_global_grid_map_gif(file_path=None, variable_name=None, title_n=None, output_gif=None, output_path=None):
    """
    Create an animated GIF from NetCDF data with consistent colormap scaling based on IQR.
    
    Parameters:
    - file_path: Path to NetCDF file (required)
    - variable_name: Variable to visualize (required)
    - title_n: Title for plot (required)
    - output_gif: Output file path (required)
    - output_path: Directory to save the output GIF (required)

    Returns:
    - Path to created GIF file, and all intermediate frames
    """
    # Validate inputs
    if None in (file_path, variable_name, title_n, output_gif, output_path):
        raise ValueError("All parameters must be provided")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # Ensure proper file extension
    if not output_gif.lower().endswith('.gif'):
        output_gif += '.gif'
    
    try:
        # Open NetCDF file
        with nc.Dataset(file_path, mode='r') as dataset:
            # Load coordinates
            lat = dataset.variables['lat'][:]
            lon = dataset.variables['lon'][:]
            # Load all data at once for consistent scaling
            all_data = np.ma.filled(dataset.variables[variable_name][:], np.nan)

            # Calculate global IQR for color scaling
            global_q1, global_q3, global_iqr = calculate_iqr(all_data)
            global_vmin = global_q1 - 1.5 * global_iqr
            global_vmax = global_q3 + 1.5 * global_iqr
            
            print(f"Creating frames with consistent color scaling based on global IQR...")
            print(f"Color scale range: {global_vmin:.2f} to {global_vmax:.2f}")

            # Prepare figure
            fig = plt.figure(figsize=(12, 6))
            frames = []
            num_time_slices = all_data.shape[0]  # Determine the number of time slices
            
            for time_slice in range(num_time_slices):
                data = all_data[time_slice]

                if np.isnan(data).all():
                    print(f"Skipping empty frame {time_slice}")
                    continue

                # Clear the current plot
                plt.clf()
                ax = plt.axes(projection=ccrs.PlateCarree())
                
                # Create plot with consistent colors
                im = ax.pcolormesh(lon, lat, data,
                                   cmap=cm.vik_r,
                                   shading='auto',
                                   transform=ccrs.PlateCarree(),
                                   vmin=global_vmin,
                                   vmax=global_vmax)
                
                # Add title and features
                year = 1976 + time_slice
                ax.set_title(f"{title_n} - Year: {year}", fontsize=14)
                ax.coastlines()
                
                # Configure grid
                # ax.set_xticks(np.arange(-180, 181, 30), crs=ccrs.PlateCarree())
                # ax.set_yticks(np.arange(-90, 91, 30), crs=ccrs.PlateCarree())
                ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)

                # Add colorbar
                cbar = plt.colorbar(im, ax=ax, pad=0.05)
                cbar.set_label(title_n)
                cbar.set_ticks(np.linspace(global_vmin, global_vmax, 5))

                # Render and store frame
                fig.canvas.draw()
                frame = np.array(fig.canvas.renderer.buffer_rgba())
                frames.append(frame)
                print(f"Created frame {time_slice + 1}/{num_time_slices}")
                # save the immediate figure
                plt.savefig(f"{output_path}/frame_{year:03d}.png", bbox_inches='tight')
                plt.close(fig)

            # Verify we have enough frames
            if len(frames) < 2:
                raise ValueError(f"Only {len(frames)} valid frames - need at least 2 for animation")
            
            # Convert frames to uint8
            frames_uint8 = [frame.astype(np.uint8) for frame in frames]
            
            # Save as animated GIF
            output_gif = os.path.join(output_path, output_gif)
            print(f"Saving {len(frames_uint8)} frames to {output_gif}...")
            imageio.mimsave(
                output_gif,
                frames_uint8,
                format='GIF',
                duration=0.5,  # 0.5 seconds per frame
                loop=0        # Infinite looping
            )
            
            # Verify output
            if not os.path.exists(output_gif):
                raise RuntimeError("Output file was not created")
            if os.path.getsize(output_gif) < 1024:
                print("Warning: Output file is unusually small")
            
            print(f"Successfully created animated GIF with {len(frames)} frames")
            return output_gif
        # Close the dataset
        dataset.close()
    
    except Exception as e:
        print(f"Error during GIF creation: {str(e)}")
        raise
    finally:
        plt.close('all')  # Close all figures to free up memory


#%% Main function
def main():

    # Step 0: get parser
    parser = argparse.ArgumentParser(description="Global Glacier Mass Balance Mapping")
    parser.add_argument("--filepath", type=str, required=True, help="Path to the input netCDF file")
    parser.add_argument("--outputpath", type=str, required=True, help="Path to the output directory")
    parser.add_argument("--variablename", type=str, default="glacier_mass_change_gt", help="Variable name to extract from the netCDF file")
    parser.add_argument("--outputgif", type=str, default="glacier_mass_change.gif", help="Name of the output GIF file")
    parser.add_argument("--titlename", type=str, default="Mass change (Gt)", help="Title for the plot")
    args = parser.parse_args()

    # Step 1: Read the netCDF file
    file_path = args.filepath
    dataset, variable_names = read_nc_file(file_path)

    if dataset is None:
        print("Failed to read netCDF file.")
        return

    # Step 2: Extract 
    mb_data = extract_variable_data(dataset, "mass_balance")

    if mb_data is None:
        print("Failed to extract mass balance data.")
        return

    # Step 3: Close the dataset
    dataset.close()

    # Step 4: Plot the glacier mass balance data
    # plot the timeseries gif
    create_global_grid_map_gif(file_path = file_path,variable_name = args.variablename,
                               title_n = args.titlename, output_gif=args.outputgif,
                                output_path=args.outputpath)


# Main code
if __name__ == "__main__":
    print("Global Glacier Mass Balance Mapping Script (Version 1.0)")
    main()

# Example usage
# python Map_Global_Glacier_MB.py --filepath "path/to/input.nc" --outputpath "path/to/output" --variablename "glacier_mass_change_gt" --outputgif "glacier_mass_change.gif" --titlename "Mass change (Gt)"
