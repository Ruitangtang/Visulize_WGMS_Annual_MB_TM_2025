# Script Name: Table_Glacier_Individual.py
# Author: Ruitang Yang, tangruiyang123@gmail.com
# Description: This script analyzes glacier data for individual regions based on individual glacier characteristics.
# Date Created: 2025-08-17
# Version: 1.0

# Import necessary libraries
import pandas as pd
import os
import argparse


# Function to analyze unobserved glaciers at the subregion level for the specified region
def analyze_unobserved_glacier_data(region_code = None, region_name = None, path_individual = None):
    """
    Analyzes glacier data from a CSV file and calculates statistics for unobserved glaciers
    for each subregion.

    Parameters:
    - region_code: The region code used in the filename (string)
    - region_name: The name of the region (string)
    - path_individual: Path to the directory containing the CSV file (string)

    Returns:
    - DataFrame with calculated statistics
    """
    # Construct the filename
    MB_meta = f"{region_code}_{region_name}_metadata.csv"
    path_meta = os.path.join(path_individual, MB_meta)

    # Create a DataFrame from the CSV data
    df = pd.read_csv(path_meta, delimiter=',')
    
    # Print the head of the DataFrame for verification
    print(f"First few rows of the data from {MB_meta}:")
    print(df.head())

    # Extract the region ID
    df['RGI_region_ID(2nd)'] = df['RGIId'].str.split('-').str[1].str.split('.').str[1].str[:2]

    # Group by RGI_region_ID(2nd) and calculate required statistics
    result = df.groupby('RGI_region_ID(2nd)').agg(
        No_glacier=('Area', 'size'),
        Area=('Area', 'sum'),
        No_unobs=('N_gla_anom_used', lambda x: (x == 'no_obs').sum()),
    ).reset_index()

    # Calculate area percentage of unobserved glaciers
    result['Area_perc(unobs)'] = (result['No_unobs'] / result['No_glacier']) * 100
    
    # Add unit to the Area and Area_perc columns in the header
    result.rename(columns={
        'Area': 'Area (km²)',
        'Area_perc(unobs)': 'Area_perc (unobs, %)'
    }, inplace=True)

    # Round the Area_perc values to 2 decimal places
    result['Area_perc (unobs, %)'] = result['Area_perc (unobs, %)'].round(2)

    # Print the final statistics
    print(f"Statistic info about unobserved glaciers in each subregion in {region_name} ({region_code}):\n")
    print(result)

    # Save the results to a CSV file
    output_csv = os.path.join(path_individual, f"{region_name}_unobserved_statistics.csv")
    result.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    return result


#%% Main function
def main():

    # Step 0: get parser
    parser = argparse.ArgumentParser(description="Table Glacier Individual Analysis")
    parser.add_argument("--regioncode", type=str, required=True, help="Region code")
    parser.add_argument("--regionname", type=str, required=True, help="Region name")
    parser.add_argument("--path_individual", type=str, required=True, help="Path to the output directory")

    args = parser.parse_args()
    # Example usage of the function
    region_code = args.regioncode
    region_name = args.regionname
    path_individual = args.path_individual

    # Step 1: Call the analysis function of unobserved glaciers
    analyze_unobserved_glacier_data(region_code, region_name, path_individual)




if __name__ == "__main__":
    # Example usage:
    print("Table Glacier Individual Analysis Script (Version 1.0)")
    main()


# Example usage of the function
# Python Table_Glacier_Individual.py --regioncode "ALA" --regionname "Alaska" --path_individual "path_to_your_directory"