# Project
Visulize_WGMS_Annual_MB_TM_2025
This repo is a package to visulize the annual MB  data from WGMS2025

# Usage

## Generate the global glacier mass change in grid 
python Map_Global_Glacier_MB.py --filepath "path/to/input.nc" --outputpath "path/to/output" --variablename "glacier_mass_change_gt" --outputgif "glacier_mass_change.gif" --titlename "Mass change (Gt)"

## Generate the unobserved glacier info as a df
Python Table_Glacier_Individual.py --regioncode "ALA" --regionname "Alaska" --path_individual "path_to_your_directory"

# Output
The output figures and Tables will be in the Folder Figures and Tables correspondingly

# File structure
Visulize_WGMS_Annual_MB_TM_2025/
│
├── Map_Global_Glacier_MB.py
├── Table_Glacier_individual.py
├── Figures/
│   ├── glacier_mass_change.gif
└── Tables/
    ├── Alaska_unobserved_statistics.csv

# License
MIT License
