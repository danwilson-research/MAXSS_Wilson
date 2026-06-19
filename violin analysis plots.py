# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:49:04 2026

@author: dw557
"""

#Script to create summary violin plots of the toal flux for storms in 2010 and 2011

# Import required packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set the colours to be used in these plots
# Define colorblind palette for cb_colors indices
cb_colors = sns.color_palette("colorblind", 11)

# Map the NEW labels (the ones in your label_mapping) to the colors
final_color_map = {
    'Full Simulation': 'k',                    # MAXSS_RUN
    'Pre-storm Ref': cb_colors[1],             # REF_RUN
    'Wind': cb_colors[7],                      # WIND_RUN
    'SST (with gradients)': cb_colors[9],       # SST_WITH_GRADIENTS_RUN (Plot 1 name)
    'SST (with Gradients)': cb_colors[9],       # SST_WITH_GRADIENTS_RUN (Plot 2 name)
    'SSS': cb_colors[4],                       # SSS_RUN
    'Salinity (SSS)': cb_colors[4],            # SSS_RUN (Plot 2 name)
    'Pressure': cb_colors[5],                  # PRESSURE_RUN
    'XCO2': cb_colors[2]                       # V_GAS_RUN
}

## Plot 1: Comparison of total flux in each model component run ##

# Load in required data
df = pd.read_csv('E:/MAXSS_working_directory/HPC_output/Spatially_integrated_fluxes/storm_component_flux_summary_2010_to_2011.csv')

# Keep 'Year' in id_vars so we can split by rows
df_melted = df.drop(columns=['Region']).melt(id_vars=['Storm', 'Year'], 
                                             var_name='Component Run', 
                                             value_name='Total Flux (Tg C)')

label_mapping = {
    'MAXSS_RUN' : 'Full Simulation', 
    'REF_RUN': 'Pre-storm Ref',
    'WIND_RUN' : 'Wind',
    'SST_WITH_GRADIENTS_RUN': 'SST (with gradients)',
    'SSS_RUN' : 'SSS',
    'V_GAS_RUN': 'XCO2',
    'PRESSURE_RUN': 'Pressure'}
                                                          
# Apply the mapping to the column
df_melted['Component Run'] = df_melted['Component Run'].map(label_mapping)        

# Set up the vertical stacked grid
g = sns.catplot(
    x='Component Run', 
    y='Total Flux (Tg C)', 
    row='Year',               
    data=df_melted, 
    kind='violin', 
    hue='Component Run',
    palette=final_color_map, 
    inner='points', 
    cut=0, 
    alpha=0.75,
    legend=False,
    height=4.5,               
    aspect=2.5                
)

# --- ADDED THIS LINE HERE ---
# Forcefully strip out Seaborn's default row/column title text templates
g.set_titles(row_template="", col_template="")

# Panel label lookup tracking index numbers to clean alphabetical tags
panel_labels = {0: 'a)', 1: 'b)'}

# Apply formatting across all vertical axes panels safely
for i, ax in enumerate(g.axes.flat):
    ax.grid(zorder=1, linestyle='--', alpha=0.5)
    ax.set_xlabel('Component Run', fontweight='bold')
    ax.set_ylabel('Total Flux (Tg C)', fontweight='bold')
    
    # Inject just 'a)' or 'b)' left-aligned without any interference
    ax.set_title(panel_labels[i], fontweight='bold', loc='left', fontsize=14)

plt.tight_layout()
plt.show()