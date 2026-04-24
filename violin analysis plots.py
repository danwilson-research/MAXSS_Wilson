# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:49:04 2026

@author: dw557
"""

#Script to create summary violin plots

#Import required packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Set the colours to be used in these plots
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
df = pd.read_csv('E:/MAXSS_working_directory/output/spatially_integrated_fluxes/storm_component_flux_summary.csv')

# Turn the run columns into a single column named 'Component Run'
# Drop the year and region columns
df_melted = df.drop(columns=['Year', 'Region']).melt(id_vars=['Storm'], 
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
                                                          
#Apply the mapping to the column
df_melted['Component Run'] = df_melted['Component Run'].map(label_mapping)        

# 3. Create the Plot
plt.figure(figsize=(12, 6))

#Set up the violin plots
sns.violinplot(x='Component Run', y='Total Flux (Tg C)', data=df_melted, legend = False, 
               hue='Component Run',palette=final_color_map, inner='points', cut=0, zorder=2, alpha = 0.75) # 'inner=points' shows each storm's value

# Set Bold Labels
plt.xlabel('Component Run', fontweight='bold')
plt.ylabel('Total Flux (Tg C)', fontweight='bold')

#add grid and xticks
plt.grid(zorder =1)
#plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()


## Plot 2 Comparison of Taylor Flux Results ##

df = pd.read_csv('E:/MAXSS_working_directory/output/taylor_decomposition_summary/taylor_decomposition_summary.csv')

# Turn the run columns into a single column named 'Taylor decomposition'
# Drop the year and region columns
df_melted = df.drop(columns=['Year', 'Region', 'SST_NoGrad_Contribution_TgC','Actual_Anomaly_TgC',
'Taylor_Sum_NoGrad_TgC', 'Taylor_Sum_WithGrad_TgC']).melt(id_vars=['Storm_ID'], 
                                                     var_name='Taylor decomposition', 
                                                     value_name='Total Flux (Tg C)')
 
#Customise plto labels
label_mapping = {
    'Wind_Contribution_TgC': 'Wind',
    'SSS_Contribution_TgC': 'Salinity (SSS)',
    'VGas_Contribution_TgC': 'XCO2',
    'Pressure_Contribution_TgC': 'Pressure',
    'SST_WithGrad_Contribution_TgC': 'SST (with Gradients)'}
                                                          
#Apply the mapping to the column
df_melted['Taylor decomposition'] = df_melted['Taylor decomposition'].map(label_mapping)                                                          
                                                          
                                                          
# 3. Create the Plot
plt.figure(figsize=(12, 6))

#Set up the violin plots
sns.violinplot(x='Taylor decomposition', y='Total Flux (Tg C)', data=df_melted, legend = False, 
               hue='Taylor decomposition', palette=final_color_map, inner='points', cut=0, zorder=2, alpha = 0.75) # 'inner=points' shows each storm's value

# Set Bold Labels
plt.xlabel('Taylor decomposition', fontweight='bold')
plt.ylabel('Total Flux (Tg C)', fontweight='bold')

#add grid and xticks
plt.grid(zorder =1)
#plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()