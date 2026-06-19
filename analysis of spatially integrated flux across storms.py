# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:33:37 2026

@author: dw557
"""

# This code is used to analyse the results of TC flux calculations run using FluxEngine
# as part of the MAXSS project.

####
# Plot the total flux across all storms and all years as panel plots
####

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
    'SSS': cb_colors[4],                       # SSS_RUN
    'Salinity (SSS)': cb_colors[4],            # SSS_RUN (Plot 2 name)
    'Pressure': cb_colors[5],                  # PRESSURE_RUN
    'XCO2': cb_colors[2]                       # V_GAS_RUN
}

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

####
# Compute Taylor decomposition results and save to .csv
####

# Import required packages
import pandas as pd
import numpy as np

# 1. Load the original raw un-melted data
file_path = 'E:/MAXSS_working_directory/HPC_output/Spatially_integrated_fluxes/storm_component_flux_summary_2010_to_2011.csv'
df = pd.read_csv(file_path)

# 2. Define the baseline reference column
baseline_col = 'REF_RUN'
total_sim_col = 'MAXSS_RUN'

# 3. Identify the separate components to evaluate
components = {
    'WIND_RUN': 'Wind Contribution',
    'SST_WITH_GRADIENTS_RUN': 'SST Contribution',
    'SSS_RUN': 'SSS Contribution',
    'PRESSURE_RUN': 'Pressure Contribution',
    'V_GAS_RUN': 'XCO2 Contribution'
}

# 4. Compute the actual total anomaly seen in the full model simulation
df['Total_Anomaly'] = df[total_sim_col] - df[baseline_col]

# 5. Taylor Decomposition: Calculate linear first-order anomaly contributions
decomposition_results = df[['Storm', 'Year', 'Region', 'Total_Anomaly']].copy()

linear_sum = 0
for run_col, label in components.items():
    # Anomaly contribution = Component Run minus the Pre-storm Baseline Reference
    decomposition_results[label] = df[run_col] - df[baseline_col]
    linear_sum += decomposition_results[label]

# 6. Calculate the Taylor Residual (Non-linear interaction term)
# Residual = Total True Anomaly - Sum of linear approximations
decomposition_results['Residual'] = decomposition_results['Total_Anomaly'] - linear_sum

# 7. Calculate Relative Contribution (%) over all storms combined

# This calculates the average contribution from each component to the total anomaly across all storms

print("=== GLOBAL MEAN TAYLOR DECOMPOSITION CONTRIBUTIONS ===")
mean_total_anomaly = decomposition_results['Total_Anomaly'].abs().mean()

for label in list(components.values()) + ['Residual']:
    # Measure contribution magnitude relative to total anomaly
    mean_comp = decomposition_results[label].mean()
    relative_pct = (decomposition_results[label].abs().mean() / mean_total_anomaly) * 100
    print(f"{label:<25} | Mean Value: {mean_comp:6.3f} Tg C | Relative Absolute Magnitude: {relative_pct:5.1f}%")

# Save results for downstream plotting or tables
decomposition_results.to_csv('E:/MAXSS_working_directory/HPC_output/Spatially_integrated_fluxes/taylor_decomposition_results.csv', index=False)

####
# Plot the Taylor decomposition results
####

# Script to create horizontal summary violin plots stacked vertically by year
# Blue if median is negative (Sink), Orange if positive (Source) 

# Import required packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load the decomposition results
df_decomp = pd.read_csv('E:/MAXSS_working_directory/HPC_output/Spatially_integrated_fluxes/taylor_decomposition_results.csv')

# Map the technical column names to short names
column_mapping = {
    'SST Contribution': 'SST',
    'XCO2 Contribution': 'XCO2',
    'Wind Contribution': 'Wind',
    'Pressure Contribution': 'Pressure',
    'SSS Contribution': 'SSS',
    'Residual': 'Residual'
}

# Melt the data into long-form
value_cols = list(column_mapping.keys())
df_melted = df_decomp.melt(
    id_vars=['Storm', 'Year'], 
    value_vars=value_cols,
    var_name='Component', 
    value_name='Flux Anomaly (Tg C)'
)
df_melted['Component'] = df_melted['Component'].map(column_mapping)

# Structural hierarchy: True physical variables first, Residual at the bottom
component_order = ['SST', 'XCO2', 'Wind', 'Pressure', 'SSS', 'Residual']

# Define the explicit binary colors
COLOR_SINK = '#2b7bba'    # Blue
COLOR_SOURCE = '#d95f02'  # Orange

# 2. Build the vertical multi-panel plot
g = sns.catplot(
    x='Flux Anomaly (Tg C)',    
    y='Component',              
    row='Year',
    data=df_melted,
    kind='violin',              
    order=component_order,
    hue='Component',            
    inner='points',             
    cut=0,                      
    alpha=0.75,
    height=4.5,
    aspect=2.2,
    sharex=True,
    legend=False                
)

# Strip default "Year = 2010" row labels
g.set_titles(row_template="", col_template="")

panel_labels = {0: 'a) 2010', 1: 'b) 2011'}
years = [2010, 2011]

# 3. Apply structural styling and binary coloring based on median
for i, ax in enumerate(g.axes.flat):
    current_year = years[i]
    
    # Isolate the current year's data subset
    df_year = df_melted[df_melted['Year'] == current_year]
    
    # Extract only the actual violin PolyCollections (filtering out inner point elements)
    violins = [c for c in ax.collections if isinstance(c, plt.matplotlib.collections.PolyCollection)]
    
    # Loop through the components in their exact plotting order
    for idx, comp in enumerate(component_order):
        # Calculate the median for this specific component and year
        comp_median = df_year[df_year['Component'] == comp]['Flux Anomaly (Tg C)'].median()
        
        # Determine color: Blue if median is negative (Sink), Orange if positive (Source)
        chosen_color = COLOR_SINK if comp_median < 0 else COLOR_SOURCE
        
        # Apply the color to the corresponding violin
        if idx < len(violins):
            violins[idx].set_facecolor(chosen_color)
            violins[idx].set_edgecolor(chosen_color)

    # Baseline grids and labels
    ax.axvline(0, color='black', linewidth=1, linestyle='-', zorder=3)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
            
    ax.set_ylabel('Driver Component', fontweight='bold')
    ax.set_title(panel_labels[i], fontweight='bold', loc='left', fontsize=13)
    
    if i == 1:
        ax.set_xlabel('Flux Anomaly (Tg C) \n[ Blue = Net Ocean Sink | Orange = Net Ocean Source ]', fontweight='bold')
    else:
        ax.set_xlabel('')

plt.tight_layout()
plt.show()

####
# Plot violin plots of total anomaly per year
####

# Script to look at the total net anomaly distribution split by year 
# with max sink, max source, mean lines, yearly cumulative impacts, and storm counts (N).

# Import required packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load the decomposition results
df_decomp = pd.read_csv('E:/MAXSS_working_directory/HPC_output/Spatially_integrated_fluxes/taylor_decomposition_results.csv')

# Pre-calculate the multi-year global total for reference if needed
global_total = df_decomp['Total_Anomaly'].sum()

# Define the explicit binary colors
COLOR_SINK = '#2b7bba'    # Blue
COLOR_SOURCE = '#d95f02'  # Orange

# 2. Build the vertical multi-panel plot for Total Anomaly
g = sns.catplot(
    x='Total_Anomaly',    
    row='Year',
    data=df_decomp,
    kind='violin',              
    inner='points',             # Shows individual storm dots cleanly
    cut=0,                      
    alpha=0.75,
    height=3.5,
    aspect=2.5,
    sharex=True                 
)

# Clear out default titles
g.set_titles(row_template="", col_template="")

panel_labels = {0: 'a) 2010 Total Anomaly Distribution', 1: 'b) 2011 Total Anomaly Distribution'}
years = [2010, 2011]

# 3. Structural styling, dynamic coloring, and annotations
for i, ax in enumerate(g.axes.flat):
    current_year = years[i]
    df_year = df_decomp[df_decomp['Year'] == current_year]
    
    # Calculate the cumulative net impact for just this specific year
    yearly_cumulative_impact = df_year['Total_Anomaly'].sum()
    
    # --- DYNAMIC COLOR LOGIC BASED ON CUMULATIVE IMPACT ---
    # Blue if net impact is negative (Sink), Orange if positive (Source)
    chosen_color = COLOR_SINK if yearly_cumulative_impact < 0 else COLOR_SOURCE
    
    # Apply color to the specific violin patch in this panel
    violins = [c for c in ax.collections if isinstance(c, plt.matplotlib.collections.PolyCollection)]
    for violin in violins:
        violin.set_facecolor(chosen_color)
        violin.set_edgecolor(chosen_color)
        
    # Baseline grid structures
    ax.axvline(0, color='black', linewidth=1.2, linestyle='-', zorder=3)
    ax.grid(axis='x', linestyle='--', alpha=0.5)
            
    ax.set_ylabel('')
    ax.set_yticks([])  
    ax.set_title(panel_labels[i], fontweight='bold', loc='left', fontsize=12)
    
    if not df_year.empty:
        # Calculate the number of storms in this specific year
        num_storms = df_year['Storm'].nunique()
        
        # --- STORM COUNT LABEL (Top Left Corner) ---
        ax.text(
            0.02, 0.93, f"n = {num_storms}",
            transform=ax.transAxes,
            color='#333333', fontweight='bold', fontsize=10,
            ha='left', va='top'
        )
        
        # Calculate the actual numerical mean dynamically
        yearly_mean = df_year['Total_Anomaly'].mean()
        
        # --- DRAW THE MEAN INDICATOR ---
        ax.axvline(yearly_mean, color='#e6194B', linewidth=1.5, linestyle='--', zorder=4)
        
        ax.text(
            yearly_mean, -0.2, f'Mean: {yearly_mean:.2f} Tg C',
            color='#e6194B', fontweight='bold', fontsize=9.5,
            ha='center', va='top'
        )
        
        # --- MAX SINK ANNOTATION (Left Side Outlier) ---
        max_sink = df_year['Total_Anomaly'].min()
        max_sink_storm = df_year.loc[df_year['Total_Anomaly'].idxmin(), 'Storm']
        
        ax.annotate(
            f'{max_sink_storm}\n({max_sink:.2f} Tg C)',
            xy=(max_sink, 0.02),
            xytext=(max_sink + 1.2, 0.25),
            arrowprops=dict(facecolor='black', shrink=0.08, width=0.8, headwidth=5),
            fontweight='bold',
            fontsize=9,
            color='#1a4a70'
        )
        
        # --- MAX SOURCE ANNOTATION (Right Side Outlier) ---
        max_source = df_year['Total_Anomaly'].max()
        max_source_storm = df_year.loc[df_year['Total_Anomaly'].idxmax(), 'Storm']
        
        if max_source > 0:
            ax.annotate(
                f'{max_source_storm}\n(+{max_source:.2f} Tg C)',
                xy=(max_source, 0.02),
                xytext=(max_source - 2.5, 0.25), 
                arrowprops=dict(facecolor='black', shrink=0.08, width=0.8, headwidth=5),
                fontweight='bold',
                fontsize=9,
                color='#d95f02' 
            )
            
        # --- YEARLY CUMULATIVE IMPACT OVERLAY (Bottom Left Box) ---
        text_string = f"Cumulative Net Impact ({current_year}):\n{yearly_cumulative_impact:.2f} Tg C"
            
        ax.text(
            0.02, 0.05, text_string, 
            transform=ax.transAxes,
            color='black', fontweight='bold', fontsize=9.5,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='gray', alpha=0.9),
            ha='left', va='bottom'
        )

    # Clean up X-axis text bounds
    if i == 1:
        ax.set_xlabel('Total Flux Anomaly per Storm (Tg C) \n[ Blue Panel = Cumulative Net Ocean Sink | Orange Panel = Cumulative Net Ocean Source ]', fontweight='bold')
    else:
        ax.set_xlabel('')

plt.tight_layout()
plt.show()