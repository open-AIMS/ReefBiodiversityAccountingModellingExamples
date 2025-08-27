############################################ Description ################################################################
# Script for plotting and exploring set of data from the Allen Coral Atlas (https://www.allencoralatlas.org/) and NOAA
# Coral Reef Watch databases (https://coralreefwatch.noaa.gov) with the intention of selection a region for a
# reef biodiversity accounting project. Details on downloading the data are in the corresponding notebook
# plot_allen_atlas_data_example.ipynb.

############################################ Load packages ##############################################################
using DataFrames
import ReefBiodiversityAccountSetup as RBAS

############################################ Load config file ###########################################################
# config.toml includes info on where the relevant Allen Atlas and NOAA data is saved (see example config.toml file)
config_file = RBAS.spatial_analysis.load_config(; config_path="config.toml")

############################################ Plot spatial datasets ######################################################
# Plot gpkgs downloaded from the Allen Coral Atlas, including benthic and geomorphic categories.

# Load spatial datasets
benthic, geomorphic, reef_extent = RBAS.spatial_analysis.load_spatial_base(config_file)

# Set box limits for loading account data
box_upper = (-16.75284, 146.15641)
box_lower = (-16.95082, 146.32396)

# Get benthic and geomorphic data within the account box
benthic = RBAS.spatial_analysis.get_geo_within_box(benthic, box_upper, box_lower)
geomorphic = RBAS.spatial_analysis.get_geo_within_box(geomorphic, box_upper, box_lower)
extent = RBAS.spatial_analysis.get_geo_within_box(reef_extent, box_upper, box_lower)

# Plot benthic and geomorphic data
fig_benthic = RBAS.plotting.spatial_map(
    benthic,
    benthic[:, :class];
    opts=Dict(:legend_name => "Benthic map",
        :color_map => :Accent_5)
)
fig_geomorphic = RBAS.plotting.spatial_map(
    geomorphic,
    geomorphic[:, :class];
    opts=Dict(:legend_name => "Geomorphic map",
        :color_map => :hawaii10)
)
fig_extent = RBAS.plotting.spatial_map(
    extent,
    extent[:, :class];
    opts=Dict(:legend_name => "Extent map",
        :color_map => :hawaii10)
)

############################################ Analyse spatial datasets ###################################################
# Manipulate spatial datasets to understand coral real estate in the region.

# Filter for Coral/Algae type
benthic_filtered = RBAS.spatial_analysis.filter_site_area(benthic)

# Get intersection of geomorphic polygons and account extent
geomorphic_ext = RBAS.spatial_analysis.multipoly_geom_intersection(
    extent, geomorphic, :class
)

# Get intersection of benthic filtered polygons and account extent
benthic_ext = RBAS.spatial_analysis.multipoly_geom_intersection(
    extent, benthic_filtered, :class
)

# Get intersection of benthic filtered polygons and geomorphic polygons
geomorphic_benthic_comb = RBAS.spatial_analysis.multipoly_geom_intersection(
    benthic_ext, geomorphic_ext, :class
)

fig_geomorphic_ext = RBAS.plotting.spatial_map(
    geomorphic_ext,
    geomorphic_ext[:, :class];
    opts=Dict(:legend_name => "Geomorphic class",
        :color_map => :hawaii10)
)

fig_benthic_filtered = RBAS.plotting.spatial_map(
    benthic_ext,
    benthic_ext[:, :class];
    opts=Dict(:legend_name => "Coral/Rock substrate",
        :color_map => :tab10)
)

fig_geomorphic_filtered = RBAS.plotting.spatial_map(
    geomorphic_benthic_comb,
    geomorphic_benthic_comb[:, :class];
    opts=Dict(:legend_name => "Geomorphic categories of Coral/Rock substrate",
        :color_map => :hawaii10)
)

# Project and add polygon areas and k_areas to gdf
geomorphic_filtered = RBAS.spatial_analysis.set_reef_k(
    geomorphic_ext, benthic_ext)

# Plot area as spatial heat map
fig_k_filtered = RBAS.plotting.spatial_map(
    geomorphic_filtered,
    geomorphic_filtered[:, :k];
    opts=Dict(:colorbar_label => "k area",
        :color_map => :lighttest)
)

##################################### Integrate other environmental datasets #############################################
# Add other datasets to the spatial dataset, including depth, turbidity and dhw information aggregated over spatial polygons.

# Extract depths from raster file
geomorphic_filtered, depths = RBAS.spatial_analysis.median_features_allen(
    geomorphic_filtered, config_file; is_depth=true
)

# Plot depths
fig_depth = RBAS.plotting.spatial_map(
    geomorphic_filtered,
    geomorphic_filtered[:, :depth_med];
    opts=Dict(:colorbar_label => "Median depth",
        :color_map => :lighttest)
)

# Extract broadscale NOAA DHWs
geomorphic_filtered, dhws = RBAS.spatial_analysis.noaa_dhw_means(geomorphic_filtered, config_file)

# Plot mean and std for dhws
fig_dhw_mean = RBAS.plotting.spatial_map(
    geomorphic_filtered,
    geomorphic_filtered[:, :dhw_hist_mean];
    opts=Dict(:colorbar_label => "Mean max dhw",
        :color_map => :lighttest)
)
fig_dhw_std = RBAS.plotting.spatial_map(
    geomorphic_filtered,
    geomorphic_filtered[:, :dhw_hist_sd];
    opts=Dict(:colorbar_label => "Std max dhw",
        :color_map => :lighttest)
)

# Allen atlas turbidity data
turb_fn = config_file["other_data"]["allen_turbid"]

# Add turbidity data from Allen Atlas
geomorphic_filtered, turbidity = RBAS.spatial_analysis.median_features_allen(
    geomorphic_filtered, turb_fn; data_name=:turb_med
)

fig_turb_med = RBAS.plotting.spatial_map(
    geomorphic_filtered,
    geomorphic_filtered[:, :turb_med];
    opts=Dict(:colorbar_label => "Turbidity mean",
        :color_map => :lighttest)
)

############################### Evaluate potential impact and control sites in the region ################################
# Use the added spatial and environmental data as criteria to evaluate site suitability for a restoration project.
# Then evaluate suitable control sites for the suggested impact sites using similarity to specified criteria.

# Add k area as a selection criteria
geomorphic_filtered[!, "k_area"] = geomorphic_filtered[:,"k"].*geomorphic_filtered[:,"area"]

# Define which columns of the dataset to use as criteria
criteria = ["class", "depth_med", "dhw_hist_mean", "turb_med", "k_area"]
geomorphic_filtered[!, "site_id"] = collect(1:size(geomorphic_filtered, 1))

# Calculate a rating for each polygon to determine suitability according to the defined criteria.
# Here only sheltered reef slope type sites are considered
impact_site_rating = RBAS.spatial_analysis.suggest_impact_sites(geomorphic_filtered, criteria, "Sheltered Reef Slope")

# Plot the site ratings as a heat map
fig_impact= RBAS.plotting.spatial_map(
    geomorphic_filtered[Int.(impact_site_rating[:, :Sites]),:],
    impact_site_rating[:, :Rating];
    opts=Dict(:colorbar_label => "Impact site rating",
        :color_map => :lighttest)
)

# Suggest control sites for a particular impact site (e.g. 5), using the criteria in geomorphic_filtered's columns
# The input [:class] specifies that the selected control sites must be of the same geomorphic class as the impact site
control_site_list = RBAS.spatial_analysis.suggest_control_sites(5, geomorphic_filtered[:, Not(:geom)], [:class]; ID_COLUMN=:site_id)
