############################################ Description ################################################################
# Script for using ReefBiodiversityAccountSetup with ADRIA data (selected impact sites from a set of ADRIA runs and
# an ADRIA domain geopackage) to suggest control sites for an intervention scenario.

############################################# Load packages #############################################################
using DataFrames, GeoDataFrames
using CSV
import ReefBiodiversityAccountSetup as RBAS

################################################ Load data ##############################################################
# Load a dataframe detailing ranked site ids selected as impact sites during a set of ADRIA runs
impact_sites_df = CSV.read(
    "path to CSV of saved site IDS in ranked order from a set of ADRIA runs")
# Load the ADRIA domain geopackage used for those runs as a geodataframe
site_data = GeoDataFrames.read(
    "path to ADRIA domain geopackage used for the set of ADRIA runs")

################################ Set variables for selecting control sites ##############################################
# Total number of control sites to select for each impact site
n_control = 3
# Number of impact sites to use in this scenario
n_impact = 5
# Columns of the geodataframe to use as criteria when selecting control sites
select_cols = ["site_id", "k_area", "habitat", "depth_med", "dhw_med"]
# Select impact site IDs from the loaded rankings dataframe
impact_sites = impact_sites_df.site_id[1:n_impact]
# Add k_area to the ADRIA domain dataframe to use as a criteria for control site selection
site_data[!, "k_area"] .= (site_data.k./100).*site_data.area

################################ Set variables for selecting control sites ##############################################
# Get indices of inpact site IDs in the ADRIA site data
impact_site_ids = findall(dropdims(any(site_data.site_id.==reshape(impact_sites,(1,n_impact)), dims=2), dims=2))
# Dataframe to store suggested control sites for each impact site
store_controls = DataFrame(fill("", n_control, length(impact_sites)), impact_sites)

########################################## Suggest control sites ########################################################
# Loop through impact sites and select n control sites for each impact site
for (imp_id_idx, imp_id) in enumerate(impact_sites)
    control_sites = RBAS.spatial_analysis.suggest_control_sites(findall(site_data.site_id .== imp_id)[1], site_data[:, select_cols], [:habitat], ID_COLUMN=:site_id)
    store_controls[1:minimum([n_control,size(control_sites,1)-1]), imp_id_idx] .= site_data.site_id[Int.(control_sites[1:minimum([n_control,size(control_sites,1)]), 1])]
end

# Save control site IDs as CSV
CSV.write("control_site_ids.csv", store_controls)
