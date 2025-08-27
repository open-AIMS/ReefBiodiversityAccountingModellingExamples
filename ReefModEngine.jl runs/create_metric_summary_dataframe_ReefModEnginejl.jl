############################################ Description ################################################################
# Script for creating a biodiversity accounting summary data table as a csv from a set of ReefModEngine.jl results.

############################################ Load packages ##############################################################
using ReefModEngine
using CSV, DataFrames
using GeoDataFrames
using NCDatasets
using Statistics

###################################### Load data summary functions ######################################################

# Updates scenario dataframe with scenario info and metrics
function update_scenario_record!(
    display_scenario_df, cell_cover, sv_reefs, reef_area, iv_num_cells, scen_count
)
    total_int_area_cm = ((100.0)^2) .* iv_num_cells
    rel_cover_sp = dropdims(cell_cover; dims=2) ./ total_int_area_cm

    display_scenario_df[scen_count, "Coral Cover"] = mean(
        rel_cover_sp[:, :, 1]; dims=[1, 2]
    )[1]
    display_scenario_df[scen_count, "Coral Cover sd"] = std(
        rel_cover_sp[:, :, 1]; dims=[1, 2]
    )[1]

    simps_d_temp = simps_D(rel_cover_sp[:, :, 2:7])
    display_scenario_df[scen_count, "Diversity"] = mean(simps_d_temp; dims=[1, 2])[1]
    display_scenario_df[scen_count, "Diversity sd"] = std(simps_d_temp; dims=[1, 2])[1]

    sv_temp = shelter_vol(sv_reefs, cell_cover, iv_num_cells)
    display_scenario_df[scen_count, "Shelter Volume"] = mean(sv_temp; dims=[1, 2])[1]
    display_scenario_df[scen_count, "Shelter Volume sd"] = std(sv_temp; dims=[1, 2])[1]

    RBCI_temp = (rel_cover_sp[:, :, 1] .+ simps_d_temp .+ sv_temp) ./ 3
    display_scenario_df[scen_count, "RBCI"] = mean(rbci_temp; dims=[1, 2])[1]
    display_scenario_df[scen_count, "RBCI sd"] = std(RBCI_temp; dims=[1, 2])[1]

    display_scenario_df[scen_count, "Reef k area m2"] = reef_area

    return nothing
end

# Updates scenario record with intervention - counterfactual data
function update_scenario_record!(
    display_scenario_df, cell_cover, sv_reefs, iv_num_cells, scen_count
)
    total_int_area_cm = ((100.0)^2) .* iv_num_cells

    rel_cover_sp_iv = dropdims(cell_cover[:, :, :, :, 2]; dims=2) ./ total_int_area_cm
    rel_cover_sp_cf = dropdims(cell_cover[:, :, :, :, 1]; dims=2) ./ total_int_area_cm

    simps_d_temp_iv = simps_D(rel_cover_sp_iv[:, :, 2:7])
    simps_d_temp_cf = simps_D(rel_cover_sp_cf[:, :, 2:7])

    sv_temp_iv = shelter_vol(
        sv_reefs[:, :, :, :, 2], cell_cover[:, :, :, :, 2], iv_num_cells
    )
    sv_temp_cf = shelter_vol(
        sv_reefs[:, :, :, :, 1], cell_cover[:, :, :, :, 1], iv_num_cells
    )

    RBCI_temp_iv = (rel_cover_sp_iv .+ simps_d_temp_iv .+ sv_temp_iv) ./ 3
    RBCI_temp_cf = (rel_cover_sp_cf .+ simps_d_temp_cf .+ sv_temp_cf) ./ 3

    RBCI_uplift = RBCI_temp_iv .- RBCI_temp_cf
    display_scenario_df[scen_count, "RBCI uplift mean"] = mean(RBCI_uplift; dims=[1, 2])[1]
    display_scenario_df[scen_count, "RBCI uplift sd"] = std(RBCI_uplift; dims=[1, 2])[1]

    return nothing
end

###################################### Load metric calculation functions ################################################

# Shelter volume metric
function shelter_vol(sv_reefs, cell_cover, iv_num_cells)
    return sum(
        dropdims(sv_reefs; dims=2) .*
        dropdims(cell_cover[:, :, :, 2:end] .* 0.0001; dims=2);
        dims=3
    ) ./ (33 .* (iv_num_cells .* (pi * (0.5)^2)))
end

# Simpson's diversity metric
function simps_D(species_rel_cover)
    sum_cover = sum(species_rel_cover; dims=3)
    simps_d = 1 .- dropdims(sum((species_rel_cover ./ sum_cover) .^ 2; dims=3); dims=3)
    simps_d[isnan.(simps_d)] .= 0.0
    return simps_d
end

################################## Load function for creating data summary table ########################################

function create_metric_summary_dataframe_reefmod(cell_log_filepath,
    scenario_par_log_filepath, domain_spatial_df,
    iv_reef, aadpt_dhw_range, n_corals_range; rme_spatial_filepath="reefmod_gbr.gpkg",
    display_years=[2025, 2030, 2035, 2040, 2045, 2050, 2055, 2059],
    save_filename="RME_summary_file.csv")

    # Load resultset file (logged cell-level data)
    ds = NCDataset(cell_log_filepath)

    # Load scenario dataframe and yearly intervention parameters dataframe
    scenario_df = CSV.read(scenario_par_log_filepath, DataFrame)
    scenario_yearly_iv_df = CSV.read(
        string(result_folder, "\\iv_yearly_scenarios.csv"), DataFrame
    )

    # Load RME spatial data
    rme_spatial = GeoDataFrames.read(rme_spatial_filepath)
    longs = rme_spatial.LON
    lats = rme_spatial.LAT

    # Get intervention years
    iv_years = unique(scenario_yearly_iv_df.year)
    n_iv_years = length(iv_years) - 1

    reef_id_order = findall(in.(rme_spatial.GBRMPA_ID, Ref(domain_spatial_df.reef_id)))
    reef_names = rme_spatial.GBRMPA_ID[reef_id_order]
    n_reefs = length(reef_names)

    reef_areas_vec = domain_spatial_df[:, "area_km2"] .* domain_spatial_df[:, "k"] .* (10^6)

    # Load cell level data on cover, shelter volume, + number of iv cells
    cell_cover = ds["reef_cell_cover_cm"][:, :, :, :, :, :]
    sv_reefs = ds["species_sv_var"][:, :, :, :, :, :]
    iv_num_cells = ds["cell_indices_counts"][:, :]

    # Total span of years to include in the dataframe
    years = collect(2025:2099)[1:size(reef_cell_cover_cm, 2)]

    # Setup storage dataframe
    n_scens = (
        (
            length(display_years) * length(n_corals_range) * length(aadpt_dhw_range) *
            n_reefs
        ) + length(display_years) * length(n_corals_range) * n_reefs
    )
    display_scenario_df = DataFrame(
        hcat(fill("", n_scens), zeros(n_scens, 18)),
        ["Reef", "Year", "Intervention",
            "Deployment Volume", "Coral Cover", "Coral Cover sd", "Diversity",
            "Diversity sd", "Shelter Volume",
            "Shelter Volume sd", "RBCI", "RBCI sd", "RBCI uplift mean", "RBCI uplift sd",
            "Estimated intervention area m2",
            "Reef k area m2", "RBCI uplift area x estimated intervention area m2", "Lat",
            "Long"]
    )

    scen_count = 0

    for yr in display_years
        # Scenario year index
        yr_id = findall(years .== yr)
        for rfs in 1:n_reefs
            # Scenario reef index
            rf_idx = rme_spatial.GBRMPA_ID .== reef_names[rfs]
            for nc in n_corals_range

                # First, counterfactual scenario
                scen_count += 1
                display_scenario_df[scen_count, "Year"] = yr
                display_scenario_df[scen_count, "Intervention"] = 0
                display_scenario_df[scen_count, "Deployment Volume"] = 0
                display_scenario_df[scen_count, "Reef"] = reef_names[rfs]
                display_scenario_df[scen_count, "Long"] = longs[rf_idx]
                display_scenario_df[scen_count, "Lat"] = lats[rf_idx]

                # Coral scenario idx (needed as RME has paired counterfactual and intervention runs)
                coral_scens = (scenario_df.n_corals .== nc)

                # Add metrics to dataframe for CF
                update_scenario_record!(display_scenario_df,
                    cell_cover[findall(coral_scens), yr_id, :, rfs, :, 1],
                    sv_reefs[findall(coral_scens), yr_id, :, rfs, :, 1],
                    reef_areas_vec[rfs],
                    iv_num_cells[findall(coral_scens), :], scen_count)

                for aa in aadpt_dhw_range
                    # Assisted adaptation level ID
                    scen_ids = findall(coral_scens .& (scenario_df.a_adapt .== aa))

                    # Update data for intervention scenario
                    scen_count += 1
                    display_scenario_df[scen_count, "Year"] = yr
                    display_scenario_df[scen_count, "Intervention"] = aa
                    display_scenario_df[scen_count, "Deployment Volume"] = nc
                    display_scenario_df[scen_count, "Reef"] = reef_names[rfs]
                    display_scenario_df[scen_count, "Long"] = longs[rf_idx]
                    display_scenario_df[scen_count, "Lat"] = lats[rf_idx]

                    # If the reef is an intervention reef, calculate actual number of corals outplanted
                    if reef_names[rfs] == iv_reef
                        iv_id = findall(
                            (scenario_df.a_adapt .== aa) .& (scenario_df.n_corals .== nc)
                        )
                        if yr <= iv_years[end]
                            # If year is within the intervention time frame, add all corals outplanted up to this year
                            actual_n_corals =
                                ((
                                    mean(
                                        scenario_yearly_iv_df[
                                            in.(
                                                scenario_yearly_iv_df[:, "intervention id"],
                                                Ref(iv_id)
                                            ),
                                            "number of corals"
                                        ]
                                    ) / n_iv_years
                                )) * findfirst(iv_years .== yr)
                        else
                            actual_n_corals = ((mean(
                                scenario_yearly_iv_df[
                                    in.(
                                        scenario_yearly_iv_df[:, "intervention id"],
                                        Ref(iv_id)
                                    ),
                                    "number of corals"
                                ]
                            )))
                        end
                        # Estimate intervention area as number of corals/5 corals/m2
                        display_scenario_df[scen_count, "Estimated intervention area m2"] =
                            actual_n_corals / 5
                    end

                    # Update metrics for intervention scenario
                    update_scenario_record!(
                        display_scenario_df,
                        cell_cover[scen_ids, yr_id, :, rfs, :, 2],
                        sv_reefs[scen_ids, yr_id, :, rfs, :, 2],
                        reef_areas_vec[rfs],
                        iv_num_cells[scen_ids, :],
                        scen_count
                    )
                    # Add intervention minus counterfactual metrics
                    update_scenario_record!(
                        display_scenario_df,
                        cell_cover[scen_ids, yr_id, :, rfs, :, :],
                        sv_reefs[scen_ids, yr_id, :, rfs, :, :],
                        iv_num_cells[scen_ids, :],
                        scen_count
                    )
                end
            end
        end
    end

    # Add column with estimated intervention area * RBCI uplift
    display_scenario_df[:, "RBCI uplift area m2 x estimated intervention area"] .=
        display_scenario_df[:, "Estimated intervention area m2"] .*
        display_scenario_df[:, "RBCI uplift mean"]
    # Save as CSV
    CSV.write(save_filename, display_scenario_df)

    return display_scenario_df
end

################################### Key filenames from ReefModEngine.jl dataset #########################################

# Filepath for netcdf which has cell-level data saved
cell_log_filepath = "cell_log_2025marchmilestoneruns.nc"
# Filepath to CSV storing info on intervention scenario parameters
scenario_par_log_filepath = "scenario_par_log_2025milestoneruns.csv"
# Folder containing RME resultset
result_folder = "ext_moore_results_2025milestoneruns"
# Dataframe containing domain of interest, reef ids, reef area and k proportion
reef_areas_file = "C:\\Users\\rcrocker\\Documents\\Github\\ReefModEngine.jl\\sandbox\\ext_moore_areas.csv"
domain_spatial_df = CSV.read(reef_areas_file, DataFrame)

####################################### Key scenarios to include in dataframe ###########################################

# Scenario parameters to include
aadpt_dhw_range = [5.0, 10.0, 20.0] # assisted adaptation ranges
n_corals_range = [1000000] # deployment volume ranges
iv_reef = ["16-071"] # Intervention reef/reefs ID

################################## Create biodiversity metric summary dataframe #########################################

# Create and save metrics summary dataframe
summ_data_frame = create_metric_summary_dataframe_reefmod(cell_log_filepath,
    scenario_par_log_filepath, domain_spatial_df,
    iv_reef, aadpt_dhw_range, n_corals_range;
    save_filename="scenario_log_biodiversity_milestone_march2025_ExtMooreRMEruns.csv")
