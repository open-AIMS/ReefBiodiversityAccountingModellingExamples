############################################ Description ################################################################
# Script for creating a biodiversity accounting summary data table as a csv from a set of ADRIA.jl results.

############################################ Load packages ##############################################################
using ADRIA
using CSV
using DataFrames, YAXArrays, GeoDataFrames
import ArchGDAL as ag
using GeoFormatTypes
using Revise, Infiltrator
using DimensionalData
import Statistics as st

###################################### Load metric calculation functions ################################################

# Relative shelter volume function
function relative_sv(rs::ADRIA.ResultSet)::YAXArray
    k_area = ADRIA.site_k_area(rs)
    SV_abs = ADRIA.metrics.absolute_shelter_volume(rs)
    max_cover::YAXArray = ADRIA.ZeroDataCube(
        (:timesteps, :species, :locations, :scenarios),
        (length(SV_abs.timesteps), 5*7,
            length(SV_abs.locations), length(SV_abs.scenarios))
    )
    max_cover[:, 7, :, :] .= 0.5 .* rs.loc_data.k'
    max_sv = ADRIA.metrics._absolute_shelter_volume(max_cover, k_area, rs.inputs)
    return YAXArray(SV_abs.axes, SV_abs.data ./ max_sv.data)
end

# Cover over each functional group (summed over size classes)
function species_cover(rs::ADRIA.ResultSet)::YAXArray
    total_taxa_cover = rs.outcomes[:relative_taxa_cover]
    n_timesteps, n_taxa, n_locs, n_scens = size(total_taxa_cover)
    total_species_cover::YAXArray = ADRIA.ZeroDataCube(
        (:timesteps, :species, :locations, :scenarios), (n_timesteps, 5, n_locs, n_scens)
    )
    for sp in 1:5
        total_species_cover[:, sp, :, :] .= Array(
            dropdims(sum(total_taxa_cover[:, (sp * 7 - 6):(sp * 7), :, :]; dims=2); dims=2)
        )
    end

    return total_species_cover
end

# Simpson's diversity function
function simps_D(rs::ADRIA.ResultSet)::YAXArray
    total_species_cover = species_cover(rs)
    sum_cover = sum(total_species_cover; dims=2)
    simps_d = dropdims(
        1 .- sum((total_species_cover ./ Array(sum_cover)) .^ 2; dims=2); dims=2
    )
    simps_d[isnan.(simps_d)] .= 0.0
    return simps_d
end

########################################### Load data summary functions #################################################

# Function to add means and stds of metrics to the table for a set of cf or intervention scenarios
function update_scenario_record!(display_scenario_df, rc, rd, sv, scen_count)
    idx_scens = scen_count:(scen_count + length(rc.locations) - 1)
    # Add selected site ids
    display_scenario_df[idx_scens, "Site"] .= rc.locations

    # Update individual metrics for selected sites
    display_scenario_df[idx_scens, "Coral Cover"] .= st.mean(rc.data; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "Coral Cover sd"] .= st.std(rc.data; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "Diversity"] .= st.mean(rd.data; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "Diversity sd"] .= st.std(rd.data; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "Shelter Volume"] .= st.mean(sv.data; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "Shelter Volume sd"] .= st.std(sv.data; dims=3)[1, :, 1]

    # Update RCI for selected sites
    display_scenario_df[idx_scens, "RCI"] .= st.mean(
        (rc.data + rd.data + sv.data) ./ 3; dims=3
    )[
        1, :, 1
    ]
    display_scenario_df[idx_scens, "RCI sd"] .= st.std(
        (rc.data + rd.data + sv.data) ./ 3; dims=3
    )[
        1, :, 1
    ]

    return nothing
end

# Function to add mean and std of difference to cf for metrics
function update_scenario_record!(
    display_scenario_df, rc_iv, rd_iv, sv_iv, rc_cf, rd_cf, sv_cf, scen_count
)
    # Caclulate RCI for counterfactual and intervention
    rci_temp_iv = (rc_iv.data + rd_iv.data + sv_iv.data) ./ 3
    rci_temp_cf = (rc_cf.data + rd_cf.data + sv_cf.data) ./ 3
    idx_scens = scen_count:(scen_count + length(rc_iv.locations) - 1)

    # Update delta RCI for selected sites
    display_scenario_df[idx_scens, "RCI uplift mean"] =
        st.mean(rci_temp_iv; dims=3)[1, :, 1] .- st.mean(rci_temp_cf; dims=3)[1, :, 1]
    display_scenario_df[idx_scens, "RCI uplift sd"] =
        st.std(rci_temp_iv; dims=3)[1, :, 1] .- st.std(rci_temp_cf; dims=3)[1, :, 1]

    return nothing
end

# Function to key spatial information for a site
function update_scenario_record!(display_scenario_df, loc_data, k_area, scen_idx)
    display_scenario_df[scen_idx, "Reef"] .= loc_data.Reef
    display_scenario_df[scen_idx, "Geomorphic zone"] .= loc_data.habitat
    display_scenario_df[scen_idx, "Site habitable area m2"] .= k_area
    return nothing
end

################################## Load function for creating data summary table ########################################

function create_biodiversity_metric_summary_df(rs_filepath, aadpt_dhw_range, n_corals_range;
    display_years=[2025, 2030, 2035, 2040, 2045, 2050, 2055, 2060, 2065, 2070], rc_ref=0.70,
    biodiversity_metric_summary_filename=[]
)
    # Load result set
    rs = ADRIA.load_results(rs_filepath)
    # Scenario dataframe
    scenario_df = rs.inputs

    # Calculate metrics for all scenarios
    rc = ADRIA.metrics.relative_cover(rs)
    sv = relative_sv(rs)
    rd = simps_D(rs)

    years = collect(rc.timesteps)
    k_area = ADRIA.site_k_area(rs) # site habitable area in m²
    n_selected = length(rc.locations) # number of locations included (all in domain)
    scen_count = 1

    # Set-up dataframe storage
    n_scens =
        (
            length(display_years)*length(n_corals_range)*length(aadpt_dhw_range)*n_selected
        )+length(display_years)*n_selected

    display_scenario_df = DataFrame("Year"=>zeros(Int64, (n_scens,)),
        "Site"=>fill("", n_scens),
        "Reef"=>fill("", n_scens),
        "Geomorphic zone"=>fill("", n_scens),
        "Intervention"=>zeros(Int64, (n_scens,)),
        "Deployment Volume"=>zeros(Int64, (n_scens,)),
        "Coral Cover"=>zeros(Float64, (n_scens,)),
        "Coral Cover sd"=>zeros(Float64, (n_scens,)),
        "Diversity"=>zeros(Float64, (n_scens,)),
        "Diversity sd"=>zeros(Float64, (n_scens,)),
        "Shelter Volume"=>zeros(Float64, (n_scens,)),
        "Shelter Volume sd"=>zeros(Float64, (n_scens,)),
        "RCI"=>zeros(Float64, (n_scens,)),
        "RCI sd"=>zeros(Float64, (n_scens,)),
        "RCI uplift mean"=>zeros(Float64, (n_scens,)),
        "RCI uplift sd"=>zeros(Float64, (n_scens,)),
        "Site habitable area m2"=>zeros(Float64, (n_scens,)),
        "Deployment area m2"=>zeros(Float64, (n_scens,)),
        "RCI uplift X deployment area m2"=>zeros(Float64, (n_scens,)),
        "deployment site flag"=>zeros(Float64, (n_scens,)),
        "site lat"=>zeros(Float64, (n_scens,)),
        "site long"=>zeros(Float64, (n_scens,)))

    # Caclualte most frequently selected sites for all intervention scenarios
    iv_scens = findall(scenario_df.guided .== 1)
    freq_rank = ADRIA.decision.selection_ranks(
        rs.ranks[:, :, :, iv_scens], :seed; desc=true
    )

    iv_seed_log = rs.seed_log # Number of corals outplanted for each scenario, timestep, location

    site_geom = rs.loc_data.geom # Site polygon geometry
    site_lat = ag.gety.(GeoDataFrames.centroid.(site_geom), 0) # Centroid latitude
    site_long = ag.getx.(GeoDataFrames.centroid.(site_geom), 0) # Centroid longitude

    for yr in display_years # For each year included in the table
        # Find indices for yr, cf scenarios
        yr_scens = findall(years .== yr)
        scen_ids = findall(scenario_df.guided .== -1)
        scen_idx = scen_count:(scen_count + n_selected - 1)

        display_scenario_df[scen_idx, "Year"] .= yr
        display_scenario_df[scen_idx, "Intervention"] .= 0
        display_scenario_df[scen_idx, "Deployment Volume"] .= 0
        display_scenario_df[scen_idx, "Deployment area m2"] .= 0
        display_scenario_df[scen_idx, "site lat"] .= site_lat
        display_scenario_df[scen_idx, "site long"] .= site_long

        # Select metrics for these scenarios
        rc_cf = rc[yr_scens, :, scen_ids] ./ rc_ref # Reference level
        rd_cf = rd[yr_scens, :, scen_ids]
        sv_cf = sv[yr_scens, :, scen_ids]

        # Add metric summaries and spatial info to the dataframe
        update_scenario_record!(display_scenario_df, rs.loc_data, k_area, scen_idx)
        update_scenario_record!(display_scenario_df, rc_cf, rd_cf, sv_cf, scen_count)

        for nc in n_corals_range # For each coral deployment volume
            coral_scens = (rs.inputs.N_seed_TA .== nc/3)

            for aa in aadpt_dhw_range # For each DHW adaptation level
                scen_count+=n_selected
                scen_idx = scen_count:(scen_count + n_selected - 1)
                display_scenario_df[scen_idx, "Year"] .= yr
                display_scenario_df[scen_idx, "Intervention"] .= aa
                display_scenario_df[scen_idx, "Deployment Volume"] .= nc
                display_scenario_df[scen_idx, "site lat"] .= site_lat
                display_scenario_df[scen_idx, "site long"] .= site_long
                display_scenario_df[scen_idx, "Deployment area m2"] .= 0.0

                # Add spatial info to dataframe
                update_scenario_record!(display_scenario_df, rs.loc_data, k_area, scen_idx)

                # Find scenario indices for the given intervention scenario
                scen_ids = findall(
                    coral_scens .& (scenario_df.a_adapt .== aa) .&
                    (scenario_df.guided .== 1)
                )

                # Find the most frequently selected sites for this intervention scenario
                freq_rank = ADRIA.decision.selection_ranks(
                    rs.ranks[:, :, :, scen_ids], :seed; desc=true
                )
                # Get intervention years which are less than or equal to the current year
                iv_years = collect(
                    years[Int64(scenario_df.seed_year_start[1])]:years[Int64(
                        scenario_df.seed_year_start[1] + scenario_df.seed_years[1]
                    )]
                )
                yr_scens_log = findall(in.(years, Ref(iv_years[iv_years .<= yr])))

                # Get number of corals outplanted in this scenario and sum over years up to current one
                iv_dep_nums = dropdims(
                    ADRIA.mean(
                        dropdims(
                            sum(iv_seed_log[yr_scens_log, :, :, scen_ids]; dims=(1, 2));
                            dims=(1, 2)
                        );
                        dims=2
                    );
                    dims=2
                )

                # Estimate deployment area as (number of corals outplanted)/(5 corals/m²)
                display_scenario_df[scen_idx, "Deployment area m2"] .= iv_dep_nums ./ 5

                # Get metrics for this intervention scenario
                rc_iv = rc[yr_scens, :, scen_ids] ./ rc_ref # reference level
                rd_iv = rd[yr_scens, :, scen_ids]
                sv_iv = sv[yr_scens, :, scen_ids]

                # Add metric summaries to table
                update_scenario_record!(
                    display_scenario_df, rc_iv, rd_iv, sv_iv, scen_count
                )
                # Add difference to cf summary to table
                update_scenario_record!(
                    display_scenario_df,
                    rc_iv,
                    rd_iv,
                    sv_iv,
                    rc_cf,
                    rd_cf,
                    sv_cf,
                    scen_count
                )
            end
        end
        scen_count+=n_selected
    end

    # Add RCI uplift * estimated deployment area to table
    display_scenario_df[:, "RCI uplift X deployment area m2"] .=
        display_scenario_df[:, "Deployment area m2"] .*
        display_scenario_df[:, "RCI uplift mean"]

    # Add deployment site flag (1 if a deployment site, 0 if not)
    display_scenario_df[:, "deployment site flag"] .=
        (display_scenario_df[:, "Deployment area m2"] .> 0.0) .* 1.0

    # If no filename provided, use the default
    if isempty(biodiversity_metric_summary_filename)
        aadapt_str = string(["$(aadapt)_" for aadapt in aadpt_dhw_range]...)
        n_corals_str = string(["$(n_c)_" for n_c in n_corals_range]...)
        biodiversity_metric_summary_filename = string(
            "biodiversity_metric_summary_n_corals_", n_corals_str, "a_adapt_",
            aadapt_str, rs.name, ".csv")
    end

    # Write to csv
    CSV.write(biodiversity_metric_summary_filename, display_scenario_df)

    # Return dataframe and filename used to save as csv
    return display_scenario_df, biodiversity_metric_summary_filename
end

###################################### Key filename from ADRIA.jl dataset ###############################################

# Results filepath
rs_filepath = "path to ADRIA resultset"

####################################### Key scenarios to include in dataframe ###########################################

aadpt_dhw_range = [0, 1, 5, 10, 15, 20] # Adaptation levels to include
n_corals_range = [200000, 500000, 1000000] # Deployment volumes to include

################################## Create biodiversity metric summary dataframe #########################################

display_scenario_df, bio_metric_sum_fn = create_biodiversity_metric_summary_df(
    rs_filepath, aadpt_dhw_range, n_corals_range
)
