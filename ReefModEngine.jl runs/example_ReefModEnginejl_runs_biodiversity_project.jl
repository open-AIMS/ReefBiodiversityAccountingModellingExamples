############################################ Description ################################################################
# Script for running ReefModEngine.jl and saving cell level data to create a biodiversity accounting table of a series of
# intervention scenarios.

############################################ Load packages ##############################################################
using ReefModEngine
using CSV, DataFrames
using NCDatasets

# Initialize RME (may take a minute or two)
init_rme("path to rme")

############################################ RME options ################################################################
set_option("thread_count", 2)  # Set to use two threads
set_option("use_fixed_seed", 1)  # Turn on use of a fixed seed value
set_option("fixed_seed", 123.0)  # Set the fixed seed value
set_option("cots_enabled", 0)  # Turn COTS off
set_option("cyclones_enabled", 0)  # Turn cyclones off
set_option("initial_set_fixed", 1)  # Fix initial cover, so not randomised at each rep
set_option("recovery_value_enabled", 0)  # Don't save recovery value to save runtime

############################################ Load reef ids ##############################################################
# Load target reefs (to run and intervene at)
# The first column is simply the row number.
# The second column is a list of target GBRMPA reef ids

# These are the reef ids we want to run (if not running the full GBR)
loc_details = CSV.read(
    "ext_moore_reefset.csv",
    DataFrame;
    header=["index_id", "reef_id"],
    types=Dict(1=>Int64, 2=>String)  # Force values to be interpreted as expected types
)
# These are the reef ids we want to intervene at
iv_loc_details = CSV.read(
    "iv_ext_moore_reefset.csv",
    DataFrame;
    header=["index_id", "reef_id"],
    types=Dict(1=>Int64, 2=>String)  # Force values to be interpreted as expected types
)

# Get fully list of reef ids as specified by ReefMod Engine
reef_id_list = reef_ids()

# Reef indices and IDs for reefset to be run
target_reef_ids = loc_details.reef_id
# Find indices of these reefs in the ReefMod list
target_reef_idx = findall(in.(reef_id_list, Ref(target_reef_ids)))
n_target_reefs = length(target_reef_idx)

# Reef indices and IDs for intervention reef
target_reef_ids_moore = iv_loc_details.reef_id
moore_idx = findall(target_reef_ids_moore .== reef_id_list)
# Find indices of this(ese) reef(s) the chosen list of reefs to run
ext_moore_idx = findall(in.(reef_id_list, Ref(target_reef_ids)))
# Specify as vector (if only one reef it still needs to be a Vector)
target_reef_ids_moore = Vector(target_reef_ids_moore)
n_iv_target_reefs = length(moore_idx)

############################################ Setup run parameters ######################################################
runs_name = "Milestone runs March 2025" # Name to associate with this set of runs
start_year = 2025
end_year = 2060
years = collect(start_year:end_year)
n_years = (end_year-start_year) + 1
RCP_scen = "SSP 2.45"  # RCP/SSP scenario to use
reps = 5 # Number of repeats: number of randomised environmental sequences to run
n_reefs = length(reef_id_list)

# Get full set of reef areas from RME
reef_area_km² = reef_areas()

# Get list of areas for the target reefs
iv_reef_areas_km² = reef_areas(target_reef_ids_moore)

# Define coral outplanting density (per m² per species)
d_density_m² = 5.0/6

# Intervention scenarios
aadpt_dhw_range = [5.0, 10.0, 20.0]
n_corals_range = [1000000]
# Std deviation for heat tolerance (needs to be set if heat tolerance mean is specified)
ht_tol_std = [2.49, 2.49, 2.49, 2.49, 2.49, 2.49]

# Climate models - in this case run just the first 3
gcm_names = [@RME gcmName(i::Cint)::Cstring for i in 1:3]

# Intervention years
iv_years = 2025:2030

############################################ Setup storage #############################################################
# Initialize result store
result_store = ResultStore(start_year, end_year)
# Set up storage for additional data logging
n_scens = length(aadpt_dhw_range)*length(n_corals_range)*length(gcm_names)
# Dataframe to save GCM, adaptation level and number of corals outplanted for each scenario
scenario_storage = DataFrame(
    hcat(fill([""], n_scens), zeros(n_scens, 2)), ["GCM", "n_corals", "a_adapt"]
)

# Storage for cell coral cover data and cell restored indices
cell_indices_store = zeros(
    Int64, length(aadpt_dhw_range)*length(n_corals_range)*length(gcm_names), reps, 200
)
cell_storage = zeros((
    length(aadpt_dhw_range)*length(n_corals_range)*length(gcm_names),
    (length(years)-1),
    reps,
    n_target_reefs,
    7,
    2
))

# Storage for shelter volume for each species
sv_storage = zeros((
    length(aadpt_dhw_range)*length(n_corals_range)*length(gcm_names),
    (length(years)-1),
    reps,
    n_target_reefs,
    6,
    2
))
species_sv = zeros(length(ext_moore_idx))
n_iv_years = length(iv_years)

# Storage for number of cells tracked in each run (may change but usually 100)
n_cells = zeros(Int64, length(aadpt_dhw_range)*length(n_corals_range)*length(gcm_names))

######################################### Loop through scenarios #######################################################

# Counter to track scenario number
scen_count_ind = 0

for a_adapt in aadpt_dhw_range
    for n_corals in n_corals_range
        for gcm in gcm_names
            @info "Starting runs"
            reset_rme()  # Reset RME to clear any previous runs

            # Create run
            @RME runCreate(
                runs_name::Cstring,
                start_year::Cint,
                end_year::Cint,
                RCP_scen::Cstring,
                gcm::Cstring,
                reps::Cint
            )::Cint

            # Create subset of reefs to run and subset to intervene on
            @RME reefSetAddFromIdList(
                "ext_moore_set"::Cstring,
                target_reef_ids::Ptr{Cstring},
                length(target_reef_ids)::Cint
            )::Cint
            @RME reefSetAddFromIdList(
                "iv_moore_set"::Cstring,
                target_reef_ids_moore::Ptr{Cstring},
                length(target_reef_ids_moore)::Cint
            )::Cint

            # Define outplanting intervention
            set_outplant_deployment!(
                "outplant_moore_iv",
                "iv_moore_set",
                n_corals,
                Int.((n_iv_years-1)*n_corals),
                iv_years[1],
                iv_years[end],
                1,
                iv_reef_areas_km²,
                fill(d_density_m², 6)
            )
            # Set DHW tolerance for intervention
            @RME ivSetOutplantHeatToleranceMeanDhw(
                "outplant_moore_iv"::Cstring, fill(a_adapt, 6)::Ptr{Cdouble}, 6::Cint
            )::Cint
            @RME ivSetOutplantHeatToleranceSdDhw(
                "outplant_moore_iv"::Cstring, ht_tol_std::Ptr{Cdouble}, 6::Cint
            )::Cint

            # Initialise run
            run_init()
            # Set reefset to run
            @RME runSetReefSet("ext_moore_set"::Cstring)::Cint

            global scen_count_ind+=1

            # Storing cell count, GCM, adaptation level and number of corals outplanted
            n_cells[scen_count_ind] = @getRME runCellCount()::Cint
            scenario_storage[scen_count_ind, "GCM"] = gcm
            scenario_storage[scen_count_ind, "n_corals"] = n_corals
            scenario_storage[scen_count_ind, "a_adapt"] = a_adapt

            for (yr_idx, yr) in enumerate(years[1:(end - 1)])
                # Process each year
                @RME runProcessYears(1::Cint)::Cint

                for rep in 1:reps
                    # Get restored cell indices
                    if yr .== iv_years[1]
                        cell_indices = zeros(
                            Float64,
                            @getRME runRestoredCellCount(
                                moore_idx[1]::Cint, rep::Cint
                            )::Cint
                        )
                        @RME runRestoredCellIndices(
                            moore_idx[1]::Cint,
                            rep::Cint,
                            cell_indices::Ptr{Cdouble},
                            length(cell_indices)::Cint
                        )::Cint
                        cell_indices_store[scen_count_ind, rep, 1:length(cell_indices)] .=
                            cell_indices
                    end

                    for (reef_idx, reef_id) in enumerate(ext_moore_idx)
                        cell_data = zeros(Float64, n_cells[scen_count_ind])
                        cell_indices = cell_indices_store[scen_count_ind, rep, :]
                        cell_indices = cell_indices[cell_indices .> 0.0]

                        # For counterfactual
                        @RME runGetCellData(
                            "coral_cm2"::Cstring,
                            reef_id::Cint,
                            0::Cint,
                            rep::Cint,
                            cell_data::Ptr{Cint},
                            length(cell_data)::Cint
                        )::Cint
                        cell_storage[scen_count_ind, yr_idx, rep, reef_idx, 1, 1] = sum(
                            cell_data[Int.(cell_indices)]
                        )
                        # For intervention
                        @RME runGetCellData(
                            "coral_cm2"::Cstring,
                            reef_id::Cint,
                            1::Cint,
                            rep::Cint,
                            cell_data::Ptr{Cint},
                            length(cell_data)::Cint
                        )::Cint
                        cell_storage[scen_count_ind, yr_idx, rep, reef_idx, 1, 2] = sum(
                            cell_data[Int.(cell_indices)]
                        )

                        for sp in 2:7
                            # Cover in each cell for each species
                            species_code = "species_$(sp-1)_cm2"
                            @RME runGetCellData(
                                species_code::Cstring,
                                reef_id::Cint,
                                0::Cint,
                                rep::Cint,
                                cell_data::Ptr{Cint},
                                length(cell_data)::Cint
                            )::Cint
                            cell_storage[scen_count_ind, yr_idx, rep, reef_idx, sp, 1] = sum(
                                cell_data[Int.(cell_indices)]
                            )
                            @RME runGetCellData(
                                species_code::Cstring,
                                reef_id::Cint,
                                1::Cint,
                                rep::Cint,
                                cell_data::Ptr{Cint},
                                length(cell_data)::Cint
                            )::Cint
                            cell_storage[scen_count_ind, yr_idx, rep, reef_idx, sp, 2] = sum(
                                cell_data[Int.(cell_indices)]
                            )
                        end
                    end

                    for sp in 1:6
                        # Shelter volume in each cell for each species
                        species_code = "sp$(sp)_shelter_volume_dm3_per_m2"
                        @RME runGetData(
                            species_code::Cstring,
                            "ext_moore_set"::Cstring,
                            0::Cint,
                            yr::Cint,
                            rep::Cint,
                            species_sv::Ptr{Cdouble},
                            length(species_sv)::Cint
                        )::Cint
                        sv_storage[scen_count_ind, yr_idx, rep, :, sp, 1] .= species_sv
                        @RME runGetData(
                            species_code::Cstring,
                            "ext_moore_set"::Cstring,
                            1::Cint,
                            yr::Cint,
                            rep::Cint,
                            species_sv::Ptr{Cdouble},
                            length(species_sv)::Cint
                        )::Cint
                        sv_storage[scen_count_ind, yr_idx, rep, :, sp, 2] .= species_sv
                    end
                end
            end

            # Collect and store results
            concat_results!(result_store, start_year, end_year, reps)
        end
    end
end

############################################# Save resultset ###########################################################

# Save regular resultset
save_result_store("ext_moore_results_milestone2025", result_store)

# Save intervention scenario parameters dataframe
CSV.write("scenario_par_log_2025milestoneruns.csv", scenario_storage)

# Save cell-level data as a netcdf
ds = NCDataset("cell_log_2025marchmilestoneruns.nc", "c")
defDim(ds, "cells", 200)
defDim(ds, "int_cf", 2)
defDim(ds, "scenarios", n_scens)
defDim(ds, "species", 7)
defDim(ds, "reef_ids", 19)
defDim(ds, "taxa", 6)
defDim(ds, "reps", reps)
defDim(ds, "iv_years", Int(length(iv_years)))
defDim(ds, "years", Int(length(years)-1))

# Save reef ids
reef_ids_var = defVar(ds, "reef_ids", String, ("reef_ids",))
reef_ids_var = target_reef_ids

# Save reef areas
reef_ids_var = defVar(ds, "reef_areas", Float32, ("reef_ids",))
reef_ids_var = reef_areas(target_reef_ids)

# Save species level shelter volume
species_sv_var = defVar(
    ds,
    "species_sv_var",
    Float32,
    ("scenarios", "years", "reps", "reef_ids", "taxa", "int_cf")
)
species_sv_var[:, :, :, :, :, :] = sv_storage

# Save cell-level species-level coral cover in cm2
cell_data_var = defVar(
    ds,
    "reef_cell_cover_cm",
    Float32,
    ("scenarios", "years", "reps", "reef_ids", "species", "int_cf")
)
cell_data_var[:, :, :, :, :, :] = cell_storage

# Save number of cells used for each reef in each scenario
cell_counts = defVar(ds, "cell_counts", Int32, ("scenarios",))
cell_counts[:] = n_cells

# Save number of cells intervened on in each scenario
cell_indices_var = defVar(ds, "cell_indices_counts", Int32, ("scenarios", "reps"))
cell_indices_var[:, :] = dropdims(sum(cell_indices_store .> 0; dims=3); dims=3)
close(ds)
