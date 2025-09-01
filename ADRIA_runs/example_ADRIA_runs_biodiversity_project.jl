############################################ Description ################################################################
# Script for running intervention scenarios ADRIA.jl to create a biodiversity accounting table of a series of
# intervention scenarios.

############################################ Load packages ##############################################################
using ADRIA
using DataFrames, YAXArrays

############################################ Load ADRIA Domain ##########################################################
path_to_domain = "C:\\Users\\rcrocker\\Documents\\Datapackages\\Moore_2025-03-18_v070_rc1"
# Load ADRIA domain (e.g. Moore)
dom = ADRIA.load_domain(
    path_to_domain,
    "45"
)

######################################### Define scenario parameters ####################################################
# Set adaptation levels to run
aadpt_dhw_range = [5, 10, 15]
# Set coral deployment levels to run
n_corals_range = [200000, 500000, 1000000]

# Extract parameters related to the coral model to hold these constant
coral_model_factors = ADRIA.model_spec(dom).fieldname[ADRIA.model_spec(dom).component.=="Coral"]

######################################### Create scenario dataframe #####################################################

# Sample 2 scenarios only to get scenario datframe structure
scens = ADRIA.sample(dom, 2^1)

# Loop through intervention scenario settings and sample
for aapt in aadpt_dhw_range
    for dep_eff in n_corals_range
        ADRIA.set_factor_bounds!(dom, :guided, ("counterfactual", "COCOSO",)) # Only run counterfactuals, unguided and guided with 1st algorithm
        ADRIA.fix_factor!(dom, coral_model_factors) # Fix coral model factors
        ADRIA.fix_factor!(dom;
            a_adapt = aapt, # Set adaptation DHW level
            N_seed_TA = dep_eff/3, # Set coral deployment level (split evenly across species)
            N_seed_CA = dep_eff/3,
            N_seed_SM = dep_eff/3,
            seed_year_start = 1, # Start seeding in 2025
            shade_year_start = 0,
            fog_year_start = 0,
            seed_deployment_freq = 1, # Choose sites to deploy at in year one and deploy at the same sites
            fog_deployment_freq = 0,
            shade_deployment_freq = 0,
            cyclone_mortality_scenario = 0,
            wave_scenario = 0,
            plan_horizon = 10, # Look 10 years into the future (for DHWs) to select sites
            seed_years = 10, # Seed for a total of 10 years
            fogging = 0.0, # No fogging or SRM
            SRM = 0.0,
            fog_years = 0.0,
            shade_years = 0.0,
            min_iv_locations = 15.0, # Choose a maximum of 15 locations to seed at
            seed_heat_stress = 1.0, # Use heat stress to choose sites
            seed_wave_stress = 0.0,
            seed_in_connectivity = 0.0,
            seed_out_connectivity = 0.0,
            seed_depth = 1.0, # Use depth to choose sites
            seed_coral_cover = 1.0) # Use available space to choose sites

        # Sample to get scenario dataframe
        scens_temp = ADRIA.sample(dom, 2^8)

        # Append scenario dataframe to others so all scenarios can be run at once
        global scens = vcat(scens, scens_temp)

    end
end
# Remove initial 2 scenarios which were sampled without settings
scens = scens[3:end,:]

######################################### Run scenarios #####################################################

rs = ADRIA.run_scenarios(dom, scens, ["45"])
