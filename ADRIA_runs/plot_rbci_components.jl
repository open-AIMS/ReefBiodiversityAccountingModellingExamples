############################################ Description ################################################################
# Script for ploting components of the RBCI from a set of ADRIA runs including taxa disaggregated coral cover vs. simpsons
# diversity, shelter volume and coral cover.

############################################ Load packages ##############################################################
using ADRIA
using CSV
using DataFrames, YAXArrays
using GLMakie, GeoMakie, GraphMakie
using DimensionalData
using Statistics

#################################### Set filepaths for ADRIA domain and resultset #######################################
rs_filepath = "C:\\Users\\rcrocker\\Documents\\Github\\ADRIA.jl\\sandbox\\march2025-milestone-runs\\Outputs\\Moore_2025-03-18_v070_rc1__RCPs_45__2025-05-08_14_12_50_575"
rs = ADRIA.load_results(rs_filepath)

######################################## Load temporal plotting function ################################################

# Function to plot data over time with quantile spread
function temporal_spread!(
    ax::Axis, data::YAXArray, years::Vector{Int64};
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)
    n_timesteps = length(years)
    data_quants = zeros(n_timesteps, 3)

    for tt ∈ 1:n_timesteps
        data_quants[tt, :] = quantile(
            data[tt, :, :], [0.35, 0.5, 0.65]
        )
    end

    plot_color = get(opts, :plot_color, :dodgerblue)

    band!(
        ax,
        years,
        vec(data_quants[:, 1]),
        vec(data_quants[:, 3]);
        color=plot_color,
        alpha=0.2
    )

    lines!(ax,
        years,
        vec(data_quants[:, 2]),
        color=plot_color,
        linewidth=3
    )
    return ax
end

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

################################ Load function for plotting taxa cover vs. diversity ####################################
# Function for plotting diversity vs. taxa diaggregated coral cover over time for a specified intervention scenario
function plot_taxa_vs_diversity(rs, a_adapt, n_corals; n_sites=10, color_vec = Makie.wong_colors())

    # Extract intervention scenarios
    scenario_df = rs.inputs
    iv_scens = findall((scenario_df.guided.==1).&(scenario_df.a_adapt.==a_adapt).&(scenario_df.N_seed_CA.== n_corals/3))

    # Get most frequently selected n_sites
    freq_rank = ADRIA.decision.selection_ranks(rs.ranks[:, :, :, iv_scens], :seed; desc=true)
    iv_sites = freq_rank[1:n_sites]

    # Calculate species cover and Simpson's diversity
    sp_cov = species_cover(rs)
    rd = simps_D(rs)

    # Years to plot
    years = collect(2025:2025+length(rd.timesteps)-1)

    # Taxa names
    species_names = ["Tabular Acropora", "Corymbose Acropora", "Corymbose non-Acropora", "Small massives", "Large massives"]

    # Plot cover for each taxa for intervention scenarios
    f1 = Figure(;resolution=(1000, 600))
    ax1 = Axis(
        f1[1, 1];
        title="Relative species cover",
        ylabel="Year",
        xlabel="Relative cover"
    )
    for sp in 1:5
        opts = Dict(:plot_color=>color_vec[sp])
        temporal_spread!(ax1, sp_cov[:, sp, iv_sites, iv_scens], years; opts=opts)
    end

    leg_els = [LineElement(color=cc, strokewidth=3) for cc in color_vec[1:5]]
    Legend(f1[1,2], leg_els, species_names, "Scenario")

    # Plot Simpson's diversity for intervention scenarios
    f2 = Figure(;resolution=(1000, 600))
    ax2 = Axis(
        f2[1, 1];
        title="Diversity",
        ylabel="Year",
        xlabel="Diversity"
    )
    opts = Dict(:plot_color=>:black)
    temporal_spread!(ax2, rd[:, iv_sites, iv_scens], years; opts=opts)

    return f1, f2

end

###################################### Load function for components of RBCI #############################################
# Function for plotting diversity vs. taxa disaggregated coral cover over time for a specified intervention scenario
function plot_metric_over_time(rs, aadpt_dhw_range, n_corals, n_sites; n_years = 35, color_vec = Makie.wong_colors())
    # Extract intervention scenarios
    scenario_df = rs.inputs
    cf_scens = findall(scenario_df.guided.==-1)
    iv_scens = findall(scenario_df.guided.==1)
    # Get most frequently selected n_sites
    freq_rank = ADRIA.decision.selection_ranks(rs.ranks[:, :, :, iv_scens], :seed; desc=true)
    iv_sites = freq_rank

    # Calculate metrics
    rc = ADRIA.relative_cover(rs)
    sd = simps_D(rs)
    rel_sv = relative_sv(rs)
    a_sv = ADRIA.metrics.absolute_shelter_volume(rs)

    # Years to plot
    years = collect(2025:2025+length(rc.timesteps)-1)

    # Legends and labels for plotting
    leg_els = [LineElement(color=cc, strokewidth=3) for cc in color_vec[1:length(aadpt_dhw_range)+1]] # Legend
    dhws_string = ["$dhw DHW" for dhw in aadpt_dhw_range]
    dhws_string = [dhws_string..., "CF"] # Adaptation scenarios
    cf_ids = []
    metrics_vec = [rc, sd, rel_sv, a_sv]
    metric_names = ["Relative cover", "Simpsons Diversity", "Relative Shelter Volume", "Absolute shelter volume dm³/m²"]
    metric_filenames = ["rel_c", "simps_d", "rel_sv", "abs_sv"]

    for (met_idx, met) in enumerate(metrics_vec) # Loop over metrics
        for (n_cor_idx, n_cor) in enumerate(n_corals) # Loop over deployment volume scenarios

            # Set up figure for each metric and deployment volume
            f = Figure(;resolution=(1000, 600))
            ax = Axis(
                f[1, 1];
                title="$(n_cor) corals per year",
                ylabel=metric_names[met_idx],
                xlabel="Year"
            )

            for (aa_idx, aadpt) in enumerate(aadpt_dhw_range)
                # Plot each adaptation scenario
                opts = Dict(:plot_color=>color_vec[aa_idx])
                iv_scens = findall((scenario_df.a_adapt.==aadpt).&(scenario_df.N_seed_CA.==n_cor/3))
                temporal_spread!(ax, met[1:n_years, iv_sites[1:n_sites[n_cor_idx]], iv_scens], years[1:n_years]; opts=opts)

            end

            # Plot same number of counterfactual scenarios
            iv_scens = findall((scenario_df.N_seed_CA.==n_cor/3))
            cf_ids = rand(cf_scens, length(iv_scens))

            opts = Dict(:plot_color=>color_vec[7])
            temporal_spread!(ax, met[1:n_years, iv_sites[1:n_sites[n_cor_idx]], cf_ids], years[1:n_years]; opts=opts)
            Legend(f[1,2], leg_els, dhws_string, "Scenario")

            # Save each figure
            save("rbci_components_$(metric_filenames[met_idx])_$(n_cor)_corals.png", f)
        end
    end
end

################################################# Plot figures ##########################################################
# Plot taxa disaggregated cover vs. Simpson's diversity for a set of intervention scenarios
a_adapt = 10.0
n_corals = 1000000
taxa_cov_plot, div_plot = plot_taxa_vs_diversity(rs, a_adapt, n_corals; n_sites = 20)

# Plot each metric for each intervention scenario over time
aadpt_dhw_range = [0, 1, 5, 10, 15, 20]
n_corals = [200000, 500000, 1000000]
n_sites = [5, 10, 20]
plot_metric_over_time(rs, aadpt_dhw_range, n_corals, n_sites)
