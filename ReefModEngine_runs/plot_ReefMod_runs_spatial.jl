import GeoDataFrames as GDF
using CSV
using DataFrames, YAXArrays
using GLMakie, GeoMakie, GraphMakie
using DimensionalData
using Statistics

function spatial_map(
    geo_df::DataFrame,
    color_vec::Vector{String};
    fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()),
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    f = Figure(; fig_opts...)
    return spatial_map(f, geo_df, color_vec; axis_opts=axis_opts, opts=opts)
end
function spatial_map(
    f::Figure,
    geo_df::DataFrame,
    color_vec::Vector{String};
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    legend_name = get(opts, :legend_name, "Class")
    color_map = get(opts, :color_map, :batlow25)

    levels = unique(color_vec)

    colors = Makie.categorical_colors(color_map, length(levels))
    legend_elems::Vector{Any} = ones(length(levels))
    spatial = GeoAxis(
        f[1, 1];
        axis_opts...
    )

    for (int_lev, lev) in enumerate(levels)
        temp_geo = geo_df[findall(color_vec .== lev), :]
        geo_data = GeoMakie.to_multipoly(temp_geo[:, :geom])

        poly!(
            spatial,
            geo_data;
            color=colors[int_lev],
            label=lev,
            transparency=true
        )
    end
    legend_elems = [PolyElement(; color=colors[int_lev]) for int_lev in 1:length(levels)]
    Legend(f[1, 2], legend_elems, levels, legend_name)

    return f
end
function spatial_map(
    geo_df::DataFrame,
    color_vec::Union{Vector{Float32},Vector{Float64}};
    fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()),
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    f = Figure(; fig_opts...)
    return spatial_map(f, geo_df, color_vec; axis_opts=axis_opts, opts=opts)
end
function spatial_map(
    f::Figure,
    geo_df::DataFrame,
    color_vec::Union{Vector{Float32},Vector{Float64}};
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    colorbar_label = get(opts, :colorbar_label, "Metric")
    color_map = get(opts, :color_map, :batlow25)
    color_range = get(opts, :color_range, (0.0, maximum(color_vec)))

    geo_data = GeoMakie.to_multipoly(geo_df[:, :geom])
    spatial = GeoAxis(
        f[1, 1];
        axis_opts...
    )
    poly!(
        spatial,
        geo_data;
        color=color_vec,
        colormap=color_map,
        colorrange=color_range
    )
    Colorbar(
        f[1, 2];
        colorrange=color_range,
        colormap=color_map,
        label=colorbar_label,
        height=Relative(0.70)
    )
    return f
end

"""
    set_figure_defaults(fig_opts::Dict{Any,Any})::Dict{Symbol,Any}

Set default figure settings for spatial figures
"""
function set_figure_defaults(
    fig_opts::Dict{Any,Any}
)::Dict{Symbol,Any}
    fig_opts[:size] = get(fig_opts, :size, (600, 700))
    return fig_opts
end

"""
    set_axis_defaults(axis_opts::Dict{Any,Any})::Dict{Symbol,Any}

Set default axis settings for spatial figures
"""
function set_axis_defaults(
    axis_opts::Dict{Any,Any}
)::Dict{Symbol,Any}
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")
    axis_opts[:dest] = get(axis_opts, :dest, "+proj=latlong +datum=WGS84")
    axis_opts[:xgridvisible] = get(axis_opts, :xgridvisible, false)
    axis_opts[:ygridvisible] = get(axis_opts, :ygridvisible, false)
    axis_opts[:yticklabelsvisible] = get(axis_opts, :yticklabelsvisible, false)
    axis_opts[:xticklabelsvisible] = get(axis_opts, :xticklabelsvisible, false)
    return axis_opts
end

"""
    plot_RBCI_spatial(RBCI_df_fn::String, gbr_spatial_fn::String, target_dhw::Float64, target_volume::Int64; RBCI_colorrange=(0.0,1.0), RBCI_uplift_colorrange=(0.0, 50.0))

Spatially plot the intervention and counterfactual RBCI and uplift
RBCI in km2 for multiple years from a biodiversity accounting summary dataframe.
"""
function plot_RBCI_spatial(RBCI_df_fn::String, gbr_spatial_fn::String,
    target_dhw::Float64, target_volume::Int64; RBCI_colorrange=(0.0,1.0),
    RBCI_uplift_colorrange=(0.0, 50.0))

    delta_RBCI_df = CSV.read(RBCI_df_fn, DataFrame)
    gbr_spatial = GDF.read(gbr_spatial_fn)

    reefset_ids = [findall(gbr_spatial.LABEL_ID.==rf)[1] for rf in lowercase.(delta_RBCI_df.Reef)]

    delta_RBCI_df[!, "reef_idx"] = reefset_ids
    delta_RBCI_df[!, "geom"] = gbr_spatial.geom[reefset_ids]
    delta_RBCI_df[!, "RBCI uplift area hectares (multiplied by reef k area X intervention area proportion)"] = delta_RBCI_df[:,"RBCI uplift area m2 x estimated intervention area"].*0.0001
    yr_slices = unique(delta_RBCI_df.Year)

    for yr_slice in yr_slices
        cf_scens = delta_RBCI_df[(delta_RBCI_df.Year.==yr_slice).&(delta_RBCI_df.Intervention.==0), :]
        iv_scens = delta_RBCI_df[(delta_RBCI_df.Year.==yr_slice).&(delta_RBCI_df.Intervention.==target_dhw).&(delta_RBCI_df[:,"Deployment Volume"].==target_volume), :]

        base_filename = "rme_spatial_plots_$(yr_slice)_dhw$(target_dhw)_ncorals$(target_volume)"

        axis_opts = set_axis_defaults(Dict())
        axis_opts[:title] = "$(yr_slice) counterfactual"
        opts = Dict(:color_range=>RBCI_colorrange, :color_map=>:roma, :colorbar_label=>"Mean RBCI")
        f1 = spatial_map(cf_scens, cf_scens[:,"RBCI"]; opts=opts, axis_opts=axis_opts)
        save(string(base_filename,"_counterfactual_RBCI", ".png"),f1)

        axis_opts[:title] = "$(yr_slice) $(target_dhw)DHW enhancement, $(target_volume) corals per year for 10 years"
        f2 = spatial_map(iv_scens, iv_scens[:,"RBCI"]; opts=opts, axis_opts=axis_opts)
        save(string(base_filename,"_intervention_RBCI", ".png"),f2)

        axis_opts[:title] = "$(yr_slice) uplift area for intervention vs. counterfactual"
        opts = Dict(:color_map=>:roma, :color_range=>RBCI_uplift_colorrange, :colorbar_label=>"RBCI uplift x k area hectares")
        f3 = spatial_map(iv_scens, iv_scens[:, "RBCI uplift area hectares (multiplied by reef k area X intervention area proportion)"]; opts=opts, axis_opts=axis_opts)
        save(string(base_filename,"_uplift_RBCI", ".png"),f3)

    end
end

# Example usage - plot RBCI spatial plots from an RME biodiversity
# accounting summary dataframe (created using 'create_metric_summary_dataframe_ReefModEngine.jl')
RBCI_df_fn = "Path to biodiversity accounting summary dataframe"
gbr_spatial_fn = "Path to GBR spatial file in RME datapackage (e.g. reefmod_gbr.gpkg)"

target_dhw = 20.0 # DHW enhancement level to plot
target_volume = 1000000 # Coral outplant volume to plot
RBCI_colorrange=(0.0, 0.55) # RBCI range to plot
RBCI_uplift_colorrange = (0.0, 60.0) # RBCI uplift range (in km2) to plot

# Create plots for each year in accounting table
plot_RBCI_spatial(RBCI_df_fn, gbr_spatial_fn,
    target_dhw, target_volume; RBCI_colorrange=RBCI_colorrange,
    RBCI_uplift_colorrange=RBCI_uplift_colorrange)
