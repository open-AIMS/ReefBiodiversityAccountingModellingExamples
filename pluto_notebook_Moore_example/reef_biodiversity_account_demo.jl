### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8d27ffec-e4d5-46bc-8218-b5a44adbb796
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(".")
	using DataFrames
	using GLMakie, GeoMakie, GraphMakie
	using GeoDataFrames
	using Distances
	using DimensionalData
	using Statistics
	import ReefBiodiversityAccountSetup as RBAS
	using PlutoUI
	import PlutoUI:combine
	
end

# ╔═╡ a2529343-55a8-43ea-b066-f928d52ee534
md"# Visualise project area"

# ╔═╡ 570139e5-cd28-4689-a2c8-0a3682f5e572
md"## Benthic and geomorphic zones"

# ╔═╡ 432ec934-ebaa-430b-9401-e64369a9ecc7
begin
	function spatial_map(
	    geo_df::DataFrame,
	    color_vec::Union{Vector{Float32},Vector{Float64}};
	    fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()),
	    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
	    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), 
		geom_col::Symbol=:geometry
	)::Figure
	    f = Figure(; fig_opts...)
	    return spatial_map(f, geo_df, color_vec; axis_opts=axis_opts, opts=opts, geom_col=geom_col)
	end
	function spatial_map(
	    f::Figure,
	    geo_df::DataFrame,
	    color_vec::Union{Vector{Float32},Vector{Float64}};
	    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
	    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
		geom_col::Symbol=:geometry
	)::Figure
	    colorbar_label = get(opts, :colorbar_label, "Metric")
	    color_map = get(opts, :color_map, :batlow25)
	    color_range = get(opts, :color_range, (0.0, maximum(color_vec)))
	
	    geo_data = GeoMakie.to_multipoly(geo_df[:, geom_col])
	    spatial = GeoAxis(
	        f[1, 1];
        	aspect=DataAspect(),
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
	function set_figure_defaults(
	    fig_opts::Dict{Any,Any}
	)::Dict{Symbol,Any}
	    fig_opts[:size] = get(fig_opts, :size, (600, 700))
	
	    return fig_opts
	end
	function set_axis_defaults(axis_opts::Dict{Any,Any})::Dict{Symbol,Any}
    	axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    	axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")
    	axis_opts[:dest] = get(axis_opts, :dest, "+proj=latlong +datum=WGS84")
		axis_opts[:xgridvisible] = get(axis_opts, :xgridvisible, false)
    	axis_opts[:ygridvisible] = get(axis_opts, :ygridvisible, false)
		axis_opts[:yticklabelsvisible] = get(axis_opts, :yticklabelsvisible, false)
		axis_opts[:xticklabelsvisible] = get(axis_opts, :xticklabelsvisible, false)
    	return axis_opts
	end

end

# ╔═╡ 4c135067-265b-4e98-bb31-a6daa9d426ad
begin
	config_file = RBAS.spatial_analysis.load_config(; config_path="C:\\Users\\rcrocker\\Documents\\Datapackages\\config.toml")
	benthic, geomorphic, extent_base = RBAS.spatial_analysis.load_spatial_base(config_file);

end

# ╔═╡ 2a5cf51c-8831-4aff-8419-3603d42030e3
md"#### Benthic zones"

# ╔═╡ 646613cb-740c-4ef7-b406-5cdd31d9dcc4
begin
	benthic_classes_fig = spatial_map(
    benthic, benthic[:, :class];
    opts=Dict(:legend_name => "Benthic class", :color_map => :Accent_6)
	)
	#save("benthic_plot.png", benthic_classes_fig)
end

# ╔═╡ 89475ad4-0102-4ea2-bbf0-497da8d935a2
md"#### Geomorphic zones"

# ╔═╡ df182f3d-29d2-4b3f-bd45-391a8a2dfee9
begin
	geo_classes_fig = spatial_map(
    geomorphic, geomorphic[:, :class];
    opts=Dict(:legend_name => "Geomorphic class", :color_map => :hawaii10)
)
	#save("geomorphic_plot.png", geo_classes_fig)
end

# ╔═╡ bfaed850-ecd4-4319-8ba4-1c48670a0a85
md"#### Habitable area (Coral/Algal/Hard rock)"

# ╔═╡ a0adc075-5289-45dc-b23e-b681d2785943
begin
	benthic_filtered = RBAS.spatial_analysis.filter_site_area(benthic)
	fig_geomorphic_ext = spatial_map(
    benthic_filtered,
    benthic_filtered[:, :class];
    opts=Dict(:legend_name => "Benthic class",
        :color_map => Makie.wong_colors())
)
	#save("rock_coral_plot.png", fig_geomorphic_ext)
end

# ╔═╡ 4a9ae439-da76-4266-98e0-849eae994020
md"#### Geomorphic zones of habitable area"

# ╔═╡ 4877ea3e-8bf9-4c08-9694-14864bfb7ad9
begin
	# Get intersection of geomorphic polygons and account extent
	geomorphic_ext = RBAS.spatial_analysis.multipoly_geom_intersection(extent_base, geomorphic, :class)

	# Get intersection of benthic filtered polygons and account extent
	benthic_ext = RBAS.spatial_analysis.multipoly_geom_intersection(extent_base, benthic_filtered, :class)

	# Get intersection of benthic filtered polygons and geomorphic polygons
	geomorphic_benthic_comb = RBAS.spatial_analysis.multipoly_geom_intersection(benthic_ext, geomorphic_ext, :class)

	fig_geomorphic_filtered = spatial_map(
    geomorphic_benthic_comb,
    geomorphic_benthic_comb[:, :class];
    opts=Dict(:legend_name => "Geomorphic class",
        :color_map => :hawaii10)
)
	#save("benthic_filtered_geomorphic_plot.png", fig_geomorphic_filtered)

end

# ╔═╡ 05aa9fcd-12ab-48db-a657-827366087583
md"## Median depth map"

# ╔═╡ 844f575d-1137-45e9-9d9a-fd84c0e71c0d
begin
	# Extract depths from raster file
	geomorphic_depths, depths = RBAS.spatial_analysis.median_features_allen(
    geomorphic, config_file; is_depth=true
)
	fig_depths = spatial_map(geomorphic_depths, geomorphic_depths[:, :depth_med]; opts=Dict(:color_map => Reverse(:viridis), :colorbar_label=>"Median depth"), geom_col=:geom)
	#save("med_depths_plot.png", fig_depths)
end

# ╔═╡ 3aa16739-3c70-4895-82e0-5341da5b9956
md"## Median max dhw"

# ╔═╡ 0e11066d-8afa-48a1-9344-c2d2b7c98223
begin
	# Extract broadscale NOAA DHWs
	geomorphic_depth_dhw, dhws = RBAS.spatial_analysis.noaa_dhw_means(geomorphic_depths, config_file)
	geomorphic_depth_dhw[!, :med_depth_5yrs] = dropdims(median(dhws.data[end-5:end, :], dims=1), dims=1)

	# Plot mean and std for dhws
	fig_dhw_mean = spatial_map(geomorphic_depth_dhw, geomorphic_depth_dhw[:, :med_depth_5yrs]; opts=Dict(:colorbar_label => "Median max dhw (2018-2023)", :color_map => :lighttest), geom_col=:geom)
	#save("med_dhw_plot.png", fig_dhw_mean)
end

# ╔═╡ e61dd618-196f-4231-b7fc-238bbff285dc
md"# Select control reefs and sites demo"

# ╔═╡ 3569c9dd-61ac-4680-a980-98dfe0b334fe
md"## Load reef data"

# ╔═╡ c410f452-cbfc-4a7b-aa33-92d7809e94ee
md"###### Load data for reefs in the region of interest."

# ╔═╡ de152600-16d0-4ad6-bafb-284b75c4c590
begin
file_name_reef_data = "C:\\Users\\rcrocker\\Documents\\Github\\ADRIA.jl\\sandbox\\Biodiversity Project\\account_region_geodata.gpkg"
reef_data = GeoDataFrames.read(file_name_reef_data)
end

# ╔═╡ deaa483c-fdfa-4f83-b06b-69e3d69ed3f8
md"## Select project reefs"

# ╔═╡ 44135aba-f0fa-4037-ae63-57069b3fba45
md"###### Select a reef to implement restoration or management activities at."

# ╔═╡ 6a1ae736-0d93-42b2-a224-6fa5387a8770
@bind impact_reef MultiCheckBox(unique(reef_data.GBR_NAME); default=["Moore Reef"])


# ╔═╡ 184dec64-c9e9-4275-8c50-b0ef6b7ef19d
begin
# Calculate composition index
composition_index = ones(length(reef_data.sum_prop_cover))
sp_cover_names = findall(occursin.("species", names(reef_data)))
sp_cover_df = reef_data[:, sp_cover_names]./reef_data.sum_prop_cover
composition_index = composition_index .- sum(Array(sp_cover_df.<0.05).*0.2, dims=2)
reef_data[!,"comp_index"] .= composition_index
end

# ╔═╡ 779aac2e-98da-40ea-81b1-2776d31e0948
md"## Control reef selection criteria"

# ╔═╡ a1d64b3c-6110-40b2-b1e9-87d333086cde
md"###### Select criteria to use to select control reefs which are similar to the selected project reef."

# ╔═╡ 69bbd34c-ee3c-4983-b619-f67fa4f8003e
begin
criteria_key = Dict("mean_dhw_SSP245"=>"Mean dhws", "sum_prop_cover"=>"Historical coral cover", "comp_index"=>"Composition", "mean_cots"=>"Mean COTs", "area_km2"=>"Area (km2)")
@bind control_criteria MultiCheckBox([criteria_key[crit] for crit in keys(criteria_key)], default=["Mean dhws"])
end

# ╔═╡ e12f8d3b-bef1-463c-962d-9110e767d12e
md"""## Reef criteria weightings
###### Select weightings to determine importance of each control reef selection criterium."""

# ╔═╡ aaf93ebb-edd6-4e65-8050-40428329a770
function weightings_input(criteria::Dict)
		
		return combine() do Child
			
			inputs = [
				md""" $(criteria[name]): $(
					Child(name, PlutoUI.Slider(0.0:0.05:1.0, show_value=true, default=1.0))
				)"""
				
				for name in keys(criteria)
			]			
			md"""
			$(inputs)
			"""

		end
end


# ╔═╡ 8fc78db0-9de1-4b58-8aa1-3f654d4650d0
@bind weightings weightings_input(criteria_key)

# ╔═╡ 02923912-1e0e-4a09-80a2-8ffa598bfa2d
md"## Select reef category constraints"

# ╔═╡ d891eed3-0b0f-4159-9562-f4223ef3309a
md"###### Select categorical criteria which must match between control and project reefs."

# ╔═╡ 4c996a9a-c754-4a5d-8606-138e41000f51
begin
	constraints = ["zone_type"=>"Management Zone", "shelf_position"=>"Shelf position"]
	@bind control_constraints MultiCheckBox(constraints, default=["zone_type"])
end

# ╔═╡ 6c279bc4-ce0c-477f-83e4-7b786c4d3661
md"## Rank control reefs"

# ╔═╡ 3f97b011-0aaf-4db1-bafc-d7c0af4221fd
md"###### Get a similarity rating to the project reef for each of the criteria and a measure of aggregate similarity."

# ╔═╡ 8ae0acff-b06f-40df-bee2-e59624126322

function normalize(x)
    return (x .- minimum(x)) ./ (maximum(x) - minimum(x))
end


# ╔═╡ 6b582d7b-300d-4ef5-bcb4-de8fcb3a4e52
function suggest_control_sites(
    impact_site_id::Int64,
    site_data::DataFrame,
    category_constraints::Union{Vector{Symbol},Vector{String}};
    weightings::Matrix{Float64}=ones(1,
        size(site_data, 2) - (1 + length(category_constraints))
    ),
    ID_COLUMN::Union{String,Symbol}=:reef_siteid,
    distance_func=chebyshev)::DataFrame
    impact_site_data = site_data[impact_site_id, :]
    constraints::Vector{Bool} = fill(true, size(site_data, 1))
    for cc in category_constraints
        sites_subset = site_data[:, cc] .== impact_site_data[cc]
        constraints = constraints .& sites_subset
    end

    criteria_df = site_data[constraints, Not(ID_COLUMN, category_constraints...)]
    distances =
        distance_func.(
            Array(impact_site_data[Not(ID_COLUMN, category_constraints...)])',
            Matrix(criteria_df)
        )
    scores = normalize(distances) .* weightings
    distances = dropdims(sum(scores; dims=2); dims=2)
	similarity = normalize(1 .- distances)

    s_order::Vector{Int64} = sortperm(distances; rev=false)
	idx = Int.(findall(constraints)[s_order])
    return DataFrame(
        hcat(idx, site_data[idx, ID_COLUMN], scores[s_order, :], similarity[s_order]),
        vcat(["Index", "Location"], names(criteria_df), ["Similarity"])
    )
end

# ╔═╡ 17aac336-a9e1-47de-bb8e-022a0a4b6121
function get_weights_vec(criteria_cols::Vector{String}, weightings::NamedTuple, n_crit_cols::Int64, n_constraints::Int64)
	weightings_vec = ones(1, n_crit_cols - (1 + n_constraints))
	
	for (w_k, k) in enumerate(criteria_cols)
		weightings_vec[w_k] = getfield(weightings, Symbol(k))
	end
	return weightings_vec
end

# ╔═╡ c88e7632-823f-4908-b55b-41bb54f76f5c
function get_criteria_columns(criteria_key::Dict, criteria::Vector{String})
	return [collect(keys(criteria_key))[collect(values(criteria_key)).==k][1] for k in criteria]
end

# ╔═╡ 0f7b3ace-e7b0-4c9c-acb0-77d919ad4ddf
begin
	impact_site_id = findall(impact_reef[1].==reef_data.GBR_NAME)[1]
	id_col = "GBR_NAME"
	criteria_cols = get_criteria_columns(criteria_key, control_criteria)
	criteria_columns = [id_col, criteria_cols..., control_constraints...]

	weightings_vec = get_weights_vec(criteria_cols, weightings, length(criteria_columns), length(control_constraints))
	control_site_breakdown = suggest_control_sites(impact_site_id, reef_data[:, criteria_columns], control_constraints; weightings=weightings_vec, ID_COLUMN=id_col)

end

# ╔═╡ 3b47f060-ac0c-4ec6-9c54-3f197889a9a7
begin
	sim_control_score = control_site_breakdown[:,"Similarity"]
	sim_control_score = (sim_control_score .- minimum(sim_control_score)) ./ (maximum(sim_control_score) .- minimum(sim_control_score))
	control_potential = -1*ones(length(reef_data.LABEL_ID))
	control_potential[control_site_breakdown[:,"Index"]] .= 1 .- sim_control_score
	fig = spatial_map(reef_data, control_potential; opts=Dict(:color_map => [:gray, :gray, :gray, :gray, :red2, :yellow, :green, :blue2], :colorbar_label=>"Similarity", :color_range => (-1.0,1.0)))

end

# ╔═╡ af05fa21-79f2-4d96-8354-e1867e1493cf
md"## Load site-level data"

# ╔═╡ 97549e30-905e-4196-a051-755cdecbdfdb
md"###### Load data for sites in the region of interest."

# ╔═╡ b9a32d5e-aea6-4c07-8cbe-4e02a7ba896e
begin
file_name_site_data = "C:\\Users\\rcrocker\\Documents\\Github\\rrap-dg-package-creation\\MooreExtDomain_2025-01-21_v070_rc1\\spatial\\MooreExtDomain_2025-01-21_v070_rc1.gpkg"
site_data = GeoDataFrames.read(file_name_site_data)

end

# ╔═╡ 83fd5744-9875-4799-93bf-4b6bbfcbe271
md"## Select project sites"

# ╔═╡ 2a48aa9a-12f6-4a68-9805-ab04c9355c10
md"###### Select sites within the project reef to implement restoration or management activities at according to a set of criteria."

# ╔═╡ 920f586b-5e0f-4da5-a8cd-88fba8eccd95
md"#### Select dominant geomorphic category for project sites"

# ╔═╡ 9a293921-8d85-47f7-abe8-5f95a0defb2a
@bind geomorphic_zone MultiCheckBox(unique(site_data.habitat); default=["Sheltered Slope"])


# ╔═╡ 8d3e1b82-4620-44bd-b2dc-52386c407482
md"#### Select criteria to use to select project sites"

# ╔═╡ a580fcda-9896-45fa-9cc0-6f61abd2edc4
begin
site_data[:, :k_area] = (site_data.k./100).*site_data.area
site_criteria_key = Dict("k_area"=>"Habitable area (km2)", "init_cover"=>"Estimated proportional cover", "depth_med"=>"Median depth (m)", "dhw_mean"=>"Mean max DHW")
@bind impact_site_criteria MultiCheckBox([site_criteria_key[crit] for crit in keys(site_criteria_key)], default=["Mean max DHW"])
end

# ╔═╡ e35cd87a-59c3-4685-86ca-d1e8f010c632
md"#### Select criteria weightings to use to select project sites"

# ╔═╡ 862f5876-5ff5-4ce5-a675-32d60963a97a
@bind impact_site_weightings weightings_input(site_criteria_key)

# ╔═╡ 848f6217-ff82-43e5-aa9b-9af36d64c221
function suggest_impact_sites(
    site_data::DataFrame,
	criteria::Vector{String},
	weights::Matrix{Float64},
	geomorphic_zone::String;
)::DataFrame
	temp_site_data = site_data[site_data.habitat.==geomorphic_zone, :]
    n_locs = size(temp_site_data, 1)
    locations = temp_site_data.site_id
	locations_inds = collect(1:n_locs)
    geomorphic_protection = zeros(n_locs)

    # geom_vec = site_data.habitat
    # geomorphic_protection[geom_vec .== "Sheltered Reef Slope"] .= 3
    # geomorphic_protection[geom_vec .== "Back Reef Slope"] .= 2
    # geomorphic_protection[geom_vec .== "Reef Slope"] .= 1
    # geomorphic_protection[geom_vec .== "Deep Lagoon"] .= 1

    # site_data[!, "geom_protect"] .= geomorphic_protection

    norm_mat = normalize(Matrix(temp_site_data[:, criteria]))
	scores = norm_mat .* weights
    scores = normalize(dropdims(sum(scores; dims=2); dims=2))

    s_order::Vector{Int64} = sortperm(scores; rev=true)

    return DataFrame(
        hcat(locations_inds[s_order], locations[s_order], scores[s_order]), ["Index", "Sites", "Rating"])
end

# ╔═╡ 8a9ff5f4-b6df-413b-9ab9-8eca97c29c60
begin

	criteria_cols_sites = get_criteria_columns(site_criteria_key, impact_site_criteria)
	weightings_impact_site = get_weights_vec(criteria_cols_sites, impact_site_weightings, length(criteria_cols_sites), -1)
	reef_idx = occursin.(split(impact_reef[1]," ")[1], site_data.Reef)
	project_site_ranks = suggest_impact_sites(site_data[reef_idx,:], criteria_cols_sites, weightings_impact_site, geomorphic_zone[1])

end

# ╔═╡ 145b2fba-4051-4c0d-95e7-18ad1a4fcedb
begin
	sites_impact_score = zeros(size(site_data,1))
	sites_impact_score[project_site_ranks.Index] .= project_site_ranks.Rating
	fig_impact_sites = spatial_map(site_data, sites_impact_score; opts=Dict(:color_map => Reverse(:batlow), :colorbar_label=>"Project site rating", :color_range => (0.0,1.0)),geom_col=:geom)

end

# ╔═╡ 9a6b1f6f-d465-4e7b-b182-0e988ebe30d3
md"## Select control sites"

# ╔═╡ 0f0ac8a3-43ae-471c-b02b-4c9fab6247cf
md"###### Select sites within the project area which are sufficiently similar to the selected project sites as control sites."

# ╔═╡ 9188f844-7f91-49cf-899b-f3ad6b526fae
md"#### Select criteria to use to select control sites"

# ╔═╡ 9557bcb0-d9b4-4632-82d7-0a16ae4e8421
begin
@bind control_site_criteria MultiCheckBox([site_criteria_key[crit] for crit in keys(site_criteria_key)], default=["Mean max DHW"])
end

# ╔═╡ 3c758242-eb4d-45e9-a81c-8f262dcdbdb1
md"#### Select weightings for control site criteria"

# ╔═╡ e61d72b7-f869-46f7-bf94-842dcdca0e0c
@bind control_site_weightings weightings_input(site_criteria_key)

# ╔═╡ bb257e2d-0470-4cd3-8388-ccd9f250f9af
md"#### Number of project sites to use"

# ╔═╡ a672759a-5f32-4305-b0f0-ef41f68c7aac
@bind n_impact NumberField(1:5, default=3)

# ╔═╡ 0d1a3f5e-1af6-4c02-bb9f-4aedc445816b
md"#### Number of control sites to choose for each project site"

# ╔═╡ 32285853-9c00-4394-94ae-54fc3f8d9b7f
@bind n_control NumberField(1:5, default=3)

# ╔═╡ 99f37fe6-73aa-4cac-86ad-311e3c5dcbef
begin
	select_controls = Dict()
	non_impact_reef = .!occursin.(split(impact_reef[1]," ")[1], site_data.Reef)
	n_impact_min = minimum([length(project_site_ranks.Index),n_impact])
 	impact_ids = project_site_ranks.Index[1:n_impact_min]
 	id_col_controls = "site_id"
 	control_criteria_cols = get_criteria_columns(site_criteria_key, control_site_criteria)
 	control_criteria_columns = [id_col_controls, "habitat", control_criteria_cols...]
 	weightings_vec_control = get_weights_vec(control_criteria_cols, control_site_weightings, length(control_criteria_columns), 1)
	for impact_idx in impact_ids
 		select_controls[string(impact_idx)] = suggest_control_sites(impact_ids[1], site_data[non_impact_reef, control_criteria_columns], ["habitat"]; weightings=weightings_vec_control, ID_COLUMN=id_col_controls)
	end

end

# ╔═╡ 002f0120-0265-4562-b44f-b6584e6ffee6
md"#### Select project site to display suggested control sites for"

# ╔═╡ eadb30dc-483c-423e-a7fa-40d8f7c3575b
@bind display_control_sites Select(project_site_ranks.Sites[1:n_impact_min])

# ╔═╡ ca4db614-bc7f-4ad7-99d0-b79e86d9d428
begin
	display_idx = project_site_ranks.Index[project_site_ranks.Sites.==display_control_sites]
	collect(keys(select_controls))
	temp_control_df = select_controls[string(display_idx[1])]
end

# ╔═╡ 740a738f-3de7-4b68-a9a7-729951c96342
begin
	sites_control_score = zeros(size(site_data,1))
	sites_control_score[temp_control_df.Index] .= temp_control_df.Similarity
	fig_control_sites = spatial_map(site_data, sites_control_score; opts=Dict(:color_map => Reverse(:batlow), :colorbar_label=>"Site similarity rating", :color_range => (0.0,1.0)),geom_col=:geom)
end

# ╔═╡ Cell order:
# ╠═8d27ffec-e4d5-46bc-8218-b5a44adbb796
# ╟─a2529343-55a8-43ea-b066-f928d52ee534
# ╟─570139e5-cd28-4689-a2c8-0a3682f5e572
# ╟─432ec934-ebaa-430b-9401-e64369a9ecc7
# ╟─4c135067-265b-4e98-bb31-a6daa9d426ad
# ╟─2a5cf51c-8831-4aff-8419-3603d42030e3
# ╟─646613cb-740c-4ef7-b406-5cdd31d9dcc4
# ╟─89475ad4-0102-4ea2-bbf0-497da8d935a2
# ╟─df182f3d-29d2-4b3f-bd45-391a8a2dfee9
# ╟─bfaed850-ecd4-4319-8ba4-1c48670a0a85
# ╟─a0adc075-5289-45dc-b23e-b681d2785943
# ╟─4a9ae439-da76-4266-98e0-849eae994020
# ╟─4877ea3e-8bf9-4c08-9694-14864bfb7ad9
# ╟─05aa9fcd-12ab-48db-a657-827366087583
# ╟─844f575d-1137-45e9-9d9a-fd84c0e71c0d
# ╟─3aa16739-3c70-4895-82e0-5341da5b9956
# ╟─0e11066d-8afa-48a1-9344-c2d2b7c98223
# ╟─e61dd618-196f-4231-b7fc-238bbff285dc
# ╟─3569c9dd-61ac-4680-a980-98dfe0b334fe
# ╟─c410f452-cbfc-4a7b-aa33-92d7809e94ee
# ╟─de152600-16d0-4ad6-bafb-284b75c4c590
# ╟─deaa483c-fdfa-4f83-b06b-69e3d69ed3f8
# ╟─44135aba-f0fa-4037-ae63-57069b3fba45
# ╟─6a1ae736-0d93-42b2-a224-6fa5387a8770
# ╟─184dec64-c9e9-4275-8c50-b0ef6b7ef19d
# ╟─779aac2e-98da-40ea-81b1-2776d31e0948
# ╟─a1d64b3c-6110-40b2-b1e9-87d333086cde
# ╟─69bbd34c-ee3c-4983-b619-f67fa4f8003e
# ╟─e12f8d3b-bef1-463c-962d-9110e767d12e
# ╟─aaf93ebb-edd6-4e65-8050-40428329a770
# ╟─8fc78db0-9de1-4b58-8aa1-3f654d4650d0
# ╟─02923912-1e0e-4a09-80a2-8ffa598bfa2d
# ╟─d891eed3-0b0f-4159-9562-f4223ef3309a
# ╟─4c996a9a-c754-4a5d-8606-138e41000f51
# ╟─6c279bc4-ce0c-477f-83e4-7b786c4d3661
# ╟─3f97b011-0aaf-4db1-bafc-d7c0af4221fd
# ╟─6b582d7b-300d-4ef5-bcb4-de8fcb3a4e52
# ╟─8ae0acff-b06f-40df-bee2-e59624126322
# ╟─17aac336-a9e1-47de-bb8e-022a0a4b6121
# ╟─c88e7632-823f-4908-b55b-41bb54f76f5c
# ╟─0f7b3ace-e7b0-4c9c-acb0-77d919ad4ddf
# ╟─3b47f060-ac0c-4ec6-9c54-3f197889a9a7
# ╟─af05fa21-79f2-4d96-8354-e1867e1493cf
# ╟─97549e30-905e-4196-a051-755cdecbdfdb
# ╟─b9a32d5e-aea6-4c07-8cbe-4e02a7ba896e
# ╟─83fd5744-9875-4799-93bf-4b6bbfcbe271
# ╟─2a48aa9a-12f6-4a68-9805-ab04c9355c10
# ╟─920f586b-5e0f-4da5-a8cd-88fba8eccd95
# ╟─9a293921-8d85-47f7-abe8-5f95a0defb2a
# ╟─8d3e1b82-4620-44bd-b2dc-52386c407482
# ╟─a580fcda-9896-45fa-9cc0-6f61abd2edc4
# ╟─e35cd87a-59c3-4685-86ca-d1e8f010c632
# ╟─862f5876-5ff5-4ce5-a675-32d60963a97a
# ╟─848f6217-ff82-43e5-aa9b-9af36d64c221
# ╟─8a9ff5f4-b6df-413b-9ab9-8eca97c29c60
# ╟─145b2fba-4051-4c0d-95e7-18ad1a4fcedb
# ╟─9a6b1f6f-d465-4e7b-b182-0e988ebe30d3
# ╟─0f0ac8a3-43ae-471c-b02b-4c9fab6247cf
# ╟─9188f844-7f91-49cf-899b-f3ad6b526fae
# ╟─9557bcb0-d9b4-4632-82d7-0a16ae4e8421
# ╟─3c758242-eb4d-45e9-a81c-8f262dcdbdb1
# ╟─e61d72b7-f869-46f7-bf94-842dcdca0e0c
# ╟─bb257e2d-0470-4cd3-8388-ccd9f250f9af
# ╟─a672759a-5f32-4305-b0f0-ef41f68c7aac
# ╟─0d1a3f5e-1af6-4c02-bb9f-4aedc445816b
# ╟─32285853-9c00-4394-94ae-54fc3f8d9b7f
# ╟─99f37fe6-73aa-4cac-86ad-311e3c5dcbef
# ╟─002f0120-0265-4562-b44f-b6584e6ffee6
# ╟─eadb30dc-483c-423e-a7fa-40d8f7c3575b
# ╟─ca4db614-bc7f-4ad7-99d0-b79e86d9d428
# ╟─740a738f-3de7-4b68-a9a7-729951c96342
