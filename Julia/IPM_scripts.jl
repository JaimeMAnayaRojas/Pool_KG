#using Distributed
# Add four processes to use for sampling.
#addprocs(4)
#@everywhere begin
using Distributions
using CSV
using DataFrames
using LinearAlgebra
using Turing
using DynamicHMC
using StatsPlots
using JLD2

#addprocs(Nchains)


function which_R(x, find)
	x[ismissing.(x)] .= 0
	return(findall(x .==find))
end

cd("$(homedir())")
cd("Dropbox/Projects_JM/FSU/Pool_manipulation/Pool_R/")
pwd()

# Statistical analyses
# As a way to better learn and use Julia, I will do the stats modeling in Turing

## First, get the data for guppies

guppy_df = DataFrame(CSV.File("data/GuppyIMP.csv"))
kill_df = DataFrame(CSV.File("data/KillifishIPM.csv"))


# use only the data for females
guppy_df= guppy_df[findall(guppy_df.Sex1 .!= "M" ),:]
kill_df= guppy_df[findall(kill_df.Sex1 .!= "M" ),:]


Ggrow_indx = which_R(guppy_df.surv, 1) # find the surviving fish
Grepr_indx = which_R(guppy_df.Repr, 1)

Kgrow_indx = which_R(kill_df.surv, 1) # find the surviving fish
Krepr_indx = which_R(kill_df.Repr, 1)

guppy_df.Location = convert(Vector{String}, guppy_df.Location)
guppy_df.stream = zeros(length(guppy_df.Location))
guppy_df.stream[findall(guppy_df.Location .== "CAI")] .= 1
guppy_df.stream[findall(guppy_df.Location .== "NAR")] .= 2
guppy_df.stream[findall(guppy_df.Location .== "QUA")] .= 3
guppy_df.stream[findall(guppy_df.Location .== "QUA2")] .= 4


kill_df.Location = convert(Vector{String}, kill_df.Location)

unique(kill_df.Location)

kill_df.stream = zeros(length(kill_df.Location))
kill_df.stream[findall(kill_df.Location .== "CAI")] .= 1
kill_df.stream[findall(kill_df.Location .== "NAR")] .= 2
kill_df.stream[findall(kill_df.Location .== "QUA")] .= 3
#kill_df.stream[findall(kill_df.Location .== "QUA2")] .= 4

guppy_df.stream = convert.(Int64, guppy_df.stream)
kill_df.stream = convert.(Int64, kill_df.stream)
unique(guppy_df.stream )

guppy_df.canopy = guppy_df.canopy .- mean(unique(guppy_df.canopy))
guppy_df.area = guppy_df.area .- mean(unique(guppy_df.area))

kill_df.canopy = kill_df.canopy .- mean(unique(kill_df.canopy))
kill_df.area = kill_df.area .- mean(unique(kill_df.area))

# make Turing model
@model function Fish_model(Surv_G, size_surv_G, NK_surv_G, stream_surv_G, canopy_surv_G, area_surv_G,
						z1_G, size_grow_G, NK_grow_G, stream_grow_G, Repr_G, canopy_grow_G, area_grow_G,
						Fec_G, size_fec_G, NK_fec_G, stream_fec_G, canopy_fec_G, area_fec_G, 
						
						Surv_K, size_surv_K, NG_surv_K, stream_surv_K, canopy_surv_K, area_surv_K,
						z1_K, size_grow_K, NG_grow_K, stream_grow_K, Repr_K, canopy_grow_K, area_grow_K,
						Fec_K, size_fec_K, NG_fec_K, stream_fec_K, canopy_fec_K, area_fec_K)
	
	## Guppy models
	# Survival
	z_surv_G = size_surv_G .- 10.0
	zNK_surv_G = z_surv_G .* NK_surv_G
	σₛ_G ~ truncated(Cauchy(0, 2), 0, Inf)
	αₛ_surv_G ~ filldist(Normal(0, σₛ_G), length(unique(stream_surv_G)))
	α_surv_G ~ Normal(0, 1)
    βₖ_surv_G ~ Normal(0, 1)
	βz_surv_G ~ Normal(0, 1)
	βc_surv_G ~ Normal(0, 1)
	βa_surv_G ~ Normal(0, 1)
	βzNK_surv_G ~ Normal(0, 1)
	p_surv_G = α_surv_G .+ βₖ_surv_G .* NK_surv_G .+ βz_surv_G .* z_surv_G .+
			 βzNK_surv_G .* zNK_surv_G .+ αₛ_surv_G[stream_surv_G] .+ βc_surv_G*canopy_surv_G .+ βa_surv_G.*area_surv_G
    Surv_G .~ BinomialLogit.(1, p_surv_G)

	# Growth

	z_grow_G = size_grow_G .- 10.0
	zNK_grow_G = z_grow_G .* NK_grow_G
	σ_grow_G ~ truncated(Cauchy(0, 2), 0, Inf)
	#σₛ_grow ~ truncated(Cauchy(0, 1), 0, Inf)
	αₛ_grow_G ~ filldist(Normal(0, σₛ_G), length(unique(stream_grow_G)))
	α_grow_G ~ Normal(10, 10)
	βₖ_grow_G ~ Normal(0, 10)
	βz_grow_G ~ Normal(0, 10)
	βzNK_grow_G ~ Normal(0, 10)
	βc_grow_G ~ Normal(0, 10)
	βa_grow_G ~ Normal(0, 10)
	μ_grow_G = α_grow_G .+ βₖ_grow_G .* NK_grow_G .+ βz_grow_G.*z_grow_G .+
			βzNK_grow_G .* zNK_grow_G .+ αₛ_grow_G[stream_grow_G] .+ βc_grow_G*canopy_grow_G .+ βa_grow_G.*area_grow_G
	z1_G .~ Normal.(μ_grow_G, σ_grow_G)


	# Reproduction
	#σₛ_rep ~ truncated(Cauchy(0, 1), 0, Inf)
	αₛ_rep_G ~ filldist(Normal(0, σₛ_G), length(unique(stream_grow_G)))
	α_rep_G ~ Normal(0, 10)
	βₖ_rep_G ~ Normal(0, 10)
	βz_rep_G ~ Normal(0, 10)
	βzNK_rep_G ~ Normal(0, 10)
	βc_rep_G ~ Normal(0, 10)
	βa_rep_G ~ Normal(0, 10)
	p_rep_G = α_rep_G .+ βₖ_rep_G .* NK_grow_G .+ βz_rep_G .*z_grow_G .+
	            βzNK_rep_G .* zNK_grow_G .+ αₛ_rep_G[stream_grow_G] .+ βc_rep_G*canopy_grow_G .+ βa_rep_G.*area_grow_G
	Repr_G .~ BinomialLogit.(1, p_rep_G)

	# Fecundity
	z_fec_G = size_fec_G .- 10.0
	zNK_fec_G = z_fec_G .* NK_fec_G
	αₛ_fec_G ~ filldist(Normal(0, σₛ_G), length(unique(stream_fec_G)))
	α_fec_G ~ Normal(0, 10)
	βₖ_fec_G ~ Normal(0, 10)
	βz_fec_G ~ Normal(0, 10)
	βzNK_fec_G ~ Normal(0, 10)
	βc_fec_G ~ Normal(0, 10)
	βa_fec_G ~ Normal(0, 10)
	λ_fec_G = α_fec_G .+ βₖ_fec_G .* NK_fec_G .+ βz_fec_G .*z_fec_G  .+
	 		βzNK_fec_G .* zNK_fec_G .+ αₛ_fec_G[stream_fec_G] .+ βc_fec_G *canopy_fec_G .+ βa_fec_G.*area_fec_G
	Fec_G .~ Poisson.(exp.(λ_fec_G))
	## Kilifish models
	# Survival
	z_surv_K = size_surv_K .- 10.0
	zNG_surv_K = z_surv_K .* NG_surv_K
	σₛ_K ~ truncated(Cauchy(0, 2), 0, Inf)
	αₛ_surv_K ~ filldist(Normal(0, σₛ_K), length(unique(stream_surv_K)))
	α_surv_K ~ Normal(0, 1)
    βₖ_surv_K ~ Normal(0, 1)
	βz_surv_K ~ Normal(0, 1)
	βc_surv_K ~ Normal(0, 1)
	βa_surv_K ~ Normal(0, 1)
	βzNG_surv_K ~ Normal(0, 1)
	p_surv_K = α_surv_K .+ βₖ_surv_K .* NG_surv_K .+ βz_surv_K .* z_surv_K .+
			 βzNG_surv_K .* zNG_surv_K .+ αₛ_surv_K[stream_surv_K] .+ βc_surv_K*canopy_surv_K .+ βa_surv_K.*area_surv_K
    Surv_K .~ BinomialLogit.(1, p_surv_K)

	# Growth

	z_grow_K = size_grow_K .- 10.0
	zNG_grow_K = z_grow_K .* NG_grow_K
	σ_grow_K ~ truncated(Cauchy(0, 2), 0, Inf)
	#σₛ_grow ~ truncated(Cauchy(0, 1), 0, Inf)
	αₛ_grow_K ~ filldist(Normal(0, σₛ_K), length(unique(stream_grow_K)))
	α_grow_K ~ Normal(10, 10)
	βₖ_grow_K ~ Normal(0, 10)
	βz_grow_K ~ Normal(0, 10)
	βzNG_grow_K ~ Normal(0, 10)
	βc_grow_K ~ Normal(0, 10)
	βa_grow_K ~ Normal(0, 10)
	μ_grow_K = α_grow_K .+ βₖ_grow_K .* NG_grow_K .+ βz_grow_K.*z_grow_K .+
			βzNG_grow_K .* zNG_grow_K .+ αₛ_grow_K[stream_grow_K] .+ βc_grow_K*canopy_grow_K .+ βa_grow_K.*area_grow_K
	z1_K .~ Normal.(μ_grow_K, σ_grow_K)


	# Reproduction
	#σₛ_rep ~ truncated(Cauchy(0, 1), 0, Inf)
	αₛ_rep_K ~ filldist(Normal(0, σₛ_K), length(unique(stream_grow_K)))
	α_rep_K ~ Normal(0, 10)
	βₖ_rep_K ~ Normal(0, 10)
	βz_rep_K ~ Normal(0, 10)
	βzNG_rep_K ~ Normal(0, 10)
	βc_rep_K ~ Normal(0, 10)
	βa_rep_K ~ Normal(0, 10)
	p_rep_K = α_rep_K .+ βₖ_rep_K .* NG_grow_K .+ βz_rep_K .*z_grow_K .+
	            βzNG_rep_K .* zNG_grow_K .+ αₛ_rep_K[stream_grow_K] .+ βc_rep_K*canopy_grow_K .+ βa_rep_K.*area_grow_K
	Repr_K .~ BinomialLogit.(1, p_rep_K)

	# Fecundity
	z_fec_K = size_fec_K .- 10.0
	zNG_fec_K = z_fec_K .* NG_fec_K
	αₛ_fec_K ~ filldist(Normal(0, σₛ_K), length(unique(stream_fec_K)))
	α_fec_K ~ Normal(0, 10)
	βₖ_fec_K ~ Normal(0, 10)
	βz_fec_K ~ Normal(0, 10)
	βzNG_fec_K ~ Normal(0, 10)
	βc_fec_K ~ Normal(0, 10)
	βa_fec_K ~ Normal(0, 10)
	λ_fec_K = α_fec_K .+ βₖ_fec_K .* NG_fec_K .+ βz_fec_K .*z_fec_K  .+
	 		βzNG_fec_K .* zNG_fec_K .+ αₛ_fec_K[stream_fec_K] .+ βc_fec_K *canopy_fec_K .+ βa_fec_K.*area_fec_K
	Fec_K .~ Poisson.(exp.(λ_fec_K))


end

model = Fish_model(guppy_df.surv, guppy_df.SL1_mm, guppy_df.NK, guppy_df.stream, guppy_df.canopy, guppy_df.area,
				guppy_df.SL2_mm[Ggrow_indx], guppy_df.SL1_mm[Ggrow_indx], guppy_df.NK[Ggrow_indx], guppy_df.stream[Ggrow_indx], guppy_df.canopy[Ggrow_indx], guppy_df.area[Ggrow_indx],
				guppy_df.Repr[Ggrow_indx],
				guppy_df.Recr[Grepr_indx], guppy_df.SL1_mm[Grepr_indx], guppy_df.NK[Grepr_indx], guppy_df.stream[Grepr_indx], guppy_df.canopy[Grepr_indx], guppy_df.area[Grepr_indx],
				kill_df.surv, kill_df.SL1_mm, kill_df.NG, kill_df.stream, kill_df.canopy, kill_df.area,
				kill_df.SL2_mm[Kgrow_indx], kill_df.SL1_mm[Kgrow_indx], kill_df.NG[Kgrow_indx], kill_df.stream[Kgrow_indx], kill_df.canopy[Kgrow_indx], kill_df.area[Kgrow_indx],
				kill_df.Repr[Kgrow_indx],
				kill_df.Recr[Krepr_indx], kill_df.SL1_mm[Krepr_indx], kill_df.NG[Krepr_indx], kill_df.stream[Krepr_indx], kill_df.canopy[Krepr_indx], kill_df.area[Krepr_indx]
				)


# # Sample four chains using multiple threads, each with 2000 samples.

@time chain = sample(model, NUTS(),  1000)
#@time chain = sample(model, NUTS(), MCMCThreads(), 4000, 4)

@save "modeL.jld2" chain
chain = load_object("modeL.jld2")


include("Julia/IPM_G.jl")

