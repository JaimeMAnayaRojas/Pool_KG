
using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall
using Plots
using Plots.Measures
using LaTeXStrings
using KernelDensity


# load my functions
include("Functions.jl")
readdir("outputs")
post_SurvG = CSV.read("outputs/Post_Survival_G.csv", DataFrame)# Statistical analyses
post_GrowthG = CSV.read("outputs/Post_Growth_G.csv", DataFrame)# Statistical analyses
post_ReprG = CSV.read("outputs/Post_Repr_G.csv", DataFrame)# Statistical analyses

post_SurvK = CSV.read("outputs/Post_Survival_K.csv", DataFrame)# Statistical analyses
post_GrowthK = CSV.read("outputs/Post_Growth_K.csv", DataFrame)# Statistical analyses
post_ReprK = CSV.read("outputs/Post_Repr_K.csv", DataFrame)# Statistical analyses


readdir("data")
DataG = CSV.read("data/GuppyIPM.csv", DataFrame);
DataK = CSV.read("data/KillifishIPM.csv", DataFrame);



a = filter(:KG => x -> x == 1, DataG)
pSG = histogram(a.SL1_mm, label = "KG", bins= 50, alpha = 0.5, 
title = "a) Guppy", titlefont = font(10), titleloc = :left)
a = filter(:NK => x -> x == 1, DataG)
histogram!(a.SL1_mm, label = "NK", bins= 50, alpha = 0.5)

a = filter(:KG => x -> x == 1, DataK)
pSK = histogram(a.SL1_mm, label = "KG", titlefont = font(10),  
bins= 50, title = "b) Killifish", titleloc = :left, alpha = 0.5)
a = filter(:NG => x -> x == 1, DataK)
histogram!(a.SL1_mm, label = "NG",  bins= 50, alpha = 0.5)

plot(pSG, pSK, layout = (2,1))
xlabel!("Size (mm)")
ylabel!("Frequency (N)")





########333


function p_link(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NK::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βz2= df.b_z2[row]
    βNK= df."b_NK"[row]
    βzNK= df."b_z.NK"[row]
    μ = α .+ βNK .* NK .+  (βz .+ βzNK .* NK) .* zc .+ βz2 .* z2  # linear predictor
    return(μ)
end



function p_linkK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    βz2= df.b_z2[row]
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βNG= df."b_NG"[row]
    βzNG= df."b_z.NG"[row]
    μ = α .+ βNG .* NG .+  (βz .+ βzNG .* NG) .* zc .+ βz2 .* z2 # linear predictor
    return(μ)
end



function G_linkK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    βz1= df.b_bz1_Intercept[row]
    βz2= df.b_bz2_Intercept[row]
    α= df.b_b0_Intercept[row]
    βNG= df."b_b0_NG"[row]
    βz1NG= df."b_bz1_NG"[row]
    βz2NG= df."b_bz2_NG"[row]
    ω = df.b_omega_Intercept[row]

    μ = zeros(length(z))

    for j in 1:length(μ)
        if z[j] < ω
            μ[j] = α .+ βNG .* NG .+ (βz1 .+ βz1NG .* NG) .* (z[j] .- ω) # linear predictor
        else
            μ[j] = α .+ βNG .* NG .+ (βz2 .+ βz2NG .* NG) .* (z[j] .- ω) # linear predictor
        end
    
    end

    
    return(μ)
end



include("Figures_Guppy.jl")
savefig("Figure-1.png")



include("Figures_Killifish.jl")
savefig("Figure-2.png")

