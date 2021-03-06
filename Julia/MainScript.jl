
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
post_SurvG = CSV.read("outputs/Post_Survival_G.csv", DataFrame)[:,2:end]# Statistical analyses
post_GrowthG = CSV.read("outputs/Post_Growth_G.csv", DataFrame)[:,2:end]# Statistical analyses
post_ReprG = CSV.read("outputs/Post_Repr_G.csv", DataFrame)[:,2:end]# Statistical analyses


# summary tables

CSV.write("outputs/Sum_SurvG.csv", Post_summary(post_SurvG))

# How much increaseing fish biomass changes survival?

names(post_SurvG)

100 * HDI(My_Logistic.(post_SurvG.b_Intercept .+ post_SurvG.b_FishBiom) .- My_Logistic.(post_SurvG.b_Intercept))
mean(My_Logistic.(post_SurvG.b_Intercept))
100 * HDI(My_Logistic.(post_SurvG.b_Intercept))

100 *mean(My_Logistic.(post_SurvG.b_Intercept .+ post_SurvG.b_FishBiom))
100 * HDI(My_Logistic.(post_SurvG.b_Intercept .+ post_SurvG.b_FishBiom))




CSV.write("outputs/Sum_GrowthG.csv", Post_summary(post_GrowthG))


a = 18 .* exp.(post_GrowthG.b_Intercept) .- 18
b =  18 .* exp.(post_GrowthG.b_Intercept .+ post_GrowthG.b_Density) .- 18
mean(b - a)
HDI(b - a)


a = 18 .* exp.(post_GrowthG.b_Intercept) .- 18
b =  18 .* exp.(post_GrowthG.b_Intercept .+ post_GrowthG.b_canopy) .- 18
mean(b - a)
HDI(b - a)



CSV.write("outputs/Sum_FecuG.csv", Post_summary(post_ReprG))




post_SurvK = CSV.read("outputs/Post_Survival_K.csv", DataFrame)[:,2:end]# Statistical analyses
post_GrowthK = CSV.read("outputs/Post_Growth_K.csv", DataFrame)[:,2:end]# Statistical analyses
post_ReprK = CSV.read("outputs/Post_Repr_K.csv", DataFrame)[:,2:end]# Statistical analyses

CSV.write("outputs/Sum_SurvK.csv", Post_summary(post_SurvK))
CSV.write("outputs/Sum_GrowthK.csv", Post_summary(post_GrowthK))
CSV.write("outputs/Sum_FecuK.csv", Post_summary(post_ReprK))


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


function G_link(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NK::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    ??= df.b_Intercept[row]
    ??z= df.b_z[row]
    ??z2= df.b_z2[row]
    ??NK= df."b_NK"[row]
    ??zNK= df."b_z.NK"[row]
    ?? = ?? .+ ??NK .* NK .+  (??z .+ ??zNK .* NK) .* zc .+ ??z2 .* z2  # linear predictor
    return(??)
end



function p_linkK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    ??z2= df.b_z2[row]
    ??= df.b_Intercept[row]
    ??z= df.b_z[row]
    ??NG= df."b_NG"[row]
    ??zNG= df."b_z.NG"[row]
    ?? = ?? .+ ??NG .* NG .+  (??z .+ ??zNG .* NG) .* zc .+ ??z2 .* z2 # linear predictor
    return(??)
end



function G_linkK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    ??z1= df.b_bz1_Intercept[row]
    ??z2= df.b_bz2_Intercept[row]
    ??= df.b_b0_Intercept[row]
    ??NG= df."b_b0_NG"[row]
    ??z1NG= df."b_bz1_NG"[row]
    ??z2NG= df."b_bz2_NG"[row]
    ?? = df.b_omega_Intercept[row]

    ?? = zeros(length(z))

    for j in 1:length(??)
        if z[j] < ??
            ??[j] = ?? .+ ??NG .* NG .+ (??z1 .+ ??z1NG .* NG) .* (z[j] .- ??) # linear predictor
        else
            ??[j] = ?? .+ ??NG .* NG .+ (??z2 .+ ??z2NG .* NG) .* (z[j] .- ??) # linear predictor
        end
    
    end

    
    return(??)
end



include("Figures_Guppy.jl")
savefig("plots/Figure-1.svg")
savefig("plots/Figure-1.png")


include("Figures_Killifish.jl")
savefig("plots/Figure-2.png")
savefig("plots/Figure-2.svg")

