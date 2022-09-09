using Distributions
using CSV
using DataFrames
using LinearAlgebra
using StatsPlots
using JLD2
using RCall
using Plots
using Plots.Measures
using CategoricalArrays

Gz = CSV.read("data/01_GuppyField_data.csv", DataFrame)
Kz = CSV.read( "data/01_KillifishField_data.csv", DataFrame)




#

println(describe(Gz))

replace!(Gz.Mark, missing=>"NaN")
replace!(Kz.Mark, missing=>"NaN")

filter!(:Mark => x -> x !="others", Gz)
filter!(:Mark => x -> x !="others", Kz)

Gz.Mark = categorical(Gz.Mark)
Kz.Mark = categorical(Kz.Mark)


replace!(Gz.SL1_mm, missing=>NaN)
replace!(Gz.SL2_mm, missing=>NaN)



replace!(Kz.SL1_mm, missing=>NaN)
replace!(Kz.SL2_mm, missing=>NaN)



Gz.Location = categorical(Gz.Location)
filter!(:Location => x -> x !="JOE", Gz)
filter!(:Location => x -> x !="JOE", Kz)

### Make data for recruitment

#

Gz.R1 .= 1
Gz.R2 .= 0

#Gz.R1[findall(Gz.SL1_mm .> 14.)] .= 0
Gz.R2[findall(Gz.Mark .== NaN)] .= 1

isna.(Gz.Mark)
Kz.R1 .= 1
Kz.R2 .= 1

Kz.R1[findall(Kz.SL1_mm .> 14.)] .= 0
Kz.R2[findall(Kz.SL2_mm .> 14.)] .= 0


println(names(Gz))

sGz = combine(groupby(Gz, [:Location, :NK]), [:R1, :R2] .=> [sum, sum] .=> [:Recr_1, :Recr_2])
sKz = combine(groupby(Kz, [:Location, :NG]), [:R1, :R2] .=> [sum, sum] .=> [:Recr_1, :Recr_2])


DataG = CSV.read("data/GuppyIPM.csv", DataFrame);
DataK = CSV.read("data/KillifishIPM.csv", DataFrame);

# Get the environmental data only
EnvDat = vcat(DataG[:, [:Location, :sp, :Pool_1, :KG, :NK, :NG, :area, :BiomassG1, :BiomassK1]],DataK[:, [:Location, :sp, :Pool_1, :KG, :NK, :NG, :area, :BiomassG1, :BiomassK1]])


EnvDat.FishBiomass = EnvDat.BiomassG1 .+ EnvDat.BiomassK1;
EnvDat.Density = EnvDat.FishBiomass ./ EnvDat.area;






EnvDat = outerjoin(EnvDat, 
combine(groupby(EnvDat, [:Location]), [:Density, :Density, :area, :area] .=> [mean, std,mean, std]), 
        on = :Location
)


EnvDat.Density_s = (EnvDat.Density .- EnvDat.Density_mean) ./ EnvDat.Density_std;  
EnvDat.Area_s = (EnvDat.area .- EnvDat.area_mean) ./ EnvDat.area_std;  

EnvDat



filter(:NG => x -> x != 1, EnvDat )

a = combine(groupby(filter(:NG => x -> x != 1, EnvDat ), [:Location, :NK]), [:Density_s, :Area_s] .=> [mean, mean] .=> [:Density, :Area])
GupRecr = outerjoin(sGz, a, on = [:Location, :NK])


a = combine(groupby(filter(:NK => x -> x != 1, EnvDat ), [:Location, :NG]), [:Density_s, :Area_s] .=> [mean, mean] .=> [:Density, :Area])

KillRecr = outerjoin(sKz, a, on = [:Location, :NG])


# Make recruitment plots for guppies and killifish

rename!(GupRecr, :NK => :Removal)
GupRecr.Sp .= "G"

rename!(KillRecr, :NG => :Removal)
KillRecr.Sp .= "K"
RecrData = vcat(GupRecr, KillRecr)

RecrData.RCh = log.(RecrData.Recr_2 ./ RecrData.Recr_1) 


include("Functions.jl")
@rput GupRecr 
@rput KillRecr

R"""
library("brms")

GupRecr$Removal = factor(GupRecr$Removal)
GupRecr$z = GupRecr$Recr_1 - 50
mG <- brm( Recr_2 ~ z * Removal + Density + Area, family = negbinomial(), GupRecr)
post_mG = posterior_samples(mG)

KillRecr$Removal = factor(KillRecr$Removal)
KillRecr$z = KillRecr$Recr_1 - 50
mK <- brm( Recr_2 ~ z * Removal + Density + Area, family = negbinomial(), KillRecr)
post_mK = posterior_samples(mK)


"""

@rget post_mG


R"""
summary(m1)
"""
println(Post_summary(postRec))


postRec



#
scatter(RecrData.Recr_1, RecrData.Recr_2, legend =:left)
plot!([0,300], [0, 300])


scatter(GupRecr.Recr_1, GupRecr.Recr_2, groups = GupRecr.NK)
plot!([0,80], [0, 80])



# Caigual
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1a = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "Before (105)",
 linewidth = 1.5, linestyle = :dot)

 length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
 histogram!(p.SL2_mm, bins = 10, alpha=0.2, color = :red, label = "End (109)",
 linewidth = 1.5, linestyle = :dash)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "Introduced (82)")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "Recovered (51)", 
        titlefont = font(10),  
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")
ylims!(0,50)


# Caigual NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1b = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "103",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "170",
        linewidth = 1.5, linestyle = :dash)
 

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "92")

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "57", 
        titlefont = font(10),  
        title = "b) Caigual-NK", titleloc = :left
)

# Naranjo KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1c = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "90",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "93",
 linewidth = 1.5, linestyle = :dash)



p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "51")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "25", 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Naranjo NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1d = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "34",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "84",
 linewidth = 1.5, linestyle = :dash)
 

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "25")

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "20", 
        titlefont = font(10),  
        title = "d) Naranjo-NK", titleloc = :left
)

## Quare 1 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1e = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "60",
 linewidth = 1.5, linestyle = :dot)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "95",
  linewidth = 1.5, linestyle = :dash)
 
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "51")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "44", 
        titlefont = font(10),  
        title = "e) Quare-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Quare NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1f = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "83",
 linewidth = 1.5, linestyle = :dot)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "132",
  linewidth = 1.5, linestyle = :dash)
 
p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "63")

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "48", 
        titlefont = font(10),  
        title = "f) Quare-NK", titleloc = :left
)


## Quare 2 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1g = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "16",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "10",
  linewidth = 1.5, linestyle = :dash)
 


p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "16")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "10", 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")
xlabel!("Size (mm)")


# Quare NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1h = histogram(p.SL1_mm, bins = 5, alpha=1, color = :white, label = "39",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "28",
 linewidth = 1.5, linestyle = :dash)


p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 10, alpha=0.5, color = :orange, label = "39")

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.7, color = :gray, label = "28", 
        titlefont = font(10),  
        title = "h) Quare 2-NK", titleloc = :left
)

xlabel!("Size (mm)")


## put everythnig together

plot(p1a, p1c, p1b, p1d, p1e, p1f, p1g, p1h, bins = 10, layout = (4,2), size = (700, 800))
xlims!((5,30))
savefig("plots/Figure_S3.png")


plot(p1a, p1c, p1b, p1d, p1e, p1f, p1g, p1h, bins = 10, layout = (4,2), size = (700, 800))
ylims!((0,60))
savefig("plots/Figure_S3b.png")



### Killifish size
# Caigual
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1ak = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "Before (12)",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.2, color = :red, label = "End (10)",
  linewidth = 1.5, linestyle = :dash)
 

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 20, alpha=0.5, color = :orange, label = "Introduced (11)")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 20, alpha=0.7, color = :gray, label = "Recovered (6)", 
        titlefont = font(10),  
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")
ylims!((0,6))


# Caigual NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Kz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
p1bk = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "18",
 linewidth = 1.5, linestyle = :dot)

 length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.2, color = :red, label = "18",
 linewidth = 1.5, linestyle = :dash)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "17")

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "10", 
        titlefont = font(10),  
        title = "b) Caigual-NK", titleloc = :left
)
ylims!((0,8))


# Naranjo KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1ck = histogram(p.SL1_mm, bins = 15, alpha=1, color = :white, label = "22",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 15, alpha=0.2, color = :red, label = "26",
  linewidth = 1.5, linestyle = :dash)
 

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 27, alpha=0.5, color = :orange, label = "22")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 27, alpha=0.5, color = :gray, label = "14", 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Naranjo NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1dk = histogram(p.SL1_mm, bins = 15, alpha=1, color = :white, label = "28",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 15, alpha=0.2, color = :red, label = "36",
  linewidth = 1.5, linestyle = :dash)
 

 p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "28")

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "14", 
        titlefont = font(10),  
        title = "d) Naranjo-NK", titleloc = :left
)



## Quare 1 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1ek = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "108",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.2, color = :red, label = "174",
  linewidth = 1.5, linestyle = :dash)
 

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)
length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "99")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "62", 
        titlefont = font(10),  
        title = "e) Quare-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Quare NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1fk = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "66",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 10, alpha=0.2, color = :red, label = "108",
  linewidth = 1.5, linestyle = :dash)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))

histogram!(p.SL1_mm, bins = 30, alpha=0.5, color = :orange, label = "57")

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 30, alpha=0.7, color = :gray, label = "33", 
        titlefont = font(10),  
        title = "f) Quare-NK", titleloc = :left
)


## Quare 2 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1gk = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "49",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "30",
  linewidth = 1.5, linestyle = :dash)
 

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 15, alpha=0.5, color = :orange, label = "49")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)


length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 15, alpha=0.7, color = :gray, label = "30", 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")
xlabel!("Size (mm)")


# Quare NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
p1hk = histogram(p.SL1_mm, bins = 10, alpha=1, color = :white, label = "25",
 linewidth = 1.5, linestyle = :dot)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 5, alpha=0.2, color = :red, label = "17",
  linewidth = 1.5, linestyle = :dash)
 

 p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)

length(p.SL1_mm) - length(findall(isnan.(p.SL1_mm)))
histogram!(p.SL1_mm, bins = 20, alpha=0.5, color = :orange, label = "25")

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)

length(p.SL2_mm) - length(findall(isnan.(p.SL2_mm)))
histogram!(p.SL2_mm, bins = 15, alpha=0.7, color = :gray, label = "17", 
        titlefont = font(10),  
        title = "h) Quare 2-NK", titleloc = :left
)

xlabel!("Size (mm)")


## put everythnig together



plot(p1ak, p1ck, p1bk, p1dk, p1ek, p1fk, p1gk, p1hk, bins = 10, layout = (4,2), size = (700, 800))
xlims!((5,100))
savefig("plots/Figure_S4.png")


plot(p1ak, p1ck, p1bk, p1dk, p1ek, p1fk, p1gk, p1hk, bins = 5, layout = (4,2), size = (700, 800))
ylims!((0,60))
xlims!((5,100))
savefig("plots/Figure_S4b.png")
