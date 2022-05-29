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


filter!(:Location => x -> x !="JOE", Gz)
filter!(:Location => x -> x !="JOE", Kz)

println(describe(Gz))

replace!(Gz.Mark, missing=>"NaN")
replace!(Kz.Mark, missing=>"NaN")

filter!(:Mark => x -> x !="others", Gz)
filter!(:Mark => x -> x !="others", Kz)




Gz.Mark = categorical(Gz.Mark)
Kz.Mark = categorical(Kz.Mark)


replace!(Gz.mass1_gr, missing=>NaN)
replace!(Gz.mass2_gr, missing=>NaN)


replace!(Kz.mass1_gr, missing=>NaN)
replace!(Kz.mass2_gr, missing=>NaN)

base = 1.5

Gz.mass1_gr = log.(base, Gz.mass1_gr)
Gz.mass2_gr = log.(base, Gz.mass2_gr)

Kz.mass1_gr = log.(base, Kz.mass1_gr)
Kz.mass2_gr = log.(base, Kz.mass2_gr)



Gz.Location = categorical(Gz.Location)



# Caigual
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1a = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = "Removed",
 linewidth = 1.5, linestyle = :dot)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.1, color = :red, label = "Recovered",
linewidth = 1.5, linestyle = :dash)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = "Introduced")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = "Recovered", legend =:topleft,
        titlefont = font(10),  
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Caigual NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1b = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "b) Caigual-NK", titleloc = :left
)

# Naranjo KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1c = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Naranjo NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1d = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "d) Naranjo-NK", titleloc = :left
)

## Quare 1 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1e = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)
length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "e) Quare-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Quare NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1f = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Gz
)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "f) Quare-NK", titleloc = :left
)


## Quare 2 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1g = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)
length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")
xlabel!("Log₁₀(Body mass, gr)")

# Quare NK

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1h = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NK, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Gz
)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "h) Quare 2-NK", titleloc = :left
)

xlabel!("Log₁₀(Body mass, gr)")

       

## put everythnig together

plot(p1a, p1b, p1c, p1d, p1e, p1f, p1g, p1h, bins = 10, layout = (4,2), size = (700, 800))
xlims!((-2.5,0))
savefig("plots/Figure_S5.png")


plot(p1a, p1b, p1c, p1d, p1e, p1f, p1g, p1h, bins = 10, layout = (4,2), size = (700, 800))
ylims!((0,21))
xlims!((-2.5,0))
savefig("plots/Figure_S5b.png")



### Killifish size
# Caigual
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1ak = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = "Before",
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = "Introduced")

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = "Recovered", 
        titlefont = font(10),  legend = :topleft,
        title = "a) Caigual-KG", titleloc = :left
)
ylabel!("Frequency (N)")
ylims!((0,6))


# Caigual NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1bk = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "CAI" && w != "NaN", Kz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "b) Caigual-NK", titleloc = :left
)
ylims!((0,6))


# Naranjo KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1ck = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "c) Naranjo-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Naranjo NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1dk = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "NAR" && w != "NaN", Kz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "d) Naranjo-NK", titleloc = :left
)

## Quare 1 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1ek = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)
length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "e) Quare-KG", titleloc = :left
)
ylabel!("Frequency (N)")


# Quare NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1fk = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA" && w != "NaN", Kz
)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "f) Quare-NK", titleloc = :left
)


## Quare 2 KG
p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1gk = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)
length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:KG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)


length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "g) Quare 2-KG", titleloc = :left
)
ylabel!("Frequency (N)")
xlabel!("Log₁₀(Body mass, gr)")


# Quare NK

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))
p1hk = histogram(p.mass1_gr, bins = 30, alpha=1, color = :white, label = false,
 linewidth = 1.5, linestyle = :dot)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)

length(p.mass1_gr) - length(findall(isnan.(p.mass1_gr)))

histogram!(p.mass1_gr, bins = 30, alpha=0.5, color = :orange, label = false)

p = filter([:NG, :Location, :Mark] => (x,y,w) -> x == 1 && 
        y == "QUA2" && w != "NaN", Kz
)

length(p.mass2_gr) - length(findall(isnan.(p.mass2_gr)))
histogram!(p.mass2_gr, bins = 30, alpha=0.7, color = :gray, label = false, 
        titlefont = font(10),  
        title = "h) Quare 2-NK", titleloc = :left
)

xlabel!("Log₁₀(Body mass, gr)")



## put everythnig together



plot(p1ak, p1bk, p1ck, p1dk, p1ek, p1fk, p1gk, p1hk, bins = 10, layout = (4,2), size = (700, 800))
xlims!((-2,1))
savefig("plots/Figure_S6.png")


plot(p1ak, p1bk, p1ck, p1dk, p1ek, p1fk, p1gk, p1hk, bins = 5, layout = (4,2), size = (700, 800))
ylims!((0,13))
xlims!((-2,1))
savefig("plots/Figure_S6b.png")
