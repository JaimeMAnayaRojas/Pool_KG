v = IPMs[1][1][:,1]


# Scripts for plots
#summ_tab = DataFrame()
df = IPMs[1][1]
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), df))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), df))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), df))

PP0 = (mapcols(x -> LOS(x,0), df))
PP0 =Vector(PP0[1,:])

PP1 = (mapcols(x -> LOS(x,1), df))
PP1 =Vector(PP1[1,:])

# p_val = (mapcols(x -> boot_p(x), Δ13C_net))
G = DataFrame(parameter = names(df), median = round.(median.(eachcol(df)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3),
pp0= PP0, pp1= PP1 
);


#append!(summ_tab, summ_tab1)

df = IPMs[2][1]
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), df))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), df))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), df))
PP0 = (mapcols(x -> LOS(x,0), df))
PP0 =Vector(PP0[1,:])

PP1 = (mapcols(x -> LOS(x,1), df))
PP1 =Vector(PP1[1,:])

# p_val = (mapcols(x -> boot_p(x), Δ13C_net))
K = DataFrame(parameter = names(df), median = round.(median.(eachcol(df)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3),
pp0= PP0, pp1= PP1 
);


sTab = [G, K]

## Make plots


G_Tab = DataFrame(Parameter = sTab[1][:,:parameter],
          median = sTab[1][:,:median], 
		  l95 = sTab[1][:,:l95], l68 = sTab[1][:,:l68],
          u68 = sTab[1][:,:u68], u95 = sTab[1][:,:u95],  
		  Species = "Guppy"
)


K_Tab = DataFrame(Parameter = sTab[2][:,:parameter],
          median = sTab[2][:,:median], 
		  l95 = sTab[2][:,:l95], l68 = sTab[2][:,:l68],
          u68 = sTab[2][:,:u68], u95 = sTab[2][:,:u95], 
		  Species = "Killifish"
)

KG_tab = vcat(G_Tab, K_Tab)

KG_tab = KG_tab[[1,2,10,11],:]

x = [0, 4.5 ]
y = [1,1]

KG_tab.X = [1,2,3,4]

pA = scatter(KG_tab.X, KG_tab.median, #yerr=[KG_tab.median - KG_tab.lc KG_tab.median - KG_tab.up],
		 markersize = 8,
        group = KG_tab.Species, label = false, markercolor = :white, m = [:circle :rect],
        xlims = (0.5,4.5), title = "a)",  titleloc = :left, titlefont = font(10) 
)
	plot!([1 ,2, NaN, 3,4  ], [KG_tab.median[1], KG_tab.median[2], NaN, KG_tab.median[3], KG_tab.median[4]], c = :gray, linewidth = 2, lab = false)
	plot!([0.5, 4.5], [1.0,1.0], c = :gray, linestyle = :dash, lab = false, linewidth = 1.5)
	plot!([2.5, 2.5], [.5,3], c = :gray,  lab = false, linewidth = 2)
	
	plot!([1, 1, NaN, 2 , 2,  NaN , 3, 3, NaN, 4,4], 
		  [KG_tab[1, :l95], KG_tab[1, :u95], NaN,
		  KG_tab[2, :l95], KG_tab[2, :u95], NaN,
		  KG_tab[3, :l95], KG_tab[3, :u95], NaN,
		  KG_tab[4, :l95], KG_tab[4, :u95]],
		   c = :black,  lab = false, linewidth = 1.5)
  	plot!([1, 1, NaN, 2 , 2,  NaN , 3, 3, NaN, 4,4], 
		  [KG_tab[1, :l68], KG_tab[1, :u68], NaN,
		  KG_tab[2, :l68], KG_tab[2, :u68], NaN,
		  KG_tab[3, :l68], KG_tab[3, :u68], NaN,
		  KG_tab[4, :l68], KG_tab[4, :u68]],
		   c = :black,  lab = false, linewidth = 4)
	scatter!(KG_tab.X, KG_tab.median, #yerr=[KG_tab.median - KG_tab.lc KG_tab.median - KG_tab.up],
		   markersize = 8,
		  group = KG_tab.Species, label = ["Guppy" "Killifish"], legend =:topleft, markercolor = :gray, m = [:circle :rect]
  	)
	

	ylims!((0.5,2.8))
	ylabel!("Fitness (λ)")
	xticks!([1,2,3,4], ["KG", "NK", "KG", "NG"], fontsize = 11)
	xlabel!("Pool treatment",   fontsize = 10)

# LTRE
## Make plots

G_Tab = DataFrame(Parameter = sTab[1][:,:parameter],
          median = sTab[1][:,:median], 
		  l95 = sTab[1][:,:l95], l68 = sTab[1][:,:l68],
          u68 = sTab[1][:,:u68], u95 = sTab[1][:,:u95],  
		  Species = "Guppy"
)


K_Tab = DataFrame(Parameter = sTab[2][:,:parameter],
          median = sTab[2][:,:median], 
		  l95 = sTab[2][:,:l95], l68 = sTab[2][:,:l68],
          u68 = sTab[2][:,:u68], u95 = sTab[2][:,:u95], 
		  Species = "Killifish"
)

KG_tab = vcat(G_Tab, K_Tab)

KG_tab = KG_tab[[4,8,5,6,7,13,17,14,15,16], :]  
KG_tab.X = ["A", "B", "C", "D", "E","A", "B", "C", "D", "E"]
df = filter(:Species => x -> x == "Guppy", KG_tab)

pB = bar(df.X, df.median, label = false, group = df.X, color = :lightgray,# [:white :red :deepskyblue2 :blue :black], 
title = "b)", titleloc=:left, titlefont = 10)
#ylabel!("Fitness contribution (Δλ)")
#xlabel!("Vital rate")
xticks!([0.5, 1.95, 3.4, 4.8, 6.25], ["All", "S(z)", "G(z, z')", "R(z)", "O(z)"], tickfontrotation = 0.5)
plot!([0, 6], [0,0], c = :black,  lab = false, linewidth = 1)

plot!([0.5, 0.5, NaN, 1.95, 1.95, NaN, 3.4, 3.4, NaN, 4.8, 4.8, NaN, 6.25, 6.25], 
[df[1, :l95], df[1, :u95], NaN,
df[2, :l95], df[2, :u95], NaN,
df[3, :l95], df[3, :u95], NaN,
df[4, :l95], df[4, :u95], NaN,
df[5, :l95], df[5, :u95]],
 c = :black,  lab = false, linewidth = 1.5)

plot!([0.5, 0.5, NaN, 1.95, 1.95, NaN, 3.4, 3.4, NaN, 4.8, 4.8, NaN, 6.25, 6.25], 
[df[1, :l68], df[1, :u68], NaN,
df[2, :l68], df[2, :u68], NaN,
df[3, :l68], df[3, :u68], NaN,
df[4, :l68], df[4, :u68], NaN,
df[5, :l68], df[5, :u68]],
 c = :black,  lab = false, linewidth = 4)
 ylabel!("∑ Δλ")


 df = filter(:Species => x -> x == "Killifish", KG_tab)
 pC = bar(df.X, df.median, label = false, group = df.X, color = :lightgray, # [:white :red :deepskyblue2 :blue :black], 
 	title = "d)", titleloc=:left, titlefont = 10)
 ylabel!("∑Δλ")
 xlabel!("Vital rate")
 xticks!([0.5, 1.95, 3.4, 4.8, 6.25], ["All", "S(z)", "G(z, z')", "R(z)", "O(z)"])
 plot!([0, 6], [0,0], c = :black,  lab = false, linewidth = 1)
 
 plot!([0.5, 0.5, NaN, 1.95, 1.95, NaN, 3.4, 3.4, NaN, 4.8, 4.8, NaN, 6.25, 6.25], 
 [df[1, :l95], df[1, :u95], NaN,
 df[2, :l95], df[2, :u95], NaN,
 df[3, :l95], df[3, :u95], NaN,
 df[4, :l95], df[4, :u95], NaN,
 df[5, :l95], df[5, :u95]],
  c = :black,  lab = false, linewidth = 1.5)
 
 plot!([0.5, 0.5, NaN, 1.95, 1.95, NaN, 3.4, 3.4, NaN, 4.8, 4.8, NaN, 6.25, 6.25], 
 [df[1, :l68], df[1, :u68], NaN,
 df[2, :l68], df[2, :u68], NaN,
 df[3, :l68], df[3, :u68], NaN,
 df[4, :l68], df[4, :u68], NaN,
 df[5, :l68], df[5, :u68]],
  c = :black,  lab = false, linewidth = 4)
 
##### Size contributions

## Guppy# size-speficif LTRE
	# nBigMatrix = 100
	min_size = 4
	max_size = 35
	size_cen = 18.0

	U= Float64(max_size)
	L=Float64(min_size)
	m = nBigMatrix
	h = (U - L)/m
	z1 =  L .+ (collect(1:m) .- 0.5) * h

vm = DataFrame(IPMs[1][2], :auto)

ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))

df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3)
);


println(df)

pD = plot(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.1, label = "S(z)", linewidth = 2)
plot!([min_size, max_size], [0,0], c = :black,  lab = false, linewidth = 0.5, linestyle = :dash)

vm= DataFrame(IPMs[1][3], :auto)
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3));

plot!(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.2, label = "G(z,z')", linewidth = 2)


vm= DataFrame(IPMs[1][4], :auto)
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3));
plot!(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.1, 
    label = "R(z)", legend= :bottomright, linewidth = 2, 
    title = "c)", 
    titleloc= :left, titilefont = 9
)
xlabel!("Guppy size (mm)")
ylabel!("Δλ")


# Killifish
# size-speficif LTRE
# nBigMatrix = 100
min_size = 2
max_size = 110
size_cen = 18.0

U= Float64(max_size)
L=Float64(min_size)
m = nBigMatrix
h = (U - L)/m
z1 =  L .+ (collect(1:m) .- 0.5) * h
z = 0
z1 = vcat(z,z1)


vm = DataFrame(IPMs[2][2], :auto)


ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))

df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3)
);

println(df)

pE = plot(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.1, label = false, linewidth = 2)
plot!([min_size, max_size], [0,0], c = :black,  lab = false, linewidth = 0.5, linestyle = :dash)

vm= DataFrame(IPMs[2][3], :auto)
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3));
println(df)
plot!(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.2, label = false, linewidth = 2)


vm= DataFrame(IPMs[2][4], :auto)
ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=3),
l99 = round.(Vector(ci99[1, :]), digits =3),
l95 = round.(Vector(ci95[1, :]), digits =3), 
l68 = round.(Vector(ci68[1, :]), digits =3),
u68 = round.(Vector(ci68[2, :]), digits =3),
u95 = round.(Vector(ci95[2, :]), digits =3),
u99 = round.(Vector(ci99[2, :]), digits =3));
plot!(df[:,:size], df[:,:median], ribbon = [df[:, :l68] df[:, :u68] ], fillalpha = 0.2, 
label = false, linewidth = 2, 
title = "e)", 
titleloc= :left, titilefont = 9
)
xlabel!("Killifish size (mm)")
ylabel!("Δλ")





#  
l = @layout [
    a 
	[grid(2,2)]
]

using Plots.Measures
plot(pA, pB, pD, pC, pE,
    layout = l, size = (600, 800), titlefont = 11, foreground_color_legend = nothing, 
    margin = [2mm 3mm]
)


