
function Surv_linkK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NG::Integer, row::Integer)
    zc = z .- size_cen
   # z2 = z.^2 .- size_cen.^2
   # βz2= df.b_z2[row]
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βNG= df."b_NG"[row]
    βzNG= df."b_z.NG"[row]
    μ = α .+ βNG .* NG .+  (βz .+ βzNG .* NG) .* zc #.+ βz2 .* z2 # linear predictor
    return(μ)
end


size_cen = 18.0

z = collect(5:0.1:95);
z1 = z;


# survival
KG = zeros(size(post_SurvK)[1], length(z));
NG = zeros(size(post_SurvK)[1], length(z));

for j in 1:size(post_SurvK)[1]
    KG[j,:] =  (Surv_linkK(post_SurvK, z, size_cen, 0, j))
    NG[j,:] =  (Surv_linkK(post_SurvK, z, size_cen, 1, j))
end

my_summary(KG)

# Plot NG
p = My_Logistic.(my_summary(KG)).*100
Fig_1A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
)

p = My_Logistic.(my_summary(NG)).*100
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NG")
#xlabel!("Initial size (mm)")
ylabel!("Survival (%)")
ylims!((0,100))
xlims!(5,95)
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, label= :false)

# estimation graphs

lab= round(LOS(post_SurvK.b_NG), digits = 2)


Fig_1A_α = plot(kde(post_SurvK.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "            Statistical test \n                KG > NG  \n $(lab)%", titleloc = :left, titlefontsize = 9)
plot!(kde(post_SurvK.b_Intercept .+ post_SurvK.b_NG), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

plot!([0,0], [0,1.5], c=:red, linewidth=2)





lab= round(LOS(post_SurvK."b_z.NG"), digits = 2)
Fig_1A_β = plot(kde(post_SurvK.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_SurvK.b_z .+ post_SurvK."b_z.NG"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope")

plot!([0,0], [0,40], c=:red, linewidth=2)



l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FA = plot(Fig_1A, Fig_1A_α, Fig_1A_β,
    layout = l
)

# growth


KG = zeros(size(post_GrowthK)[1], length(z));
NG = zeros(size(post_GrowthK)[1], length(z));


for i in 1:size(post_GrowthK)[1]
    KG[i,:] =  (G_linkK(post_GrowthK, z, size_cen, 0, i))
    NG[i,:] =  (G_linkK(post_GrowthK, z, size_cen, 1, i))
end

# Plot NG
p = my_summary(KG)
Fig_2A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "b)", titleloc = :left, legend = :bottomright    
)

p = my_summary(NG)
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NG")
#xlabel!("Initial size (mm)")
α = mean(post_GrowthK.b_omega_Intercept)

plot!([α,α], [-0.12,0.7], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Growth ln(z₁/z)")
xlims!(5,95)
#ylims!(-0.3,0.3)


DataK.growth = log.(DataK.SL2_mm ./ DataK.SL1_mm) 
scatter!(DataK.SL1_mm, DataK.growth, groups = DataK.NG, c= [:lightskyblue, :red], alpha = 0.8, label =false)


xlabel!("Killifish size (mm)")
    # estimation graphs


names(post_GrowthK)

lab= round(LOS(post_GrowthK.b_b0_NG), digits = 2)

Fig_2A_α = plot(kde(post_GrowthK.b_b0_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_GrowthK.b_b0_Intercept .+ post_GrowthK.b_b0_NG), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

lab= round(LOS(post_GrowthK.b_bz1_NG), digits=2)
Fig_2A_β = plot(kde(post_GrowthK.b_bz1_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_GrowthK.b_bz1_Intercept .+ post_GrowthK.b_bz1_NG), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)

α = round(α, digits = 2)
xlabel!("Slope 1")


lab= round(LOS(post_GrowthK.b_bz2_NG), digits=2)
Fig_2A_c = plot(kde(post_GrowthK.b_bz2_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n $(lab)%", titlefontsize = 9, titleloc = :left, ticks =:false)
plot!(kde(post_GrowthK.b_bz2_Intercept .+ post_GrowthK.b_bz2_NG), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)

α = round(α, digits = 2)
xlabel!("Slope 2")


l = @layout [
    a [b{0.4h}
       grid(1,2){0.4h}]
]

FB = plot(Fig_2A, Fig_2A_α, Fig_2A_β, Fig_2A_c,
    layout = l
)



# Reproduction

# estimation graphs
function Rep_linkK(df::AbstractDataFrame, z::AbstractVector, 
      size_cen::AbstractFloat, NG::Integer, row::Integer)
    zc = z .- size_cen
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    μ = α .+ βz .* zc  # linear predictor
    return(μ)
end
    

KG = zeros(size(post_ReprK)[1], length(z));
NG = zeros(size(post_ReprK)[1], length(z));



for i in 1:size(post_ReprK)[1]
    KG[i,:] =  (Rep_linkK(post_ReprK, z, size_cen, 0, i))
    NG[i,:] =  (Rep_linkK(post_ReprK, z, size_cen, 1, i))
end


# Plot NG
p = exp.(my_summary(KG))
Fig_3A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "c)", titleloc = :left, legend = :bottomright, c = :gray    
)

# p = exp.(my_summary(NG))
# plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NG")
# #xlabel!("Initial size (mm)")
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, legend = :false)
ylabel!("Offspring (N)")

scatter!(DataK.SL1_mm, DataK.Recr, groups = DataK.NG, c= [:lightskyblue, :red], alpha = 0.8, label =false)
ylims!((-1,20))
xlims!(5,90)



Fig_3A_α = plot(kde(post_ReprK.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, c = :gray,   ticks =:false)
xlabel!("Intercept")
plot!([0,0], [0,0.4], c=:red, linewidth=2)


Fig_3A_β = plot(kde(post_ReprK.b_z), fillrange=-0, fillalpha=0.25, legend= :false, c =:gray,  ticks =:false)
xlabel!("Slope")
plot!([0,0], [0,15], c=:red, linewidth=2)


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]


FC = plot(Fig_3A, Fig_3A_α, Fig_3A_β,
    layout = l
)


plot(FA, FB, FC, layout = (3,1), size = (500, 700))



