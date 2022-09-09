# size clases
size_cen = 18.0
z =  collect(5:0.1:30)
z1 = z
# survival
KG = zeros(size(post_SurvG)[1], length(z));
NK = zeros(size(post_SurvG)[1], length(z));



function p_link(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NK::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βz2= df.b_z2[row]
    βNK= df."b_NK"[row]
    βzNK= df."b_z.NK"[row]
   # βz2NK= df."b_z2.NK"[row]
    μ = α .+ βNK .* NK .+  (βz .+ βzNK .* NK) .* zc .+ βz2 .* z2  # linear predictor
    return(μ)
end




function sG_link(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NK::Integer, row::Integer)
    zc = z .- size_cen
    z2 = z.^2 .- size_cen.^2
    α= df.b_Intercept[row]
    βz= df.b_z[row]
   
    βNK= df."b_NK"[row]
    βzNK= df."b_z.NK"[row]
   
    μ = α .+ βNK .* NK .+  (βz .+ βzNK .* NK) .* zc # linear predictor
    return(μ)
end


for i in 1:size(post_SurvG)[1]
    KG[i,:] =  (sG_link(post_SurvG, z, size_cen, 0, i))
    NK[i,:] =  (sG_link(post_SurvG, z, size_cen, 1, i))
end

my_summary(KG)

# Plot NG
p = My_Logistic.(my_summary(KG))
Fig_1A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
)

p = My_Logistic.(my_summary(NK))
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")
ylabel!("Survival")
ylims!((-0.01,1.01))
xlims!(8,30)
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, label= :false)

    # estimation graphs




scatter!(DataG.SL1_mm, DataG.surv , groups = DataG.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false, markersize = 3)


lab= round(LOS(post_SurvG.b_NK), digits = 2)


Fig_1A_α = plot(kde(post_SurvG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "            Statistical Test  \n   probability that NK > KG:  \n \n $(lab)%", titleloc = :left, titlefontsize = 9)
plot!(kde(post_SurvG.b_Intercept .+ post_SurvG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")
plot!([0,0], [0,1.1], c=:red, linewidth=2)


lab= round(LOS(post_SurvG."b_z.NK"), digits = 2)
Fig_1A_β = plot(kde(post_SurvG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_SurvG.b_z .+ post_SurvG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope (z)")
plot!([0,0], [0,12], c=:red, linewidth=2)




l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FA = plot(Fig_1A, Fig_1A_α, Fig_1A_β,
    layout = l
)

# growth
KG = zeros(size(post_GrowthG)[1], length(z));
NK = zeros(size(post_GrowthG)[1], length(z));



for i in 1:size(post_GrowthG)[1]
    KG[i,:] =  (p_link(post_GrowthG, z, size_cen, 0, i))
    NK[i,:] =  (p_link(post_GrowthG, z, size_cen, 1, i))
end


# Plot NG
p = my_summary(KG)
Fig_2A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "b)", titleloc = :left, legend = :bottomright    
)

p = my_summary(NK)
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")


plot!([18,18], [-0.2,maximum(p.u68)], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Growth ln(z₁/z)")
#ylims!((5,30))
xlims!(8,30)


df = filter(:Sex2 => x -> x != "M", subset(DataG, :Sex2 .=> ByRow(!ismissing)))

df.growth = log.(df.SL2_mm ./ df.SL1_mm)
scatter!(df.SL1_mm, df.growth, groups = df.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false, markersize = 3)

    # estimation graphs



lab= round(LOS(post_GrowthG.b_NK), digits = 2)

Fig_2A_α = plot(kde(post_GrowthG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_GrowthG.b_Intercept .+ post_GrowthG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")
plot!([0,0], [0,10], c=:red, linewidth=2)

lab= round(LOS(post_GrowthG."b_z.NK"), digits = 2)
Fig_2A_β = plot(kde(post_GrowthG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_GrowthG.b_z .+ post_GrowthG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope (z)")
plot!([0,0], [0,45], c=:red, linewidth=2)



l = @layout [
    a [b{0.4h}
       c{0.4h}]
]

FB = plot(Fig_2A, Fig_2A_α, Fig_2A_β,
    layout = l
)



# Reproduction


# function Rep_link(df::AbstractDataFrame, z::AbstractVector, 
#     size_cen::AbstractFloat, NK::Integer, row::Integer)
#     zc = z .- size_cen
#    # z2 = z.^2 .- size_cen.^2
#     α= df.b_Intercept[row]
#     βz= df.b_z[row]
#     #βz2= df.b_z2[row]
#     βNK= df."b_NK"[row]
#     βzNK= df."b_z.NK"[row]
#     μ = α .+ βNK .* NK .+  (βz .+ βzNK .* NK) .* zc  # linear predictor
#     return(μ)
# end


KG = zeros(size(post_ReprG)[1], length(z));
NK = zeros(size(post_ReprG)[1], length(z));



for i in 1:size(post_ReprG)[1]
    KG[i,:] =  (p_link(post_ReprG, z, size_cen, 0, i))
    NK[i,:] =  (p_link(post_ReprG, z, size_cen, 1, i))
end


# Plot NG
p = exp.(my_summary(KG))
Fig_3A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "c)", titleloc = :left, legend = :bottomright    
)

p = exp.(my_summary(NK))
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")
plot!([18,18], [0,100], linestyle = :dash, colour = :gray, legend = :false)

ylabel!("Offspring (N)")

scatter!(df.SL1_mm, df.Recr, groups = df.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false, markersize = 3)


ylims!((-1,20))

xlims!(8,30)

xlabel!("Guppy size (mm)")

    # estimation graphs



lab= round(LOS(post_ReprG.b_NK), digits = 2)

Fig_3A_α = plot(kde(post_ReprG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
title = "\n $(lab)%", titlefontsize = 9, titleloc = :left)
plot!(kde(post_ReprG.b_Intercept .+ post_ReprG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
xlabel!("Intercept")

plot!([0,0], [0,1.2], c=:red, linewidth=2)

lab= round(LOS(post_ReprG."b_z.NK"), digits = 2)
Fig_3A_β = plot(kde(post_ReprG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
plot!(kde(post_ReprG.b_z .+ post_ReprG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
xlabel!("Slope (z)")
plot!([0,0], [0,2.2], c=:red, linewidth=2)


l = @layout [
    a [b{0.4h}
       c{0.4h}]
]


FC = plot(Fig_3A, Fig_3A_α, Fig_3A_β,
    layout = l
)


plot(FA, FB, FC, layout = (3,1), size = (500, 700), leftmargin = 5mm, rightmargin = 5mm, grid = :false)


