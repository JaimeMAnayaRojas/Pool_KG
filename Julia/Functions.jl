## Funcitons for analyses


function HDI(samples; credible_mass=0.95)
	# Computes highest density interval from a sample of representative values,
	# estimated as the shortest credible interval
	# Takes Arguments posterior_samples (samples from posterior) and credible mass (normally .95)
	# Originally from https://stackoverflow.com/questions/22284502/highest-posterior-density-region-and-central-credible-region
	# Adapted to Julialang
	sorted_points = sort(samples)
	ciIdxInc = Int(ceil(credible_mass * length(sorted_points)))
	nCIs = length(sorted_points) - ciIdxInc
	ciWidth = repeat([0.0],nCIs)
	for i in range(1, stop=nCIs)
		ciWidth[i] = sorted_points[i + ciIdxInc] - sorted_points[i]
	end
	HDImin = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)]
	HDImax = sorted_points[findfirst(isequal(minimum(ciWidth)),ciWidth)+ciIdxInc]
	return([HDImin, HDImax])
end


function LOS(v, b = 0)
	return 100*length(findall(v .> b)) ./length(v)
end


function Post_summary(df, digits = 3)
    #dfvm= DataFrame(df, :auto)
    ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), df))
    ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), df))
    ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), df))
    PostP = Vector(((mapcols(x -> LOS(x), df)))[1,:])
    df = DataFrame(Parameters = names(df), median = round.(median.(eachcol(df)), digits=digits),
    l99 = round.(Vector(ci99[1, :]), digits =digits),
    l95 = round.(Vector(ci95[1, :]), digits =digits), 
    l68 = round.(Vector(ci68[1, :]), digits =digits),
    u68 = round.(Vector(ci68[2, :]), digits =digits),
    u95 = round.(Vector(ci95[2, :]), digits =digits),
    u99 = round.(Vector(ci99[2, :]), digits =digits),
    PP = round.(PostP, digits =digits))
   

    return df
end


function my_summary(df, digits = 3)
    vm= DataFrame(df, :auto)
    ci68 = (mapcols(x -> HDI(x, credible_mass=0.68), vm))
    ci95 = (mapcols(x -> HDI(x, credible_mass=0.95), vm))
    ci99 = (mapcols(x -> HDI(x, credible_mass=0.997), vm))
    df = DataFrame(size = z1, median = round.(median.(eachcol(vm)), digits=digits),
    l99 = round.(Vector(ci99[1, :]), digits =digits),
    l95 = round.(Vector(ci95[1, :]), digits =digits), 
    l68 = round.(Vector(ci68[1, :]), digits =digits),
    u68 = round.(Vector(ci68[2, :]), digits =digits),
    u95 = round.(Vector(ci95[2, :]), digits =digits),
    u99 = round.(Vector(ci99[2, :]), digits =digits))

    return df
end


function My_Logit(z)
    l = log(z / (1 - (z)))
    return l 
end 




function My_Logistic(x)
    l = 1/(1+exp(-x))
    return l 
end 
  


### Plot recruitment function


# # size clases
# size_cen = 50.0
# z =  collect(5:1:300)
# z1 = z

# KG = zeros(size(post_mG)[1], length(z));
# NK = zeros(size(post_mG)[1], length(z));



function Recr_link(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, NK::Integer, row::Integer)
    zc = z .- size_cen
    α= df.b_Intercept[row]
    βz= df.b_z[row]
    βNK= df.b_Removal1[row]
    βzNK= df."b_z:Removal1"[row]
   # βz2NK= df."b_z2.NK"[row]
    μ = α .+ βNK .* NK .+  (βz .+ βzNK .* NK) .* zc # linear predictor
    return(exp.(μ))
end



# for i in 1:size(post_mG)[1]
#     KG[i,:] =  (Recr_link(post_mG, z, size_cen, 0, i))
#     NK[i,:] =  (Recr_link(post_mG, z, size_cen, 1, i))
# end

# my_summary(KG)

# # Plot NG
# p = My_Logistic.(my_summary(KG))
# Fig_1A =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
#     linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
# )

# p = My_Logistic.(my_summary(NK))
# plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NK")
# #xlabel!("Initial size (mm)")
# ylabel!("Survival")
# ylims!((-0.01,1.01))
# xlims!(8,30)
# plot!([18,18], [0,100], linestyle = :dash, colour = :gray, label= :false)

    # estimation graphs




# scatter!(DataG.SL1_mm, DataG.surv , groups = DataG.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false, markersize = 3)


# lab= round(LOS(post_SurvG.b_NK), digits = 2)


# Fig_1A_α = plot(kde(post_SurvG.b_Intercept), fillrange=-0, fillalpha=0.25, legend= :false, 
# title = "            Statistical Test  \n   probability that NK > KG:  \n \n $(lab)%", titleloc = :left, titlefontsize = 9)
# plot!(kde(post_SurvG.b_Intercept .+ post_SurvG.b_NK), fillrange=-0, fillalpha=0.25, legend= :false, ticks =:false)
# xlabel!("Intercept")
# plot!([0,0], [0,1.1], c=:red, linewidth=2)


# lab= round(LOS(post_SurvG."b_z.NK"), digits = 2)
# Fig_1A_β = plot(kde(post_SurvG.b_z), fillrange=-0, fillalpha=0.25, legend= :false, 
# title= "\n  $(lab)%", titlefontsize = 9, titleloc= :left, ticks =:false)
# plot!(kde(post_SurvG.b_z .+ post_SurvG."b_z.NK"), fillrange=-0, fillalpha=0.25, legend= :false,  ticks =:false)
# xlabel!("Slope (z)")
# plot!([0,0], [0,12], c=:red, linewidth=2)




# l = @layout [
#     a [b{0.4h}
#        c{0.4h}]
# ]

# FA = plot(Fig_1A, Fig_1A_α, Fig_1A_β,
#     layout = l
# )
