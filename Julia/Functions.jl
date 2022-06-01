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
  
