
post = CSV.read("outputs/Posteriors.csv", DataFrame)# Statistical analyses

DataG = CSV.read("data/GuppyIMP.csv", DataFrame);
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





# size clases
# Guppy
# nBigMatrix = 100
min_size = 4
max_size = 35
size_cen = 18.0

U= Float64(max_size)
L=Float64(min_size)
m = nBigMatrix
h = (U - L)/m
z1 =  L .+ (collect(1:m) .- 0.5) * h
z = z1
meshpts = z1


##########33


z1k =  2 .+ (collect(1:m) .- 0.5) .* ((110 - 2)/m)

########333


pars_GR_G = select(post,
:Intercept_survG => :"α_surv",
:b_z_survG => :"βz_surv",

:b_canopy_survG => :"β_canopy_surv",

:Intercept_growG => :"α_grow",
:b_z_growG => :"βz_grow",

:b_canopy_growG => :"β_canopy_grow",
:sigma_growG => :σ_grow,

:Intercept_recrG => :"α_fec",
:b_z_recrG => :"βz_fec",

:b_canopy_recrG => :"β_canopy_fec"
)		


pars_NR_G = select(post,
:Intercept_survG => :"α_surv",
:b_z_survG => :"βz_surv",

:b_canopy_survG => :"β_canopy_surv",

:Intercept_growG => :"α_grow",
:b_z_growG => :"βz_grow",

:b_canopy_growG => :"β_canopy_grow",
:sigma_growG => :σ_grow,

:Intercept_recrG => :"α_fec",
:b_z_recrG => :"βz_fec",
:b_BiomassK_recrG => :"β_BiomassK_fec",
:b_canopy_recrG => :"β_canopy_fec"


)

pars_NR_G.α_surv = pars_NR_G.α_surv .+ post.b_NK_survG
pars_NR_G.βz_surv = pars_NR_G.βz_surv .+ post.b_zNK_survG 
pars_NR_G.α_grow = pars_NR_G.α_grow .+ post.b_NK_growG
pars_NR_G.βz_grow = pars_NR_G.βz_grow .+ post.b_zNK_growG 
pars_NR_G.α_fec = pars_NR_G.α_fec .+ post.b_NK_recrG
pars_NR_G.βz_fec = pars_NR_G.βz_fec .+ post.b_zNK_recrG 


## Killifish

println(names(post))
pars_GR_K = select(post,
:Intercept_survK => :"α_surv",
:b_z_survK => :"βz_surv",

:b_canopy_survK => :"β_canopy_surv",

:Intercept_growK => :"α_grow",
:b_z_growK => :"βz_grow",
:b_z2_growK => :"βz2_grow",

:b_canopy_growK => :"β_canopy_grow",
:sigma_growK => ByRow(x-> sqrt(x)) =>:σ_grow,

:Intercept_recrK => :"α_fec",
:b_z_recrK => :"βz_fec",

:b_canopy_recrK => :"β_canopy_fec"


)



pars_NG_K = select(post,
:Intercept_survK => :"α_surv",
:b_z_survK => :"βz_surv",

:b_canopy_survK => :"β_canopy_surv",

:Intercept_growK => :"α_grow",
:b_z_growK => :"βz_grow",
:b_z2_growK => :"βz2_grow",


:b_canopy_growK => :"β_canopy_grow",
:sigma_growK => ByRow(x-> sqrt(x)) =>:σ_grow,

:Intercept_recrK => :"α_fec",
:b_z_recrK => :"βz_fec",

:b_canopy_recrK => :"β_canopy_fec"
)



pars_NG_K.α_surv = pars_NG_K.α_surv .+ post.b_NG_survK
pars_NG_K.βz_surv = pars_NG_K.βz_surv .+ post.b_zNG_survK 
pars_NG_K.α_grow = pars_NG_K.α_grow .+ post.b_NG_growK
pars_NG_K.βz_grow = pars_NG_K.βz_grow .+ post.b_zNG_growK 
pars_NG_K.α_fec = pars_NG_K.α_fec .+ post.b_NG_recrK
pars_NG_K.βz_fec = pars_NG_K.βz_fec .+ post.b_zNG_recrK 

## Surival guppy

function s_z(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    α= df.α_surv[row]
    β= df.βz_surv[row]
    
    linear_p = α .+ β * (z .- size_cen)  # linear predictor
    p = 1 ./(1 .+ exp.(-linear_p))
    p = diagm(p)
    return(p)
end


KG = zeros(size(pars_GR_G)[1], length(z))
NK = zeros(size(pars_NR_G)[1], length(z))
for i in 1:size(pars_GR_G)[1]
    KG[i,:] =  sum(eachcol(s_z(pars_GR_G, z, size_cen, i)))
    NK[i,:] =  sum(eachcol(s_z(pars_NR_G, z, size_cen, i)))
end

p = My_Logistic.(my_summary(My_Logit.(KG))).*100
pS1a =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "a)", titleloc = :left, legend = :bottomright    
    )
p = My_Logistic.(my_summary(My_Logit.(NK))).*100
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NK")
#xlabel!("Initial size (mm)")
ylabel!("Guppy \n Survival (%)")
ylims!((0,100))
xlims!(5,30)

######################################################################

KG = zeros(size(pars_GR_K)[1], length(z))
NG = zeros(size(pars_NG_K)[1], length(z))

for i in 1:size(pars_NG_K)[1]
    KG[i,:] =  sum(eachcol(s_z(pars_GR_K, z1k, size_cen, i)))
    NG[i,:] =  sum(eachcol(s_z(pars_NG_K, z1k, size_cen, i)))
end

p = My_Logistic.(my_summary(My_Logit.(KG))).*100
pS1d =plot(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 5, label = "KG", title = "d)", titleloc = :left, legend = :bottomright    
    )
p = My_Logistic.(my_summary(My_Logit.(NG))).*100
plot!(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), linewidth = 5, label = "NG")
xlabel!("Initial size (mm)")
ylabel!("Killifish \n Survival (%)")
ylims!((0,100))
xlims!(5,100)
#scatter!(z, DataG.surv.*100)
my_summary(NK)



### Growth plot


function g_z1z(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    α= df.α_grow[row]
    β= df.βz_grow[row]
    σ= df.σ_grow[row]

   p_den_grow = zeros(size(z)[1],size(z)[1])
    μ = α .+ β * (z .- size_cen ) 
    for i in 1:nBigMatrix
        p_den_grow[:,i] = pdf.(Normal(μ[i], σ), z1).*h
    end
    return(μ)
end

KG = zeros(size(pars_GR_G)[1], length(z))
NK = zeros(size(pars_NR_G)[1], length(z))


for i in 1:size(pars_NR_G)[1]
    KG[i,:] =  sum(eachcol(g_z1z(pars_GR_G, z, z,size_cen, i)))
    NK[i,:] =  sum(eachcol(g_z1z(pars_NR_G, z, z, size_cen, i)))
end

p = my_summary(KG)
pS1b =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 3, title = "b)", titleloc = :left, label = false  
    )
p = my_summary(NK)
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
linewidth = 3, label = false)

# xlabel!("Guppy initial size \n (mm)")
ylabel!("Final size \n (mm)")

scatter!(DataG.SL1_mm, DataG.SL2_mm, groups = DataG.NK, c= [:lightskyblue, :red], alpha = 0.8, label =false)
xlims!(5,32)
ylims!(5,32)

# Killifish

KG = zeros(size(pars_GR_K)[1], length(z))
NG = zeros(size(pars_NG_K)[1], length(z))


function g_z1zK2(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    α= df.α_grow[row]
    β= df.βz_grow[row]
    β2= df.βz2_grow[row]
    σ= df.σ_grow[row]

   p_den_grow = zeros(size(z)[1],size(z)[1])
    μ = α .+ β * (z .- size_cen ) .+ β2 * (z.^2 .- size_cen^2 ) 
    for i in 1:nBigMatrix
        p_den_grow[:,i] = pdf.(Normal(μ[i], σ), z1).*h
    end
    return(μ)
end


for i in 1:size(pars_NG_K)[1]
    KG[i,:] =  sum(eachcol(g_z1zK2(pars_GR_K, z1k, z1k, size_cen, i)))
    NG[i,:] =  sum(eachcol(g_z1zK2(pars_NG_K, z1k, z1k, size_cen, i)))
end

p = my_summary(KG)

pS1e =plot(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 3, title = "e)", titleloc = :left, label = false 
)



p = my_summary(NG)
plot!(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
linewidth = 3, label = false)

xlabel!("Initial size \n (mm)")
ylabel!("Final size \n (mm)")

scatter!(DataK.SL1_mm, DataK.SL2_mm, groups = DataK.NG, 
c = [:lightskyblue, :red], alpha = 0.8, label = false)

ylims!(5,100)
xlims!(5,100)

## Fecundity
function pr_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    α= df.α_fec[row]
    β= df.βz_fec[row]
    

    linear_p = α .+ β * (z .- size_cen)     # linear predictor
    p = exp.(linear_p)#*(1/2)
    p = diagm(p)
    return(p)
end


KG = zeros(size(pars_GR_G)[1], length(z));
NK = zeros(size(pars_NR_G)[1], length(z));


for i in 1:size(pars_NR_G)[1]
    KG[i,:] =  sum(eachcol(pr_z(pars_GR_G, z, size_cen, i)))
    NK[i,:] =  sum(eachcol(pr_z(pars_NR_G, z,  size_cen, i)))
end

p = my_summary(KG)

pS1c =plot(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 3, label = false, title = "c)", titleloc = :left    
)

p = my_summary(NK)
plot!(z, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
linewidth = 3, label = false)

# xlabel!("Initial size \n (mm)")
ylabel!("Offspring (N)")
scatter!(DataG.SL1_mm, DataG.Recr, groups = DataG.NK, c = [:lightskyblue, :red], alpha = 0.8, label = false)

xlims!(5,30)

# Killifish

KG = zeros(size(pars_GR_K)[1], length(z))
NG = zeros(size(pars_NG_K)[1], length(z))

for i in 1:size(pars_NG_K)[1]
    KG[i,:] =  sum(eachcol(pr_z(pars_GR_K, z1k, size_cen, i)))
    NG[i,:] =  sum(eachcol(pr_z(pars_NG_K, z1k, size_cen, i)))
end

p = my_summary(KG)

pS1f =plot(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
    linewidth = 3, label = false, title = "f)", titleloc = :left    
)


p = my_summary(NG)
plot!(z1k, p[:,:median], ribbon = (p.median .- p.l68, p.u68 .- p.median), 
linewidth = 3, label = false)

ylabel!("Offspring (N)")
xlabel!("Initial size \n (mm)")

scatter!(DataK.SL1_mm, DataK.Recr, groups = DataK.NG, 
c = [:lightskyblue, :red], alpha = 0.8, label = false)

xlims!(5,100)
plot(pS1a, pS1b, pS1c, pS1d, pS1e,pS1f, layout = (2,3), size = (900, 600), margin = [2mm 3mm])