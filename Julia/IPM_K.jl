# these dataframes are from previous analyses in stan
function Killifish_IPM(post; nBigMatrix = 100, min_size = 2, max_size = 110, 
  size_cen = 18.0)
	
  # nBigMatrix = 100
  # min_size = 2.0
  # max_size = 110.0
  U= max_size
  L= min_size
  m = nBigMatrix
  h = (U - L)/m

  z1 =  L .+ (collect(1:m) .- 0.5) .* h
  z1 = round.(z1, digits = 6)
  z = z1
  meshpts = z1

  size_cen = (18.0)

  # Get the parameters

  first(post,6)
  post.sigma_growK
  pars_GR_K = select(post,
      :Intercept_survK => :"α_surv",
      :b_z_survK => :"βz_surv",
      :b_area_survK => :"β_area_surv",
      :b_canopy_survK => :"β_canopy_surv",

      :Intercept_growK => :"α_grow",
      :b_z_growK => :"βz_grow",
      :b_area_growK => :"β_area_grow",
      :b_canopy_growK => :"β_canopy_grow",
      :sigma_growK => ByRow(x-> sqrt(x)) =>:σ_grow,

      :Intercept_recrK => :"α_fec",
      :b_z_recrK => :"βz_fec",
      :b_area_recrK => :"β_area_fec",
      :b_canopy_recrK => :"β_canopy_fec"


  )



  pars_NG_K = select(post,
      :Intercept_survK => :"α_surv",
      :b_z_survK => :"βz_surv",
      :b_area_survK => :"β_area_surv",
      :b_canopy_survK => :"β_canopy_surv",

      :Intercept_growK => :"α_grow",
      :b_z_growK => :"βz_grow",
     
      :b_area_growK => :"β_area_grow",
      :b_canopy_growK => :"β_canopy_grow",
      :sigma_growK => ByRow(x-> sqrt(x)) =>:σ_grow,

      :Intercept_recrK => :"α_fec",
      :b_z_recrK => :"βz_fec",
      :b_area_recrK => :"β_area_fec",
      :b_canopy_recrK => :"β_canopy_fec"
  )



  pars_NG_K.α_surv = pars_NG_K.α_surv .+ post.b_NG_survK
  pars_NG_K.βz_surv = pars_NG_K.βz_surv .+ post.b_zNG_survK 
  pars_NG_K.α_grow = pars_NG_K.α_grow .+ post.b_NG_growK
  pars_NG_K.βz_grow = pars_NG_K.βz_grow .+ post.b_zNG_growK 
  pars_NG_K.α_fec = pars_NG_K.α_fec .+ post.b_NG_recrK
  pars_NG_K.βz_fec = pars_NG_K.βz_fec .+ post.b_zNG_recrK 



  function g_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    α= df.α_grow[row]
    β= df.βz_grow[row]
   
    σ= df.σ_grow[row]
    p_den_grow = zeros(size(z)[1],size(z)[1])
    μ = ((α .+ β .* (z .- size_cen ) .- (z))./2 .+ z) # average growth in two weeks
    #μ = round.(μ, digits = 10)
    
    for i in 1:nBigMatrix
        
      p_den_grow[:,i] =   (pdf.(Normal(μ[i], σ), z1).*h)
      
    end
    matex = zeros(nBigMatrix+1, nBigMatrix+1)
    matex[2:end,2:end] = p_den_grow
    return(matex)
  end


  #@time g = g_z1zK(pars_GR_K, z1, z, size_cen, row)


  ## Surival function


  function s_zK(df::AbstractDataFrame, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    α= df.α_surv[row]
    β= df.βz_surv[row]
   
    linear_p = α .+ β * (z .- size_cen)      # linear predictor
    p = 1 ./(1 .+ exp.(-linear_p))
    p = diagm(sqrt.(p))
    matex = zeros(nBigMatrix+1, nBigMatrix+1)
    matex[2:end,2:end] = p
    return(matex)
  end



  ## Reproduction function, logistic regression
  # function pb_z(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
  #   α= df.α_rep[row]
  #   β= df.βz_rep[row]
  #   linear_p = α .+ β * (z .- size_cen)       # linear predictor
  #   p = 1 ./(1 .+ exp.(-linear_p))
  #   p = diagm(p)
  #   return(p)
  # end



  ## Recruitment function (N.B - from birth in spring to first summer), logistic regression
  function pr_zK(df::AbstractDataFrame, z::AbstractVector, size_cen::AbstractFloat, row::Integer)
    α= df.α_fec[row]
    β= df.βz_fec[row]
    linear_p = α .+ β * (z .- size_cen)   # linear predictor  # linear predictor
    p = exp.(linear_p).*(1/2)
    #p = diagm(p)
    matex = zeros(nBigMatrix+1, nBigMatrix+1)
    matex[1,2:end] = p
    return(matex)

  end

  ## Recruit size function
  function c_z1zK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, 
    size_cen::AbstractFloat, row::Integer)
    #α= df.rcz_int[row]
    #β= df.rcz_z[row]
    #βa= df.β_area_rcz[row]
    #βc= df.β_canopy_rcz[row]
    σ= 0.6 #pars.rcz_sd[row]
    p_den_grow = zeros(nBigMatrix+1,nBigMatrix+1)

    μ = 4.4 .+ 0#*(z .- size_cen)
      #   for i in 1:nBigMatrix
      #     p_df = pdf.(Normal(μ[i], σ), z1)*h
      #     for j in 1:nBigMatrix
      #       p_den_grow[j,i] = p_df[j]
      #     end
      #   end
    p_den_grow[2:end, 1] =  pdf.(Normal(μ, σ), z).*h
    return(p_den_grow)
  end
 

  ##----------------------------------------------------
  ## Functions to build IPM kernels P, F, and K



  function mk_KK(df::AbstractDataFrame, z1::AbstractVector, z::AbstractVector, size_cen::AbstractFloat, 
    row::Integer, V::AbstractFloat)
      F = pr_zK(df, z, size_cen, row) * s_zK(df, z, size_cen, row)
      P = g_z1zK(df, z, z, size_cen, row) * s_zK(df, z, size_cen, row)
      A = (P + F) + (V .* c_z1zK(df, z1, z, size_cen, row))
      K = A * A
      out = Dict("K"=> K, "meshpts" => z, "P" => P, "F"=> F, "A" => A)
      return(out)
  end

  # Matrices for vital rates
  Ksurv_mat = zeros(size(post)[1], nBigMatrix+1)
  Kgrow_mat = zeros(size(post)[1], nBigMatrix+1)
  Krep_mat = zeros(size(post)[1], nBigMatrix+1)
  Kfec_mat = zeros(size(post)[1], nBigMatrix+1)
  Krcz_mat = zeros(size(post)[1], nBigMatrix+1)

  ## make DataFrame to store the results


  Kres_IPM = DataFrame(zeros(size(post)[1], 15), :auto)


  Kres_IPM = select(Kres_IPM, :x1 => "lam_GR", :x2 => "lam_NG", :x3 => "delta_lam",
            :x4 => "sum_lam_eff", :x5 => "grow_con", :x6 => "fec_con", 
            :x7 => "rcz_con", :x8 => "sur_con",
            :x9 => "sum_con")



  for row in 1:size(pars_GR_K)[1]
      # Make projection kernels
    IPM_GR = mk_KK(pars_GR_K, z1, z, size_cen, row, 0.7)
    IPM_NG = mk_KK(pars_NG_K, z1, z, size_cen, row, 0.7)
      # calculate the population growth rate (λ)

    vv_GR = eigen(IPM_GR["K"])
    vv_NG = eigen(IPM_NG["K"])

    λ_GR = real(vv_GR.values[end])
    λ_NG = real(vv_NG.values[end])
    Kres_IPM.lam_GR[row] = λ_GR
    Kres_IPM.lam_NG[row] = λ_NG
      λ_NG -λ_GR

      ## calculate average matrix
    K_avg = (IPM_NG["K"] + IPM_GR["K"])./2
    
    vv_avg = eigen(K_avg)

      # Normalize stable size distribution
    ω_GR = real(vv_GR.vectors[:, end]) ./ sum(real(vv_GR.vectors[:, end]))
    ω_NG = real(vv_NG.vectors[:, end]) ./  sum(real(vv_NG.vectors[:, end]))
    ω_avg = real(vv_avg.vectors[:, end]) ./ sum(real(vv_avg.vectors[:, end]))

      # Reproductive value

    a_avg = eigen(K_avg')
    v_avg = real(a_avg.vectors[:, end]) 
    v_avg = v_avg ./ dot(transpose(v_avg), ω_avg)


      ## Sensitivity matrix

    sens_avg = v_avg * ω_avg'  # this equivalent to outer() in R
    ΔK = IPM_NG["K"] .- IPM_GR["K"]
    λ_eff = ΔK .* sens_avg
    Δλ = sum(λ_eff)
    

    Kres_IPM.sum_lam_eff[row] = Δλ
    Kres_IPM.delta_lam[row] = λ_NG - λ_GR

    ## Make life-response table

    #one_mat = ones(nBigMatrix+1, nBigMatrix+1)

      # Function differences
    Δ_grow = g_z1zK(pars_NG_K, z1, z, size_cen, row) .- g_z1zK(pars_GR_K, z1, z, size_cen, row)
    #Δ_rep = one_mat*(pb_z(pars_NG_K, z, size_cen, row) - pb_z(pars_GR_K, z, size_cen, row))
    Δ_fec = (pr_zK(pars_NG_K, z, size_cen, row) .- pr_zK(pars_GR_K, z, size_cen, row))
    Δ_rcz = (c_z1zK(pars_NG_K, z1, z, size_cen, row) .- c_z1zK(pars_GR_K, z1, z, size_cen, row))
    Δ_sur = (s_zK(pars_NG_K, z, size_cen, row) .- s_zK(pars_GR_K, z, size_cen, row))

      # Function averages
    grow_avg = (g_z1zK(pars_NG_K, z1, z, size_cen, row) + g_z1zK(pars_GR_K, z1, z, size_cen, row))/2
    #rep_avg = (one_mat*(pb_z(pars_NG_K, z, size_cen, row) + pb_z(pars_GR_K, z, size_cen, row)))/2
    fec_avg = ((pr_zK(pars_NG_K, z, size_cen, row) + pr_zK(pars_GR_K, z, size_cen, row)))/2
    rcz_avg = ((c_z1zK(pars_NG_K, z1, z, size_cen, row) + c_z1zK(pars_GR_K, z1, z, size_cen, row)))/2
    sur_avg = ((s_zK(pars_NG_K, z, size_cen, row) + s_zK(pars_GR_K, z, size_cen, row)))/2

      # derivates
    I_mat = diagm(ones(nBigMatrix+1))


    𝛿_grow = kron(transpose(sur_avg), I_mat)
    𝛿_sur  = kron(transpose(I_mat),grow_avg) +  kron(transpose(I_mat),(fec_avg * rcz_avg))
    𝛿_fec  = kron(transpose(sur_avg), rcz_avg)
    𝛿_rcz  = kron(transpose(fec_avg * sur_avg), I_mat)

    𝛿A = kron(I_mat, K_avg) + kron(transpose(K_avg), I_mat)

    λ_grow = Δ_grow .* reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_grow), (nBigMatrix+1,nBigMatrix+1))
    λ_fec = Δ_fec .*  reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_fec), (nBigMatrix+1,nBigMatrix+1))
    λ_rcz = Δ_rcz .*  reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_rcz), (nBigMatrix+1,nBigMatrix+1)) 
    λ_sur = Δ_sur .* reshape((([transpose(sens_avg[:]) ; ] * 𝛿A )* 𝛿_sur), (nBigMatrix+1,nBigMatrix+1)) 



      # to put in a DataFrame
    sur_con = sum(λ_sur)
    grow_con = sum(λ_grow)
    #rep_con = sum(λ_rep)
    fec_con = sum(λ_fec)
    rcz_con = sum(λ_rcz)
    sum_con = sur_con + grow_con + fec_con + rcz_con
    Kres_IPM.sum_con[row] = sum_con


    Kres_IPM.sur_con[row] = sur_con
    Kres_IPM.grow_con[row] = grow_con
    #Kres_IPM.rep_con[row] = rep_con
    Kres_IPM.fec_con[row] = fec_con
    Kres_IPM.rcz_con[row] = rcz_con


    Ksurv_mat[row, : ] =  sum.(eachcol(λ_sur)) 
    Kgrow_mat[row, : ] =  sum.(eachcol(λ_grow)) 
    #Krep_mat[row, : ] =  sum.(eachcol(λ_rep)) 
    Kfec_mat[row, : ] =  sum.(eachcol(λ_fec)) 
    Krcz_mat[row, : ] =  sum.(eachcol(λ_rcz)) 

  end
  return [Kres_IPM, Ksurv_mat, Kgrow_mat, Kfec_mat, Krcz_mat]

end