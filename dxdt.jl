function dxdt!(du,u,p,t)
    # dx/dt function for the ODEs of the liver model.
    # Takes standard Julia Differential Equations inputs of du, u, p, and t 
    # Returns du (dx/dt) as a Vector{Float64}
    #

    # Extensions:
    #   Possible performance enhancement via StaticArrays

    ## Parameters, de-aliased:
    klipase_clear       = p[1] # [days]
    sens_nefa_uptake    = p[2] # [dimensionless]
    sens_betaox_dnl     = p[3] # [dimensionless]
    kuptake_er          = p[4] # [1/day]
    kuptake_liver_tg    = p[5] # [1/day]
    ksynth_cyt_tg       = p[6] # [1/mM^2*day]
    klipo_cyt_tg        = p[7] # [1/day]
    ksynth_er_tg        = p[8] # [1/mM^2*day]
    kbetaox             = p[9] # [1/day]
    emax_vldl_prod      = p[10] # [mmol-TG/day]
    ec50_vldl_prod      = p[11] # [mM-TG-ER]
    chylo_basal_flux    = p[12] # [mmol-TG/day]
    dnl_basal_flux      = p[13] # [mmol-FA/day]
    nefa_uptake_flux    = p[14] # [mmol-FA/day]
    vd_tg_p             = p[15] # [L]
    vd_cyt              = p[16] # [L]
    vd_er               = p[17] # [L]
    scale_chylo         = p[18] # [dimensionless]
    scale_dnl           = p[19] # [dimensionless]
    scale_nefa_uptake   = p[20] # [dimensionless]
    scale_tg_ester      = p[21] # [dimensionless]
    scale_vldl_prod     = p[22] # [dimensionless]
    fa_cy_basal         = p[23] # [mM-FA-cyt], needed for feedback loop
    tg_cy_basal         = p[24] # [mM-TG-cyt]
    fa_er_basal         = p[25] # [mM-FA-ER]
    tg_er_basal         = p[26] # [mM-TG-ER]
    tg_p_basal          = p[27] # [mM-TG-plasma]
    
    ## Concentrations, de-aliased:
    fa_cy = u[1] # [mM-cytosol]
    tg_cy = u[2] # [mM-cytosol]
    fa_er = u[3] # [mM-ER]
    tg_er = u[4] # [mM-ER]
    tg_p  = u[5] # [mM-plasma]

    ## Derived values that are "optional" to use:
    # Approximate a maximum function, to avoid an occasional error message:
    h(k,x,minx) = 1/ (1+exp(-2*k*(x- minx)))
    denom = (fa_cy - fa_cy_basal/10)*h(20,fa_cy,fa_cy_basal/10) + fa_cy_basal/10
    nefa_uptake_feedback = (fa_cy_basal/denom)^sens_nefa_uptake # [dimensionless], fix sens_nefa_uptake = 0 to ignore
    
    ## Fluxes:
    nefa_uptake     = nefa_uptake_flux*nefa_uptake_feedback*scale_nefa_uptake # [mmol-FA-cyt/day]
    er_uptake       = kuptake_er*fa_cy                                              # [mM-FA-cyt/day]
    tg_cy_ester     = ksynth_cyt_tg *fa_cy^3*scale_tg_ester                         # [mM-FA-cyt/day]
    tg_cy_lipo      = klipo_cyt_tg*tg_cy                                            # [mM-TG-cyt/day]
    tg_p_liver_uptake = kuptake_liver_tg*tg_p                                       # [mM-TG-plasma/day]
    tg_er_ester     = ksynth_er_tg*fa_er^3                                          # [mM-FA-ER/day]
    tg_er_lipo      = 0.0                                                           # [mM-TG-ER/day], Assumption*
    tg_er_vldl_prod = emax_vldl_prod*tg_er/(ec50_vldl_prod+tg_er)*scale_vldl_prod   # [mmol-TG/day]
    tg_p_clear      = klipase_clear *tg_p                                           # [mM-TG-plasma/day]
    chylo_prod      = chylo_basal_flux*scale_chylo                                  # [mmol-TG/day]
    dnl             = dnl_basal_flux*scale_dnl                                      # [mmol-FA/day]

    beta_ox_feedback = (dnl_basal_flux/dnl)^sens_betaox_dnl                         # [dimensionless]

    beta_ox         = kbetaox*fa_cy*beta_ox_feedback                                # [mM-FA-cyt/day]

    ## dxdt as fluxes:
    dfa_cy  = (nefa_uptake/vd_cyt - er_uptake - tg_cy_ester + 3*tg_cy_lipo 
                    - beta_ox + 3*tg_p_liver_uptake*vd_tg_p/vd_cyt +
                    dnl/vd_cyt)                                             # [mM-cytosol/day]
    dtg_cy  = tg_cy_ester/3 - tg_cy_lipo                                    # [mM-cytosol/day]
    dfa_er  = er_uptake*vd_cyt/vd_er - tg_er_ester + 3*tg_er_lipo           # [mM-ER/day]
    dtg_er  = tg_er_ester/3 - tg_er_lipo - tg_er_vldl_prod/vd_er            # [mM-ER/day]
    dtg_p   = (tg_er_vldl_prod/vd_tg_p - tg_p_liver_uptake - 
                    tg_p_clear + chylo_prod/vd_tg_p)                        # [mM-plasma/day]

    # Re-package:
    du[1] = dfa_cy      # [mM-cytosol/day]
    du[2] = dtg_cy      # [mM-cytosol/day]
    du[3] = dfa_er      # [mM-ER/day]
    du[4] = dtg_er      # [mM-ER/day]
    du[5] = dtg_p       # [mM-plasma/day]
    
    return du
 
end # dxdt!

# *TG ER Lipolysis was assumed to low relative to the flux through synthesis and VLDL release, 
# given that ATGL is usually referred to as a cytosolic lipase (Zimmermann et al. Science. 2004).