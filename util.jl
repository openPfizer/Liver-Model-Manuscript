# Handful of useful, but misc. functions that helped along the way.

using StatsBase, Statistics, LinearAlgebra, Distributions

function getLTGprct(u::Vector{Float64},p)
    # Converts mM --> % liver fat on volumetric basis.
    # Future extension: error check that vd_cyt and vd_er are correct with p-vector.
    rho_tg = 0.9*1000.0 # [g/L], includes /cm^3 --> /L conversion
    mw_tg = 772
    vd_cyt = p[16] # [L]
    vd_er = p[17] # [L]
    hepat_vol = 0.95*0.7089 # [L], estimate of non-LTG hepatocyte volume, see Pluto notebook: derived_parameters.jl
    ltg_mass = (u[2]*vd_cyt+u[4]*vd_er)*mw_tg/1000.0 # [g-TG], includes mg-->g conversion
    ltg_vol = ltg_mass/rho_tg # [L-TG]
    ltg_prct = 100*ltg_vol/(ltg_vol + hepat_vol) # [%]
    return ltg_prct
end

function getLTGprct(u::Matrix{Float64},p)
    # Converts mM --> % liver fat on volumetric basis.
    # Future extension: error check that vd_cyt and vd_er are correct with p-vector.
    rho_tg = 0.9*1000.0 # [g/L], includes /cm^3 --> /L conversion
    mw_tg = 772
    vd_cyt = p[16] # [L]
    vd_er = p[17] # [L]
    hepat_vol = vd_cyt + vd_er # [L]
    ltg_mass = (u[2,:] .* vd_cyt .+ u[4,:] .* vd_er) .* mw_tg/1000.0 # [g-TG], includes mg-->g conversion
    ltg_vol = ltg_mass ./ rho_tg # [L-TG]
    ltg_prct = 100 .* ltg_vol ./ (ltg_vol .+ hepat_vol) # [%]
    return ltg_prct
end

function joint_tg_dist(frac_fast_tg)
    # Creates a MvLogNormal distribution from the data below based on the mean and covariance of the data.
   
    # Data from captured from:
    # Anna Kotronen, Jukka Westerbacka, Robert Bergholm, Kirsi H. Pietiläinen, Hannele Yki-Järvinen, 
    # Liver Fat in the Metabolic Syndrome, 
    # The Journal of Clinical Endocrinology & Metabolism, 
    # Volume 92, Issue 9, 1 September 2007, Pages 3490–3497, 
    # https://doi.org/10.1210/jc.2007-0482
    x = [1.98E-01	6.78E-01
    4.93E-01	4.12E-01
    3.93E-01	6.47E-01
    3.86E-01	7.71E-01
    4.84E-01	5.30E-01
    4.92E-01	5.78E-01
    5.00E-01	6.24E-01
    5.00E-01	6.98E-01
    4.91E-01	7.91E-01
    4.91E-01	8.64E-01
    4.82E-01	9.43E-01
    4.90E-01	1.02E+00
    4.73E-01	1.15E+00
    4.89E-01	1.27E+00
    4.80E-01	1.72E+00
    8.46E-01	1.81E+00
    8.64E-01	1.10E+00
    9.64E-01	3.47E-01
    1.01E+00	4.69E-01
    1.01E+00	5.25E-01
    9.94E-01	6.17E-01
    1.23E+00	2.06E+00
    1.49E+00	2.33E+00
    4.29E+00	1.00E+01
    6.70E+00	9.92E+00
    3.89E+00	3.72E+00
    1.93E+00	2.19E+00
    2.45E+00	2.22E+00
    3.07E+00	2.25E+00
    2.91E+00	2.40E+00
    3.01E+00	2.90E+00
    2.00E+00	1.91E+00
    1.90E+00	1.77E+00
    1.93E+00	1.55E+00
    1.66E+00	1.52E+00
    1.42E+00	1.54E+00
    1.47E+00	1.25E+00
    1.57E+00	1.26E+00
    1.45E+00	9.71E-01
    1.45E+00	8.78E-01
    1.50E+00	7.75E-01
    1.48E+00	6.92E-01
    1.48E+00	6.18E-01
    1.48E+00	5.39E-01
    2.01E+00	5.06E-01
    1.94E+00	6.11E-01
    1.94E+00	7.10E-01
    1.94E+00	7.85E-01
    1.97E+00	8.79E-01
    1.94E+00	9.96E-01
    1.90E+00	1.09E+00
    1.94E+00	1.22E+00
    1.90E+00	1.38E+00
    2.42E+00	1.65E+00
    2.42E+00	1.02E+00
    2.92E+00	1.69E+00
    2.97E+00	1.47E+00
    9.93E-01	7.27E-01
    9.93E-01	8.13E-01
    9.75E-01	8.99E-01
    9.92E-01	9.57E-01
    9.73E-01	1.21E+00
    9.73E-01	1.29E+00
    9.89E-01	1.47E+00
    1.01E+00	1.58E+00
    9.88E-01	1.75E+00
    9.55E-01	1.86E+00
    2.87E+00	1.30E+00
    2.87E+00	1.22E+00
    2.93E+00	1.09E+00
    2.83E+00	1.06E+00
    2.98E+00	9.37E-01
    2.88E+00	8.48E-01
    2.93E+00	7.48E-01
    3.36E+00	7.12E-01
    3.86E+00	7.68E-01
    4.43E+00	7.49E-01
    5.00E+00	7.04E-01
    5.93E+00	6.62E-01
    5.83E+00	7.69E-01
    5.93E+00	8.19E-01
    5.08E+00	9.04E-01
    4.50E+00	8.71E-01
    3.92E+00	9.99E-01
    3.66E+00	9.15E-01
    3.30E+00	9.61E-01
    3.41E+00	1.09E+00
    3.97E+00	1.87E+00
    4.63E+00	2.20E+00
    4.87E+00	2.56E+00
    4.87E+00	2.69E+00
    5.79E+00	2.66E+00
    5.79E+00	2.18E+00
    6.32E+00	1.92E+00
    5.70E+00	1.97E+00
    5.90E+00	1.80E+00
    5.80E+00	1.65E+00
    4.81E+00	1.46E+00
    4.34E+00	1.38E+00
    4.50E+00	1.09E+00
    5.43E+00	1.12E+00
    5.34E+00	1.08E+00
    5.91E+00	1.22E+00
    5.92E+00	9.76E-01
    6.91E+00	9.17E-01
    6.80E+00	8.29E-01
    8.82E+00	5.09E-01
    9.28E+00	6.71E-01
    1.06E+01	7.90E-01
    9.59E+00	8.73E-01
    1.03E+01	1.02E+00
    1.12E+01	1.09E+00
    1.35E+01	9.54E-01
    1.84E+01	8.43E-01
    2.47E+01	7.26E-01
    2.23E+01	4.51E-01
    1.64E+01	5.39E+00
    3.88E+01	6.69E+00
    3.33E+01	4.59E+00
    2.21E+01	3.89E+00
    1.62E+01	4.14E+00
    1.46E+01	3.39E+00
    1.37E+01	2.54E+00
    2.10E+01	2.71E+00
    3.34E+01	2.45E+00
    4.32E+01	2.30E+00
    3.71E+01	1.50E+00
    3.85E+01	1.02E+00
    2.73E+01	8.76E-01
    2.51E+01	1.04E+00
    2.64E+01	1.10E+00
    2.11E+01	1.03E+00
    1.84E+01	1.10E+00
    1.97E+01	1.11E+00
    2.22E+01	1.20E+00
    2.11E+01	1.37E+00
    1.93E+01	1.31E+00
    1.52E+01	1.15E+00
    1.42E+01	1.29E+00
    1.37E+01	1.41E+00
    1.60E+01	1.41E+00
    1.80E+01	1.48E+00
    2.26E+01	1.68E+00
    2.72E+01	1.98E+00
    2.91E+01	2.08E+00
    3.24E+01	1.56E+00
    3.52E+01	1.77E+00
    2.82E+01	1.68E+00
    2.77E+01	1.54E+00
    2.59E+01	1.43E+00
    2.68E+01	1.36E+00
    2.25E+01	1.91E+00
    2.21E+01	2.05E+00
    2.00E+01	1.90E+00
    1.74E+01	1.77E+00
    8.06E+00	1.19E+00
    8.48E+00	1.30E+00
    9.57E+00	1.16E+00
    1.06E+01	1.27E+00
    1.16E+01	1.30E+00
    9.40E+00	1.29E+00
    1.08E+01	1.48E+00
    7.38E+00	1.53E+00
    6.22E+00	1.48E+00
    7.25E+00	1.76E+00
    7.37E+00	2.05E+00
    8.17E+00	2.23E+00
    7.24E+00	2.56E+00
    9.55E+00	1.78E+00
    1.06E+01	1.92E+00
    1.15E+01	1.81E+00
    1.26E+01	1.79E+00
    1.23E+01	1.95E+00
    1.24E+01	1.58E+00
    1.35E+01	1.52E+00
    3.47E+00	1.24E+00] # Liver Fat [%], TG [mM]

    x[:,2] = x[:,2]./frac_fast_tg
    xl = log.(x)
    mu = vec(mean(xl,dims=1))
    cov = StatsBase.cov(xl)
    d_joint = MvLogNormal(mu,cov) # Multivariate log-normal distribution
    
    return d_joint
    
    end
    
    function errorellipse(d::Distribution)
        # Returns a 2D error ellipse at the 95% prediction interval
        # Future exensions:
        #   Expansion to 3D+ (projections)
        #   Assumes distribution is MvLogNormal

        nvar = 100000 # number of points to estimate ellipse using
        nplot = nvar # number of points for plotting
        r = log.(rand(d,nvar)) # draw variates
        rcov = StatsBase.cov(r')
        reigval = LinearAlgebra.eigvals(rcov)
        reigvec = LinearAlgebra.eigvecs(rcov)

        reigval_max = findmax(reigval)[1]
        reigvec_max = reigvec[:,findmax(reigval)[2]]
        reigval_min = findmin(reigval)[1]
        reigvec_min = reigvec[:,findmin(reigval)[2]]
        rangle = atan(reigvec_max[2],reigvec_max[1])

        if (rangle < 0)
            rangle = rangle + 2*pi
        end

        rmean = mean(r,dims=2)

        chisquare_val = 2.4477
        theta_grid = LinRange(0,2*pi,nplot)
        phi = rangle
        X0 = rmean[1]
        Y0 = rmean[2]
        a = chisquare_val*sqrt(reigval_max)
        b = chisquare_val*sqrt(reigval_min)

        R = [cos(phi) sin(phi);
            -sin(phi) cos(phi)]

        ellipsexy = Matrix{Float64}(undef,nplot,2)
        ellipsexy[:,1] = a*cos.(theta_grid)
        ellipsexy[:,2] = b*sin.(theta_grid)

        ellipsexy = ellipsexy * R
        ellipsexy[:,1] = ellipsexy[:,1] .+ X0
        ellipsexy[:,2] = ellipsexy[:,2] .+ Y0

        return ellipsexy

    end


    function LogNormFromUT(MeanX::Float64,LV::Float64,HV::Float64)
        # Univariate LogNormal from untransformed values
        
        sigma = log(HV/LV)/3.2897072539029435
        #VarX = StdX^2
        mu = log(MeanX) - sigma^2/2#log(MeanX^2/sqrt(VarX + MeanX^2))
        #sigma = sqrt(log(VarX/MeanX^2+1))

        return LogNormal(mu,sigma)
    end