#=
Specialized Acceptance-Rejection sampling algorithm for virtual patients, 
based on original algorithm from: Allen et al. 2016. In Matlab.

=#

using SpecialFunctions, NearestNeighbors, Distributions
using CairoMakie, Makie.GeometryBasics, Colors

function select_vps(data_dist::Distribution, pp_obs::Matrix)
    npp = size(pp_obs,1)
    nobs = size(pp_obs,2)
    selecti = Vector{Bool}(undef,npp)
    pp_pdf = Vector{Float64}(undef,npp)
    for ii = 1:npp
        pp_pdf[ii] = pdf(data_dist,pp_obs[ii,:])
    end

    k = 5 # 5 nearest NearestNeighbors
    kdtree = KDTree(pp_obs')

    idxs, dists = knn(kdtree,pp_obs',k,true)
    pp_dists = median.(dists)
    pp_rho = k./nsphere_vol(Float64(nobs),pp_dists)
    pp_incld = pp_pdf./pp_rho

    # Select the VPs:
    pp_incld_scaled = pp_incld./maximum(pp_incld)
    for ii = 1:npp
        if pp_incld_scaled[ii] > rand()
            selecti[ii] = 1
        else
            selecti[ii] = 0
        end
    end

    return selecti, pp_incld
end

# CALCULATE N-DIMENSIONAL SPHERE VOLUME
function nsphere_vol(n::Float64,radius::Float64)::Float64
    # https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    vol = pi^(n/2)/gamma(n/2+1)*radius^n
    return vol
end

function nsphere_vol(n::Float64,radius::Vector{Float64})::Vector{Float64}
    # https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    vol = pi.^(n/2)./gamma(n/2+1) .* radius.^n
    return vol
end

# HELPER PLOTTING FUNCTIONS
function plot_pprob(p::Vector{Float64})
    Fig_P = Figure()
    g1 = Fig_P[1,1] = GridLayout()
    g2 = Fig_P[1,2] = GridLayout()

    hg1 = Axis(g1[1,1], xlabel = "prob. incld.")
    density!(hg1,p)

    hg2 = Axis(g2[1,1],xlabel = "log(prob. incld.)")
    density!(hg2,log.(p))
    current_figure()
end

# Generates Figure 2 in the manuscript, but included here since it can be made more general for any future Vpop.
# Presently many of the labels are hard-coded for the particular liver model though.
function plot_select_vs_rand(dist::Distribution,pp_obs::Matrix,
        isel::Vector{Bool},fig_name::String,
        size_pt)
    npp = size(pp_obs,1)

    iNAFLD = (isel .&& pp_obs[:,1] .> 5)

    ellipxy = errorellipse(dist)

    fac_red = Int64(round(minimum([1,npp/1000]))) # Arbitrary reduction factor just for plotting
    t = 1:fac_red:npp   
    pp_obs_plot = pp_obs[t,:]

    r = rand(dist,2000) # Draw some random variates from the actual distribution
    fig_compare = Figure(resolution = size_pt)

    labels = ["All VPs","NAFLD VPs","Data"]

    ga = fig_compare[1,1] = GridLayout()
    gb = fig_compare[2,1] = GridLayout()
    gc = fig_compare[2,2] = GridLayout()
    gl = fig_compare[1,2] = GridLayout()
    
    axga = Axis(ga[1,1],limits=(-3,5,0,0.8))
    axgb = Axis(gb[1,1], xlabel="log(Liver Fat, %)", 
        ylabel="log(Mean Plasma TG, mM)", 
        limits = (-3,5,-1.5,2.5))
    axgc = Axis(gc[1,1],limits=(0,0.8,-1.5,2.5))
    #poly!(axga,Point2f[(log(5),0),(5,0),(5,0.8),(log(5),0.8)],color = (:red,0.1))
    density!(axga,log.(r[1,:]),color=(:black,0),strokecolor=:black,strokewidth=2)
    density!(axga,log.(pp_obs[isel,1]),color=(:blue,0),strokecolor=:blue,strokewidth=2)
    density!(axga,log.(pp_obs[iNAFLD,1]),color=(:red,0),strokecolor=:red,strokewidth=2)
    
    #poly!(axgb,Point2f[(log(5),-1.5),(5,-1.5),(5,2.5),(log(5),2.5)],color = (:red,0.1))
    lines!(axgb,ellipxy[:,1],ellipxy[:,2],label = labels[3],color=:black,linewidth=2)
    scatter!(axgb,log.(pp_obs[isel,1]),log.(pp_obs[isel,2]), label = labels[1], color=(:blue))
    scatter!(axgb,log.(pp_obs[iNAFLD,1]),log.(pp_obs[iNAFLD,2]), color = (:red,1),label = labels[2])
    
    density!(axgc,log.(r[2,:]),direction=:y,color=(:black,0),strokecolor=:black,strokewidth=2)
    density!(axgc,log.(pp_obs[isel,2]),direction=:y,color=(:blue,0),strokecolor=:blue,strokewidth=2)
    density!(axgc,log.(pp_obs[iNAFLD,2]),direction=:y,color=(:red,0),strokecolor=:red,strokewidth=2)
 
    leg = Legend(fig_compare[1,2],axgb)
    hidedecorations!(axga,grid=false)
    hidedecorations!(axgc,grid=false)
    leg.tellheight = true
    rowgap!(fig_compare.layout,1,Fixed(10))
    colgap!(fig_compare.layout,1,Fixed(10))
    save(fig_name,fig_compare, pt_per_unit = 1)
end