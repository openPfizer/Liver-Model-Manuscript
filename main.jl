#=
Main simulation script, running this script will generate all the figures from the manuscript.
By default, the generation of new virtual patients is off, and they are loaded from a file.
Due to the stochastic nature of virtual patient generation, if you re-generate, your results will be similar
to the manuscript but may be slightly qunatitatively different.

For efficiency we split out "setup.jl" from this script. Please invoke setup.jl first (include("setup.jl")) to
ensure you have the correct project environment running.
=#

# Dependencies:
using DifferentialEquations, CairoMakie, DataFrames, XLSX, Distributions, Statistics, CSV
using Trapz, LsqFit, ReadStatTables, DiffEqCallbacks
using FileIO, JLD2, Colors

include("util.jl")
include("dxdt.jl")
include("mh.jl")
include("mh_sim.jl")
include("select_vps.jl")

MAX_PP = 100_000 # Number of plausible patients
NUM_SPECIES = 5 # Number of equations in model
NUM_PARAM_FIT = 14 # Number of parameters to search in M-H
FRAC_FAST_TG = 0.8 # Future extension: parameterize separately (fitting parameter?)
MH_LOG_FIT = true  # Flag for setting M-H to use log boundaries for parameters
MH_FIT = true # Flag for doing the PP fitting vs. loading existing, setting to true is a much longer run.

# Read in parameters and their ranges:
df = CSV.read("parameters_pluto.csv",DataFrame) #DataFrame(XLSX.readtable("parameters.xlsx", "parameters")...) #
NUM_PARAM = size(df)[1] # Total number of parameters, only the first NUM_PARAM_FIT are fitted in M-H algorithm

# Initial parameter vector:
p = Vector{Float64}(undef,NUM_PARAM)
plb = Vector{Float64}(undef,NUM_PARAM)
pub = Vector{Float64}(undef,NUM_PARAM)
# Parse dataframe:
p .= df.TV[1:NUM_PARAM] # [various], most rate constants are 1/day, most fluxes are mmols/day
plb .= df.LV[1:NUM_PARAM]
pub .= df.HV[1:NUM_PARAM]

# Initial species vector:
x0 = Vector{Float64}(undef,NUM_SPECIES)
xlb = Vector{Float64}(undef,NUM_SPECIES)
xub = Vector{Float64}(undef,NUM_SPECIES)
x0 .= df.TV[(NUM_PARAM-NUM_SPECIES+1):NUM_PARAM] # [mM], generally
xlb .= df.LV[(NUM_PARAM-NUM_SPECIES+1):NUM_PARAM] # [mM], generally
xub .= df.HV[(NUM_PARAM-NUM_SPECIES+1):NUM_PARAM] # [mM], generally

# Implement a Metropolis-Hasting algorithm, we will functionalize this latter:
d_joint = joint_tg_dist(FRAC_FAST_TG) # retrieve the joint distribution

# Solve to baseline:
tspan = (0.0,10000.0) # [days]
cb = TerminateSteadyState(1e-8,1e-6) # callback to stop at steadystate, more efficient for M-H search

# Check SS, if needed (debugging only)
prob1 = ODEProblem(dxdt!,x0,tspan,p,callback=cb)
sol1 = solve(prob1,Rodas4())

######
# Metropolis-Hastings-based selection of plausible patients (PP):
F_ODE(a) = ODEProblem(dxdt!,x0,tspan,a,callback=cb)
F_score(a,OdeF) = mh_sim(a,OdeF,xlb,xub,d_joint)

if MH_FIT
    # Perform the full plausible patient search, this step is slow
    pp_p = mh(MAX_PP, NUM_PARAM_FIT, MH_LOG_FIT,
        p, plb, pub,
        F_score, F_ODE)
    FileIO.save("save_pp.jld2","pp_p",pp_p)
else
    # Skip the M-H plausible patient generation (DEFAULT)
    pp_p = FileIO.load("save_pp.jld2","pp_p")
end

# Re-simulate each new PP
pp_ss = Matrix{Float64}(undef,NUM_SPECIES,MAX_PP)
pp_obs = Matrix{Float64}(undef,MAX_PP,2) # [%, mM-plasma TG], observables as Array
for ii = 1:MAX_PP
    local pp_sol = solve(F_ODE(pp_p[:,ii]),Rodas4())
    ltg_ss = getLTGprct(pp_sol[:,end],pp_p[:,ii])
    pp_ss[:,ii] = pp_sol[:,end]
    pp_obs[ii,:] = [ltg_ss,pp_sol[5,end]]
end
selecti, pp_incld = select_vps(d_joint,pp_obs) # Select a Virtual Population

# Flag NAFLD vs. Non-NAFLD patients:
iNonNAFLD = ((pp_obs[:,1] .<= 5.0) .& (selecti))
iNAFLD = ((pp_obs[:,1] .> 5.0) .& (selecti))
NUM_NAFLD = sum(iNAFLD)
pp_NAFLD_p = pp_p[:,iNAFLD]
pp_NAFLD_ss = pp_ss[:,iNAFLD]

# PIOGLITAZONE
#######################################################################################
nefa_scalar = 1 - 0.275
p0 = p
pp_pio_ss_ltg_pcb = Vector{Float64}(undef,NUM_NAFLD)
for ii = 1:NUM_NAFLD
    ptmp = pp_NAFLD_p[:,ii]
    ptmp[20] *= nefa_scalar
    prob = ODEProblem(dxdt!,pp_NAFLD_ss[:,ii],(0.0,24.0*7),ptmp)
    sol = solve(prob,Rodas4())
    ltg_tmp_init = getLTGprct(sol[:,1],ptmp)
    ltg_tmp_final = getLTGprct(sol[:,:],ptmp)
    global pp_pio_ss_ltg_pcb[ii] =  100 * (ltg_tmp_final[end]/ ltg_tmp_init - 1)
end

# DIET
#######################################################################################
dnl_scale = 0.448
nefa_scale = 0.89
chylo_scale = 0.8
pp_diet_ss_ltg_pcb = Vector{Float64}(undef,NUM_NAFLD)
for ii = 1:NUM_NAFLD
    ptmp = pp_NAFLD_p[:,ii]
    ptmp[18] *= chylo_scale
    ptmp[19] *= dnl_scale
    ptmp[20] *= nefa_scale
    prob = ODEProblem(dxdt!,pp_NAFLD_ss[:,ii],(0.0,26.0*7),ptmp)
    pp_sol = solve(prob,Rodas4())
    ltg_tmp_init = getLTGprct(pp_sol[:,1],ptmp)
    ltg_tmp_final = getLTGprct(pp_sol[:,:],ptmp)
    global pp_diet_ss_ltg_pcb[ii] = 100 * (ltg_tmp_final[end]/ ltg_tmp_init - 1)
end

# SENSITIVITY ANALYSIS
#######################################################################################
sa_scale = LinRange(1.0,0.05,20) # Scalars for simple SA to loop over
p_idx_scale = [19,20,21,22] # Hard-coded indices for SA, can be more elegant
p_scale_inv = [false,false,false,true] # Inverts the scalar to "activation"
pp_ltg_SA = Matrix{Float64}(undef,size(p_idx_scale,1),size(sa_scale,1))
pp_fa_SA = Matrix{Float64}(undef,size(p_idx_scale,1),size(sa_scale,1))

ltg_pcb_tmp = Vector{Float64}(undef,NUM_NAFLD)
fa_pcb_tmp = Vector{Float64}(undef,NUM_NAFLD)

sa_ltg_pcb_vpop = Array{Float64}(undef,size(p_idx_scale,1),size(sa_scale,1),NUM_NAFLD)

for ii = 1:size(p_idx_scale,1)
    for jj = 1:size(sa_scale,1)
        for kk = 1:NUM_NAFLD
            ptmp = pp_NAFLD_p[:,kk]
            if p_scale_inv[ii] == true
                ptmp[p_idx_scale[ii]] *= 1/sa_scale[jj]
            else
                ptmp[p_idx_scale[ii]] *= sa_scale[jj]
            end
            prob = ODEProblem(dxdt!,pp_NAFLD_ss[:,kk],(0.0,26.0*7),ptmp)
            pp_sol = solve(prob,Rodas4())
            ltg_tmp_init = getLTGprct(pp_sol[:,1],ptmp)
            ltg_tmp_final = getLTGprct(pp_sol[:,:],ptmp)
            global ltg_pcb_tmp[kk] = 100 * (ltg_tmp_final[end]/ ltg_tmp_init - 1)
            global fa_pcb_tmp[kk] = 100 * (pp_sol[4,end]/pp_sol[4,1] - 1)
            global sa_ltg_pcb_vpop[ii,jj,kk] = ltg_pcb_tmp[kk]
        end
        global pp_ltg_SA[ii,jj] = mean(ltg_pcb_tmp)
        global pp_fa_SA[ii,jj] = mean(fa_pcb_tmp)
        
    end
end

# PLOTTING
#######################################################################################

# FIGURE 2: PLAUSIBLE PATIENTS
plot_select_vs_rand(d_joint,pp_obs,selecti,"Figure02.png")

colors = distinguishable_colors(6,[RGB(1,1,1), RGB(0,0,0)], dropseed=true)

# FIGURE 3: PIOGLITAZONE
data_pio_mean = -54 # [%], mean-delta% Liver Fat
data_pio_eb = 5.4 # [%], std-delta% Liver Fat

fig3 = Figure()
ha3 = fig3[1,1] = GridLayout()
axha3 = Axis(ha3[1,1], ylabel = "Δ% Liver Fat", 
    xticks = ([1,2],["Observed","Simulated"]),
    limits = (0.5,2.5,-100,0),
    xaxisposition = :top)
barplot!(axha3,[1.0,2.0],[data_pio_mean,mean(pp_pio_ss_ltg_pcb)],color = colors[1:2])
errorbars!(axha3,[1.0,2.0],[data_pio_mean,mean(pp_pio_ss_ltg_pcb)],
    [data_pio_eb,std(pp_pio_ss_ltg_pcb)],
    whiskerwidth = 15)
save("Figure03.png",fig3)

# FIGURE 4: DIET
data_diet_mean = -44.6 # [%], mean-delta% liver fat on 
data_diet_eb =  19.3 # [%] std-delta% liver fat on diet

fig4 = Figure()
ha4 = fig4[1,1] = GridLayout()
axha4 = Axis(ha4[1,1], ylabel = "Δ% Liver Fat", 
    xticks = ([1,2],["Observed","Simulated"]),
    limits = (0.5,2.5,-100,0),
    xaxisposition = :top)
    barplot!(axha4,[1.0,2.0],[data_diet_mean,mean(pp_diet_ss_ltg_pcb)],color = colors[1:2])
    errorbars!(axha4,[1.0,2.0],[data_diet_mean,mean(pp_diet_ss_ltg_pcb)],
        [data_diet_eb,std(pp_diet_ss_ltg_pcb)],
        whiskerwidth = 15)
save("Figure04.png", fig4)

# FIGURE 5: SA
fig5 = Figure()
labels5 = ["DNL", "NEFA Uptake", "Esterification", "VLDL Synth. (Inversely Plotted)"]
ha5 = fig5[1,1] = GridLayout()
axha5 = Axis(ha5[1,1], ylabel = "Δ% Liver Fat", xlabel = "Mean Inhibition [%]", 
    limits=(0,100,-100,0),xticks=[0,25,50,75,100],
    yticks=[-100,-75,-50,-25,0])
lines!(axha5,100.0 .- 100 .* sa_scale,pp_ltg_SA[1,:],label = labels5[1],linewidth=3,color=colors[3])
lines!(axha5,100.0 .- 100 .* sa_scale,pp_ltg_SA[2,:],label = labels5[2],linewidth=3,color=colors[4])
lines!(axha5,100.0 .- 100 .* sa_scale,pp_ltg_SA[3,:],label = labels5[3],linewidth=3,color=colors[5])
lines!(axha5,100.0 .- 100 .* sa_scale,pp_ltg_SA[4,:],label = labels5[4],linewidth=3,color=colors[6])
axislegend(framevisible=false,position=:lb)
save("Figure05.png", fig5)

# Figure S1: Pioglitazone and Diet with Box-and-Whiskers
figS1 = Figure()
haS11 = figS1[1,1] = GridLayout()
haS12 = figS1[2,1] = GridLayout()

axhaS11 = Axis(haS11[1,1], ylabel = "Δ% Liver Fat", 
    xticks = ([1,2],["Observed","Simulated"]),
    limits = (0.5,2.5,-100,0),
    xaxisposition = :top,
    title = "Pioglitazone")
barplot!(axhaS11,[1],[data_pio_mean],color = colors[1])
errorbars!(axhaS11,[1],[data_pio_mean],[data_pio_eb],whiskerwidth = 15)
violin!(axhaS11,2.0 .* ones(size(pp_pio_ss_ltg_pcb,1)),pp_pio_ss_ltg_pcb, 
    color = (colors[2]))

axhaS12 = Axis(haS12[1,1], ylabel = "Δ% Liver Fat", 
    xticks = ([1,2],["Observed","Simulated"]),
    limits = (0.5,2.5,-100,0),
    xaxisposition = :top,
    title = "Diet")
barplot!(axhaS12,[1.0],[data_diet_mean],color = colors[1])
errorbars!(axhaS12,[1.0],[data_diet_mean],
    [data_diet_eb],
    whiskerwidth = 15)
violin!(axhaS12,2.0 .* ones(size(pp_diet_ss_ltg_pcb,1)),pp_diet_ss_ltg_pcb, 
    color = (colors[2]))

noto_sans = assetpath("fonts", "NotoSans-Regular.ttf")
noto_sans_bold = assetpath("fonts", "NotoSans-Bold.ttf")   
for (label, layout) in zip(["A", "B"], [haS11, haS12])
    Label(layout[1, 1, TopLeft()], label,
    textsize = 26,
    font = noto_sans_bold,
    padding = (0, 5, 5, 0),
    halign = :left)
end

save("FigureS1.png", figS1)

# Figure S2: SA distribution
figS2 = Figure(fontsize = 12)
h = figS2[1:5,1:5] = GridLayout()

rcount=1
ccount=1
for jj = 1:size(sa_scale,1)

    global axh = Axis(h[rcount,ccount],xlabel = "Δ% Liver Fat",ylabel="Vpop PDF",
        limits = (-100,0,0,0.3),title=string(round(100.0 - 100.0*sa_scale[jj]),"%"))
    tmp = sa_ltg_pcb_vpop[:,jj,:]
    for ii = 1:size(p_idx_scale,1)
        hist!(axh,tmp[ii,:],label=labels5[ii],normalization=:pdf)
    end
    
    global ccount += 1
    if ccount > 5
        global ccount = 1
        global rcount += 1
    end
end
Legend(figS2[6,3],axh,orientation=:horizontal)

save("FigureS2.png", figS2)
