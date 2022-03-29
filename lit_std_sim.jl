#=
Quick calculations on the literature target values.

Estimate of the percent change and variance of the percent change for our
target literature values is done via simulation using an assumed rho (correlation coefficient).
Change rho below with different prior information.

=#

include("lognormal_util.jl")

rho = 0.9 # correlation coefficient
nvar = 10000 # number of variates, should be reasonably large

# Belfort et al (Pioglitazone)
eofx1 = [21.3,10.2] # [%], E[X]
stdofx1 = [2.98805,2.29084] # [%], Sqrt(Var[X]))
r1 = lognormal_std_sim(eofx1,stdofx1,rho,nvar); # Perform simulation

# Haufe et al (Diet)
eofx2 = [16.7,9.3] # [%], E[X]
stdofx2 = [12.9,8.3] # [%], Sqrt(Var[X])
r2 = lognormal_std_sim(eofx2,stdofx2,rho,nvar); # Perform simulation

# Misc. calculation for the change in fat mass in Haufe et al.:
tot_adipose_frac = (2.04+7.9)/(2.56+9.3) # [fractional change from bsln]
lipolysis_reduction = (tot_adipose_frac)^(2/3) # [% change from bsln]