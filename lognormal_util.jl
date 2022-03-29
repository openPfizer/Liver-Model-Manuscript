#=
Useful helper functions for working with LogNormal distributions and simulating % change variates.
Not used in main model simulation, the output from lit_std_sim.jl was hard-coded into main.jl.
=#

using Statistics, Distributions

# Formulas from Wikipedia:
function mu_for_mean(meanX::Vector{Float64},varX::Vector{Float64})::Vector{Float64}
    mu = log.(meanX.^2 ./ (sqrt.(varX .+ meanX.^2)))
    return mu
end

function sigma_for_var(meanX::Vector{Float64},varX::Vector{Float64})::Vector{Float64}
    sigmasq = log.(1 .+ varX./meanX.^2)
    sigma = sqrt.(sigmasq)
    return sigma
end

function lognormal_std_sim(eofx::Vector{Float64}, 
            stdofx::Vector{Float64}, rho::Float64, 
            nvar::Int64)::Vector{Float64}

# Parameters of the lognormal distribution:
muX = mu_for_mean(eofx,stdofx.^2)
sigmaX = sigma_for_var(eofx,stdofx.^2)
covX = [sigmaX[1]^2 rho*sigmaX[1]*sigmaX[2];rho*sigmaX[1]*sigmaX[2] sigmaX[2]^2]

d = MvLogNormal(muX,covX)

# Perform simulations:
r = rand(d,nvar)
rdp = 100 .* (r[2,:] .- r[1,:])./r[1,:]

return rdp

end