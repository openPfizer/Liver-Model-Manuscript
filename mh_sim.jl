#=
Helper function for the M-H search.

This function needs to know how to take an input parameter vector and ODE function
    and return a score the M-H algorithm to make its next move on.
=#

function mh_sim(p::Vector,OdeF::Function,
    xlb::Vector,xub::Vector,d_obs::Distribution)

    sol = solve(OdeF(p),Rodas4()) # Simulate the model, the timing etc. should be pre-definied in OdeF
    sim_obs = GetObs(sol,p) # Get the observables from the solution construct
    
    # Score the solution, if any X is out-of-bounds, default to eps()
    if all(sol[:,end] .>= xlb) && all(sol[:,end] .<= xub)
        q = pdf(d_obs,sim_obs)
    else
        q = eps()
    end

    return q
end

function GetObs(sol,p::Vector)
    # Helper function that extracts the relevant part of the solution structure, simple in our case.
    sim_obs = [Float64.(getLTGprct(sol[:,end],p)[1]),sol[5,end]]
    return sim_obs
end
