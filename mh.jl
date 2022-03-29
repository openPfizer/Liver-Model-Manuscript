#=
Specialized implementation of a M-H algorithm for generating plausible patients.

Based on: Rieger et al. 2018. (Matlab implementation)

Note: we put *some of* the pieces in place for multi-threading below, but right now this algorithm only runs on a single thread.
Similar to StaticArrays and the ODE solver, there is likely a performance boost possible, but unneeded at the moment.

=#
using Base.Threads

function mh(n_pp::Integer, n_param_fit::Integer, log_flag::Bool,
        p::Vector{Float64}, plb::Vector{Float64},pub::Vector{Float64}, 
        Sim_and_Score::Function, OdeF::Function)
    
    nth = 1 # nthreads()
    n_param = size(p,1) # total number of parameters
    global p_pp = Matrix{Float64}(undef,n_param,1)

   for t in 1:nth
        mh_count = 0 # Loop counter
        pp_count = 0 # Number of PP counter
        q0 = 0.000 # Initial proposal score
        q1 = 0.000 # Next proposal score
        p0 = p # Initial parameter vector
        p1 = p # Next parameter vector

        if (mod(n_pp,nth) == 0)
            n_pp_th = n_pp/nth
        else
            if (t == nth)
                n_pp_th = n_pp - (nth-1)*floor(n_pp/nth)
            else
                n_pp_th = floor(n_pp/nth)
            end
        end
        n_pp_th = Int64(n_pp_th)

        p_pp_th = Matrix{Float64}(undef,n_param,n_pp_th)

        while (pp_count < n_pp_th)
            mh_count += 1
        
            # First time through we need to run two guesses:
            if mh_count == 1

                # Solve initial parameter proposal:
                p0 = Draw_P(p,plb,pub,n_param_fit,log_flag)
                q0 = Sim_and_Score(p0,OdeF)
            end

            p1 = Draw_P(p0,plb,pub,n_param_fit,log_flag)

            q1 = Sim_and_Score(p1,OdeF)

            # Core of the M-H algorithm, probalistic acceptance of the next guess:
            # We add onto q0 in the denominator to avoid an odd corner case where we were accepting a lot of
            # low probability PPs, or in other words q1 must be "meaningfully" better than q0 or we are not interested.
            if (q1/(1.1*q0 + 1000.0*eps()) > rand())
                # Accept the next solution, record a PP in various ways:
                pp_count += 1 # iterate the number of PPs
                p_pp_th[:,pp_count] = p1
                q0 = q1 # Reset for the next iteration
                p0 = p1 # Reset for the next iteration
            end  
        end

        global p_pp = hcat(p_pp,p_pp_th) # Save our PPs, this extra step seemed necessary if we ever wanted multi-threading.
    end
    p_pp = p_pp[:,2:end] # remove the first dummy row

    return p_pp
end

function Draw_P(p::Vector{Float64},plb::Vector{Float64},pub::Vector{Float64}, 
        n_param_fit::Integer,log_flag::Bool)

        #=
        Helper function for choosing a new P-vector based on the current P-vector.
        Presently, the new P-vector is chosen as a uncorrelated normal distribution
        centered on the current P-vector. This choice is not the only possible way of making the next decision
        and requires a choice of variance for those normal distributions. The variance is hard-coded below, but
        other choices may be more problem appropriate.
        =#

        n_param = size(p,1) # Total number of parameters
        
        # Avoids occasional problem where the unfitted parameters don't have an assigned range:
        p_tmp = p[1:n_param_fit] # Future extension: more flexible mapping can be added
        plb_tmp = plb[1:n_param_fit]
        pub_tmp = pub[1:n_param_fit]

        EofX = p_tmp
        VarofX = ((pub_tmp .- plb_tmp)./10).^2 # This piece can be customized some, hard-coded for now
        mu = log.(EofX.^2 ./ sqrt.(VarofX .+ EofX.^2))
        sigma = sqrt.(log.(1 .+ VarofX./EofX.^2))

        for ii = 1:n_param_fit
            
            if log_flag
                dtmp = truncated(LogNormal(mu[ii],sigma[ii]),plb_tmp[ii],pub_tmp[ii])
            else
                dtmp = truncated(Normal(EofX[ii],sqrt(VarofX[ii])),plb_tmp[ii],pub_tmp[ii])
                
            end
            global p_tmp[ii] = quantile(dtmp,rand())

        end
        
        # Assign the return p-vector:
        p_new = Vector{Float64}(undef,n_param)
        p_new = p # Reassign the full parameter vector, including all constant entries
        p_new[1:n_param_fit] = p_tmp # Overwrite the drawn entries with new parameter values
        
        # Debugging, flag NaNs. This case occurs when the bounds for p are not set properly:
        if any(isequal.(p_new,NaN))
            println(p_new)
        end

        return p_new
end
