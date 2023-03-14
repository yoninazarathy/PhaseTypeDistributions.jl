
#temp



function very_crude_c_solver(y::Float64, i::Int, j::Int, k::Int, maph::MAPHDist)
    quadgk(u -> (maph.α * exp(maph.T*u))[i] * (exp(maph.T*(y-u))*maph.T0[:,k])[j] , 0, y, rtol=1e-8) |> first
end




function poi_bound(λ, ε)
    K = 0
    pdf = (λ^K/factorial(K))*exp(-λ)

    while pdf < 1-ε
        K += 1
        pdf += (λ^K/factorial(K))*exp(-λ)
    end

    return K

end

function poi_prob(k, λ)

    return (λ^k/factorial(k))*exp(-λ)

end



function uniformizaition_integral(y::Float64, row_vec, col_vec, maph::MAPHDist)
    r = maximum(abs.(diag(maph.T)))
    P = I + maph.T./r
    R = poi_bound(r*y,0.1)
    

    column_vecs = zeros(length(col_vec),R)

    column_vecs[:,1] = col_vec
    for i = 2:R
        column_vecs[:,i] = P*column_vecs[:,i-1]
    end

    row_vecs = zeros(R,length(row_vec))
    row_vecs[1,:] = (poi_prob(R+1,r*y)*row_vec)

    for i = 2:R
        row_vecs[i,:] = row_vecs[i-1,:]'*P + poi_prob(i+1,r*y)*row_vec
    end

    row_vecs = row_vecs[end:-1:1,1:1:end]
    integral = 0
    for i = 1:R
        integral += row_vecs[i,:]'*column_vecs[:,i]/r
    end

    return integral
end




"""
Compute the expected value of the sufficient stats. (QQQQ maybe change name of function)
"""
function sufficient_stats(  observation::SingleObs, 
                            maph::MAPHDist; 
                            c_solver = very_crude_c_solver)::MAPHSufficientStats

    stats = MAPHSufficientStats(maph)

    m, n = model_size(maph)
    
    a(y::Float64) = maph.α * exp(maph.T*y)
    b(y::Float64, k::Int) = exp(maph.T*y) * maph.T0[:,k]
    c(y::Float64, i::Int, j::Int, k::Int) = very_crude_c_solver(y, i, j, k, maph)


    # D = Diagonal(maph.T)
    # PT = I-inv(Diagonal(maph.T))*maph.T
    # PT0 = -inv(Diagonal(maph.T))*maph.T0
    # A = inv(I-PT)*PT0
    # PA = maph.α*A

    EB(y::Float64, i::Int, k::Int) = maph.α[i] * b(y, k)[i] / (maph.α * b(y, k))
    EZ(y::Float64, i::Int, k::Int) = c(y, i, i, k) / (maph.α * b(y,k))
    ENT(y::Float64, i::Int, j::Int, k::Int) = i != j ? maph.T[i,j] .* c(y, i, j, k) / (maph.α * b(y,k)) : zeros(m)
    ENA(y::Float64, i::Int, k::Int) = a(y)[i] * maph.T0[i,k] / (maph.α * b(y,k))


    ### stop here
    stats.B = [sum([EB(observation.y, i, j) for j = 1:q]) for i =1:p]
    stats.Z = [sum([EZ(observation.y,i,j,i) for j =1:q]) for i = 1:p]

    for i = 1:p
        for k = (q+1):(q+p)
            V = sum([ENT(observation.y,i,j,k-q) for j in 1:q])
            stats.N[i,k] = V[k-q]
        end

        for j = 1:q
            stats.N[i,j] = ENA(observation.y,i,j)
        end
    end

    return stats
end

sufficient_stats(data::MAPHObsData, maph::MAPHDist; c_solver = very_crude_c_solver)::MAPHSufficientStats = mean([sufficient_stats(d, maph) for d in data])
