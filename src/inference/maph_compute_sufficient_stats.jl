
#temp
function very_crude_c_solver(y::Float64,i::Int,j::Int,maph::MAPHDist)
    quadgk(u -> (maph.α*exp(maph.T*u))[i]*exp(maph.T*(y-u))*maph.T0[:,j] , 0, y, rtol=1e-8) |> first
end

function poi_bound(λ, ε)
    K = 0
    pdf = (λ^K/factorial(K))*exp(-λ)

    while cdf < 1-ε
        K += 1
        pdf += (λ^K/factorial(K))*exp(-λ)
    end
    
    return K

end

function poi_prob(k, λ)

    return (λ^k/factorial(k))*exp(-λ)

end



function uniformizaition_integral(y::Float64, v1, v2, maph::MAPHDist)
    r = maximum(abs.(diag(maph.T)))
    P = I + maph.T./r
    R = poi_bound(r*y,0.1)


    column_vs = zeros(length(v1),R)

    column_vs[:,1] = v1
    for i = 2:R
        column_vs[:,i] = P*column_vs[:,i-1]
    end

    row_vs = zeros(R,length(v2))
    row_vs[1,:] = (poi_prob(R+1,r*y)*v2)'

    for i = 2:R
        row_vs[i,:] = row_vs[i-1,:]'*P + poi_prob(i+1,r*y)*v2
    end

    row_v = reverse(row_v)
    integral = 0
    for i = 1:R
        intergral += column_v[i]*row_v[i]/r
    end

    return integral
end


# function uniformizaition_solver(y::Float64,i::Int,j::Int,maph::MAPHDist)
#     return uniformizaition_integral(y, )



function sufficient_stats(observation::SingleObs, maph::MAPHDist; c_solver = very_crude_c_solver)::MAPHSufficientStats
    stats = MAPHSufficientStats(maph)

    p, q = model_size(maph)

    a(y::Float64) = maph.α*exp(maph.T*y)
    b(y::Float64,j::Int) = exp(maph.T*y)*maph.T0[:,j]
    c(y::Float64,i::Int,j::Int) = very_crude_c_solver(y,i,j,maph)

    D = Diagonal(maph.T)
    PT = I-inv(Diagonal(maph.T))*maph.T
    PT0 = -inv(Diagonal(maph.T))*maph.T0
    A = inv(I-PT)*PT0
    PA = maph.α*A

    EB(y::Float64, i::Int, j::Int) = maph.α[i] * b(y,j)[i]*PA[j]/ (maph.α*b(y,j))
    EZ(y::Float64, i::Int, j::Int) = c(y,i,j)[i]*PA[j]/(maph.α*b(y,j))
    ENT(y::Float64,i::Int,k::Int,j::Int) = i !=k ? maph.T[i,:].*c(y,i,j)*PA[j]/(maph.α*b(y,j)) : zeros(p)
    ENA(y::Float64,i::Int,j::Int) = PA[j]*a(y)[i]*maph.T0[i,j]/(maph.α*b(y,j))

    stats.B = [sum([EB(observation.y, i, j) for j = 1:q]) for i =1:p]
    stats.Z = [sum([EZ(observation.y,i,j) for j =1:q]) for i = 1:p]

    for i = 1:p
        for k = (q+1):(q+p)
            V = sum([ENT(observation.y,i,k-q,j) for j in 1:q])
            stats.N[i,k] = V[k-q]
        end

        for j = 1:q
            stats.N[i,j] = ENA(observation.y,i,j)
        end
    end

    return stats
end

sufficient_stats(data::MAPHObsData, maph::MAPHDist; c_solver = very_crude_c_solver)::MAPHSufficientStats = mean([sufficient_stats(d, maph) for d in data])
