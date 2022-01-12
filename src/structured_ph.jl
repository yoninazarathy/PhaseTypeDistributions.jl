"""

Returns parameters of a two phase hyper-exponential fitting a mean and an SCV.

"""
function hyper_exp_init(mean_desired::Float64, scv_desired::Float64)::PHDist
    scv_desired < 1.0 && error("SCV must be greater than 1")
    μ1 = 1/(scv_desired+1) #mean parameter 
    p2 = (scv_desired-1)/(scv_desired+1+2/(μ1^2)-4/μ1)
    μ2 = (1-p2)/(1-p2/μ1) #mean parameter
    α = zeros(2)'
<<<<<<< HEAD
    α[1] = p
    α[2] = 1-p

    T = zeros(2,2)
    T[1,1] = -1/μ1
    T[2,2] = -1/μ2
    # @show p/μ1+(1-p)/μ2,μ1,μ2,m1, 2*(1-μ1)/(2-μ1*(scv+1)), -α*(mean)*T*ones(2,1)
    return PHDist(α, (1/mean)*T^(-1))
=======
    α[1] = 1-p2
    α[2] = p2

    mean_temp = (1-p2)*μ1 + p2*μ2
    T = Diagonal([-1/μ1, -1/μ2])

    return PHDist(α, (mean_temp/mean_desired)*T)
>>>>>>> d71935519d27b6b8e7403df7fab131c832b45a11
end

 

"""

Returns parameters of a hypo-exponential (generalized erlang) dist which is a sum of n exponentials with the last one different

"""
function hypo_exp_init(mean::Float64,scv::Float64)::PHDist

    scv ≥ 1.0 && error("SCV must be less than 1")

    n = Int(ceil(1/scv))

    ν1 = n/(1+sqrt((n-1)*(n*scv-1)))
    ν2 = ν1*(n-1)/(ν1-1)

    α = zeros(n)'
    α[1] = 1
    T = zeros(n,n)
    T[1,1] = -ν1
    T[1,2] = ν1

    for i = 2:(n-1)
        T[i,i] = -ν2
        T[i,i+1] = ν2
    end

    T[n,n] = -ν2

    return PHDist(α, (1/mean)*T) 
end

