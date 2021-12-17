function test_init()
    out = MAPHDist(4,[0.3,0.5,0.1,0.08,0.02],[2.0,3.1,2.3,4.5,0.2],[1.1,2.3,2.0,0.33,1.5])
    @show model_size(out)
    data = [rand(out) for _ in 1:10^5]
    pi_est = [count((x)->x.a == i, data)/length(data) for i in 1:5]
    μ_est = [mean(first.(filter((x)->x.a == i, data))) for i in 1:5] 
    @show pi_est
    @show μ_est
    return out
end
