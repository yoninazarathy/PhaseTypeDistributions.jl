function E_step(all_obs::Vector{SingleObservation}, maph::MAPHDist)
    data_length, sorted_stats = compute_sorted_stats(all_obs, maph)

    sum_stats_per_state = map(ss -> sum(ss), sorted_stats)
    @show sum_stats_per_state
    total_stats = sum(sum_stats_per_state)

    return data_length, sum_stats_per_state, total_stats
end

function M_step(data_length::Int, stats::Vector{MAPHSufficientStats}, sum_stats::MAPHSufficientStats, maph::MAPHDist)
    m, n = model_size(maph)
    number_times_exit = [sum_stats.N[i] + sum(sum_stats.M[i,:]) for i ∈ 1:m]
    number_times_exit_per_state = [[stats[k].N[i] + sum(stats[k].M[i, :]) for i ∈ 1:m] for k ∈ 1:n]
    new_α =  sum_stats.B ./data_length
    new_q = [number_times_exit[i] / sum_stats.Z[i] for i ∈ 1:m]
    new_R = [stats[k].B[i] / sum_stats.B[i] for i = 1:m, k = 1:n]
    new_U = [stats[1].M[i,j] / number_times_exit_per_state[1][i] for i ∈ 1:m, j ∈ 1:m ]
    @assert sum(new_α) ≈ 1.0
    @assert sum(new_α .* new_ρ) ≈ 1.0
    @assert all(new_q .≥ 0)  "we have $new_q"
    maph.α = reshape(new_α, (1,m)) 

    if !satisfies_constraint_U(new_R, new_U)
        @info "we are going to project a new U matrix here since the fitted one does not statisfy the constraint"
        new_U = project_U_step(R, U)
    end

    return MAPHDist(new_α, new_q, new_R, new_U)

end

function project_U_step(R::Matrix, target_U::Matrix)
    m, n = size(R)

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, u[i=1:m, j = 1:m] >= 0)
    @variable(model, d[i=1:m, j =1:m] >= 0)
    
    @constraint(model, [i=1:m, j = 1:m], d[i,j] >= u[i,j] - target_U[i,j])
    @constraint(model, [i=1:m, j = 1:m], d[i,j] >= -(u[i,j] - target_U[i,j]))
    @constraint(model, [i=1:m, k = 1:n], sum((R[j,k]/R[j,1]) * u[i,j] for j = 1:m) <= (R[i,k]/R[i,1]))

    @objective(model, Min, sum(d[i,j] for i in 1:m, j in 1:m))
    optimize!(model)

    model_status = termination_status(model)
    uvals = [value(u[i,j]) for i in 1:m, j in 1:m]
    U = 
    if model_status == MOI.OPTIMAL
        @info "Optimal solution found!"
        uvals
    elseif model_status in [MOI.FEASIBLE_POINT, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL]
        @info "Suboptimal solution found. Extracting best known values."
        uvals
    elseif model_status == MOI.INFEASIBLE 
        @info "The problem is infeasible. No feasible solution exists."
        nothing  # No valid solution
    elseif model_status == MOI.UNBOUNDED
        @info "The problem is unbounded. The objective function can decrease indefinitely."
        nothing
    elseif model_status == MOI.INFEASIBLE_OR_UNBOUNDED
        @info "The problem is either infeasible or unbounded."
        nothing
    else
        @info "Optimization did not fully converge. Status: $model_status"
        nothing  # Handle unexpected cases
    end
    return U

end


function EM_fit(all_obs::Vector{SingleObservation}, maph::MAPHDist, iterations::Int = 20)
    maph_out = maph
    for _ ∈ 1:iterations
        @time data_length, stats, s_stats = E_step(all_obs, maph)
        @time maph_out = M_step(data_length, stats, s_stats, maph)
    end
    return maph_out
end

   