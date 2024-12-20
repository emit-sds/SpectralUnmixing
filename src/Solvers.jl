#  Copyright 2022 California Institute of Technology
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
# Author: Philip G. Brodrick, philip.g.brodrick@jpl.nasa.gov

using JuMP
using NLopt
using LinearAlgebra

"""
    opt_solve(A::Matrix{Float64}, b::Vector{Float64}, x0::Vector{Float64},
    lb::Vector{Float64}, ub::Vector{Float64})

Solve a nonlinear least squares problem, finding the vector `x` subject to lower (`lb`)
and upper (`ub`) bounds, that minimizes the residual ||Ax - b||² using the NLopt library's
implementation of the SLSQP algorithm.

# Arguments
- `A::Matrix{Float64}`: Coefficient matrix of size (m, n).
- `b::Vector{Float64}`: Target vector of size (m,).
- `x0::Vector{Float64}`: Initial guess vector of size (n,) for the optimization variables.
Values will be clipped to the range [0, 1].
- `lb::Vector{Float64}`: Vector of size (n,) specifying the lower bounds for the
optimization variables.
- `ub::Vector{Float64}`: Vector of size (n,) specifying the upper bounds for the
optimization variables.

# Returns
- A tuple containing:
  - `x`: The optimized values of the variables `x` (vector of size (n,)).
  - `mse_opt`: The exponential of the objective function value at the optimum.
"""
function opt_solve(A::Matrix{Float64}, b::Vector{Float64}, x0::Vector{Float64},
    lb::Vector{Float64}, ub::Vector{Float64})

    #mle = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    mle = Model(NLopt.Optimizer)
    set_optimizer_attribute(mle, "algorithm", :LD_SLSQP) #LD_SLSQP

    x0[x0.<0] .= 0
    x0[x0.>1] .= 1

    @variable(mle, lb[i] <= x[i=1:length(x0)] <= ub[i])
    #@constraint(mle, sum(x) == 1)
    for n in 1:length(x0)
        set_start_value(x[n], x0[n])
    end

    #@NLexpression(mle, mse ,  sum( (b .- sum(x[i] .* A[:,i] for i in 1:length(x) )).^2    ) )
    @expression(mle, mse, sum((b .- sum(x[i] .* A[:, i] for i in 1:length(x))) .^ 2))

    @NLobjective(mle, Min, mse)

    JuMP.optimize!(mle)
    return value.(mle[:x]), exp(objective_value(mle))
end

"""
    dolsq(A, b; method::String="default")

Solve a least squares problem, finding the vector `x` that minimizes the residual
||Ax - b||².

# Arguments
- `A`: Coefficient matrix of size (m, n).
- `b`: Target vector of size (m,).
- `method`: An optional keyword argument specifying the solver method.
Available options are:
    - `"default"`: Solve the system using the backslash operator (`\`).
    - `"pinv"`: Use the pseudoinverse of `A` to compute the solution.
    - `"qr"`: Use QR decomposition to solve the least squares problem.

# Returns
- `x`: Vector of size (n,) that minimizes the least squares residual.
"""
function dolsq(A, b; method::String="default")
    if method == "default"
        x = A \ b
    elseif method == "pinv"
        x = pinv(A) * b
    elseif method == "qr"
        #Q,R = qr(A)
        #x = inv(R)*Q'*b
        qrA = qr(A)
        x = qrA \ b
    end
    return x
end

"""
    bvls(A, b, x_lsq, lb, ub, tol::Float64, max_iter::Int64, verbose::Int64,
         inverse_method::String)

Solve a bounded value least squares problem, finding the vector `x` subject to lower (`lb`)
and upper (`ub`) bounds, that minimizes the residual ||Ax - b||².

- The function employs an iterative approach, adjusting the solution based on the specified
  bounds.
- Convergence is determined by whether KKT optimality condition (see
[`compute_kkt_optimality`](@ref)) is within the specified tolerance.

# Arguments
- `A`: Coefficient matrix of size (m, n).
- `b`: Target vector of size (m,).
- `x_lsq`: Initial guess vector of size (n,) for the optimization variables.
- `lb`: Vector of size (n,) specifying the lower bounds for the
optimization variables.
- `ub`: Vector of size (n,) specifying the upper bounds for the
optimization variables.
- `tol::Float64`: Tolerance for convergence of the optimization.
- `max_iter::Int64`: Maximum number of iterations allowed. If not specified, defaults to n.
- `verbose::Int64`: Verbosity of output (0, 1, or 2).
- `inverse_method::String`: Least squares solver method, see [`dolsq`](@ref) for options.

# Returns
- A tuple containing:
    - `x`: Vector of size (n,) that minimizes the least squares residual.
    - `cost`: The final cost function value of the minimized sum of squared residuals.

# Notes
- Reference:  https://www.stat.berkeley.edu/~stark/Preprints/bvls.pdf
"""
function bvls(A, b, x_lsq, lb, ub, tol::Float64, max_iter::Int64, verbose::Int64,
    inverse_method::String)

    n_iter = 0
    m, n = size(A)

    x = x_lsq
    on_bound = zeros(n)

    mask = x .<= lb
    x[mask] = lb[mask]
    on_bound[mask] .= -1

    mask = x .>= ub
    x[mask] = ub[mask]
    on_bound[mask] .= 1

    free_set = on_bound .== 0
    active_set = .!free_set
    free_set = (1:size(free_set)[1])[free_set.!=0]

    r = A * x - b
    cost = 0.5 * dot(r, r)
    initial_cost = cost
    g = A' * r

    cost_change = nothing
    step_norm = nothing
    iteration = 0

    while size(free_set)[1] > 0
        if verbose == 2
            optimality = compute_kkt_optimality(g, on_bound)
        end

        iteration += 1
        x_free_old = x[free_set]

        A_free = A[:, free_set]
        b_free = b - A * (x .* active_set)
        z = dolsq(A_free, b_free, method=inverse_method)

        lbv = z .< lb[free_set]
        ubv = z .> ub[free_set]

        v = lbv .| ubv

        if any(lbv)
            ind = free_set[lbv]
            x[ind] = lb[ind]
            active_set[ind] .= true
            on_bound[ind] .= -1
        end

        if any(ubv)
            ind = free_set[ubv]
            x[ind] = ub[ind]
            active_set[ind] .= true
            on_bound[ind] .= 1
        end

        ind = free_set[.!v]
        x[ind] = z[.!v]

        r = A * x - b
        cost_new = 0.5 * dot(r, r)
        cost_change = cost - cost_new
        cost = cost_new
        g = A' * r
        step_norm = sum((x[free_set] .- x_free_old) .^ 2)

        if any(v)
            free_set = free_set[.!v]
        else
            break
        end
    end

    if isnothing(max_iter)
        max_iter = n
    end
    max_iter += iteration

    termination_status = nothing

    optimality = compute_kkt_optimality(g, on_bound)
    for iteration in iteration:max_iter
        if optimality < tol
            termination_status = 1
        end

        if !isnothing(termination_status)
            break
        end

        move_to_free = argmax(g .* on_bound)
        on_bound[move_to_free] = 0

        x_free = copy(x)
        x_free_old = copy(x)
        while true

            free_set = on_bound .== 0
            sum(free_set)
            active_set = .!free_set
            free_set = (1:size(free_set)[1])[free_set.!=0]

            x_free = x[free_set]
            x_free_old = copy(x_free)
            lb_free = lb[free_set]
            ub_free = ub[free_set]

            A_free = A[:, free_set]
            b_free = b - A * (x .* active_set)
            z = dolsq(A_free, b_free, method=inverse_method)

            lbv = (1:size(free_set)[1])[z.<lb_free]
            ubv = (1:size(free_set)[1])[z.>ub_free]
            v = cat(lbv, ubv, dims=1)

            if size(v)[1] > 0
                alphas = cat(
                    lb_free[lbv] - x_free[lbv],
                    ub_free[ubv] - x_free[ubv],
                    dims=1
                ) ./ (z[v] - x_free[v])

                i = argmin(alphas)
                i_free = v[i]
                alpha = alphas[i]

                x_free .*= (1 .- alpha)
                x_free .+= (alpha .* z)
                x[free_set] = x_free

                vsize = size(lbv)
                if i <= size(lbv)[1]
                    on_bound[free_set[i_free]] = -1
                else
                    on_bound[free_set[i_free]] = 1
                end
            else
                x_free = z
                x[free_set] = x_free
                @goto start
            end
        end #while
        @label start
        step_norm = sum((x_free .- x_free_old) .^ 2)

        r = A * x - b
        cost_new = 0.5 * dot(r, r)
        cost_change = cost - cost_new

        combo = tol * cost
        if cost_change < tol * cost
            termination_status = 2
        end
        cost = cost_new

        g = A' * r
        optimality = compute_kkt_optimality(g, on_bound)
    end #iteration

    if isnothing(termination_status)
        termination_status = 0
    end

    x[x.<1e-5] .= 0
    return x, cost
end

"""
    compute_kkt_optimality(g::Vector{Float64}, on_bound::Vector)

Computes the Karush-Kuhn-Tucker (KKT) optimality condition value for a given gradient
vector and a vector indicating which variables are on bounds.

- Returns a value representing the maximum KKT condition across all variables

# Arguments
- `g`: A vector of size (n,) representing the gradient of the objective function
  at the current point.
- `on_bound`: A vector of size (n,) indicating the status of each variable:
  - `-1` if the variable is at its lower bound,
  - `1` if the variable is at its upper bound,
  - `0` if the variable is free (not constrained).
"""
function compute_kkt_optimality(g::Vector{Float64}, on_bound::Vector)
    g_kkt = g .* on_bound
    free_set = on_bound .== 0
    g_kkt[free_set] = broadcast(abs, g[free_set])

    return maximum(g_kkt)
end
