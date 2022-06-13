using LinearAlgebra
using DifferentialEquations.OrdinaryDiffEq

## 运动轨迹微分方程
# u[1]=r; u[2]=θ; u[3]=ṙ; u[4]=θ̇; u[5]=real(ϕ); u[6]=imag(ϕ); u[7]=无吸收ϕ
function CM( du, u, p, t)
    μ, k0, Ecm = p
    dϕ = k0*sqrt(complex(1-cenV(u[1])/Ecm)) * hypot(u[3],u[1]*u[4])
    @inbounds begin
    du[1] = u[3]
    du[2] = u[4]
    du[3] = u[4]^2 * u[1] - dcenV(u[1])/μ
    du[4] = -2*u[3] *u[4] /u[1]
    du[5] = real(dϕ)
    du[6] = imag(dϕ)
    # du[7] = k0*sqrt(1-realcenV(u[1])/Ecm) * hypot(u[3],u[1]*u[4])
    end
end

function solveCMproblem(u0,tspan,p)
    solve(ODEProblem(CM,u0,tspan,p),TsitPap8();
        dtmax = 0.5/(k0*v0),
        reltol = 1e-8, abstol = 1e-8,
        adaptive = true)
end

function theta(sol)
    u = sol[end]
    return u[2] + atan(u[1]*u[4]/u[3])
end

function phase(sol, k0)
    u = sol[end]
    return mod2pi(u[5] - k0 * u[1] * u[3]/sqrt(u[3]^2+(u[1]*u[4])^2))
end

# ## 统计散射截面
#
# N = 500
# dθ = π/N
# θlist = LinRange(dθ,π,N)
# σlist = zeros(N)
#
# function
#
# function crosssection(sols,b,Δb,)
#     u = [sol.u for sol in sols]
#     θ = [π .- abs.(u[2][end] .-π) for u in u]
#     for i in 1:n
#         j = div(θ[i],dθ)+1
#         σ[j] =
#
# end
