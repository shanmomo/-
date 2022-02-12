using LinearAlgebra
using DifferentialEquations.OrdinaryDiffEq
using SpecialFunctions
using HypergeometricFunctions
using GSL

## Initial Condition for Radial Function
function ul_ini!(ρmax::Float64, l::Int; ρ0tolp1 = 0.3)
    u0 = ComplexF64[ρ0tolp1, 1.]
    ρspan = (ρ0tolp1*(l+1), ρmax)
    return u0, ρspan
end

## The Radial Differential Equation of ul
function ul_prob!(du, u, paras, ρ)
    l, η, k, Vc, Vn = paras
    @inbounds du[1] = u[2]
    @inbounds du[2] = ((l^2 + l) / (ρ^2) + Vc(ρ) + Vn(ρ/k) - 1.0) * u[1]
end

## Match the Boundary Condition
function ul_match!(sol, l, η)    # Boundary Condition & Sl & Modify ul
    # mc = length(sol.t)
    # while (abs(sol[mc][1] / sol[mc][2]) |> (x -> (x + x^-1))) > 4.0
    #     mc = mc-1
    # end                      # Avoid u->0 or u'->0
    H, Hp = coulomb_H(sol.t[end], l, η)
    ul, ulp = sol.u[end]
    # equation: A*u = ρ( h2 + Sl*h1 )
    A, Sl = [-2im * ul H; -2im * ulp Hp] \ conj.([H,Hp])
    sol.u .*= A
    sol.k .*= A
    return Sl
end

## Core of The Program
function scatt(ρmax, lmax, η, k, Vc, Vn ; relρ0=0.3, Sltol=1e-8, tstep=0.4, l_adaptive)
    ls = []
    σcs = []
    Sls = []
    sols = []
    for l = 0:lmax
        push!(ls,l)
        push!(σcs,σc_cal(l,η))
        u0, ρspan = ul_ini!(ρmax, l; ρ0tolp1 = relρ0)
        ul_sol = solve(
            ODEProblem(ul_prob!, u0, ρspan, [l, η, k, Vc, Vn]),
            TsitPap8();
            dt = tstep,
            adaptive = false
        )
        Sl = ul_match!(ul_sol, l, η)
        push!(Sls,Sl)
        push!(sols, ul_sol)
        if (abs(1-Sl)<Sltol) & (l_adaptive)
            break
        end
    end
    return ls, σcs, Sls, sols
end

## Scattering Amplitudes and Cross Sections

function fcal(k, Sls, σcs, η)
    ls = collect(0:length(Sls)-1)
    lm = length(Sls) - 1
    ntemp = ((2im * k)^-1) .* (2ls .+ 1) .* exp.(2im.*σcs) .* (Sls .- 1.0)
    # ctemp = ((2im * k)^-1) .* (2ls .+ 1) .* (exp.(2im.*σcs) .- 1)
    # ttemp = ((2im * k)^-1) .* (2ls .+ 1) .* (Sls.*exp.(2im.*σcs) .- 1)
    fnlT(θ) = ntemp .* PLv(cos(θ),lm)
    fnlF(θ) = ntemp.*((1/2).*PLv(cos(θ),lm).-(1im/π).*QLv(cos(θ),lm))
    fnlN(θ) = ntemp.*((1/2).*PLv(cos(θ),lm).+(1im/π).*QLv(cos(θ),lm))
    # fclT(θ) = ctemp .* PLv(cos(θ),lm)
    # fclF(θ) = ctemp.*((1/2).*PLv(cos(θ),lm).-(1im/π).*QLv(cos(θ),lm))
    # fclN(θ) = ctemp.*((1/2).*PLv(cos(θ),lm).+(1im/π).*QLv(cos(θ),lm))
    # ftlT(θ) = ttemp .* PLv(cos(θ),lm)
    # ftlF(θ) = ttemp.*((1/2).*PLv(cos(θ),lm).-(1im/π).*QLv(cos(θ),lm))
    # ftlN(θ) = ttemp.*((1/2).*PLv(cos(θ),lm).+(1im/π).*QLv(cos(θ),lm))

    # Rutherford cross section
    σR(θ) = (η/(2*k*sin(θ/2)^2))^2
    # Nuclear scattering amplitude
    fn(θ; sumrange=1:lm+1) = sum(fnlT(θ)[sumrange])
    fnF(θ; sumrange=1:lm+1) = sum(fnlF(θ)[sumrange])
    fnN(θ; sumrange=1:lm+1) = sum(fnlN(θ)[sumrange])
    # Coulomb scattering amplitude
    fc(θ) = -η/(2k*sin(θ/2)^2)*exp(-1im*η*log(sin(θ/2)^2)+2im*σcs[1])
    S_term(θ) = _₂F₁(1,complex(1,η),complex(2,η),sin(θ/2)^2)
    fcN(θ) = fc(θ) * ( (1-exp(-2*π*η)).^(-1) -
                    (1im/2π) * sin(θ/2)^complex(2,2η) * S_term(θ) )
    fcF(θ) = fc(θ) * (-exp(-2*π*η)) * ( (1-exp(-2*π*η)).^(-1) +
                    (1im/2π) * sin(θ/2)^complex(2,2η) * S_term(θ) )
    # Total amplitude
    ft(θ) = fn(θ) + fc(θ)

    return fn,fnF,fnN,fc,fcF,fcN,ft,σR
end

## Interpolation Functions

# Radial Partial Wave Function
function ul(sol,ρ,l)
    if ρ < sol.t[1]
        return (ρ/sol.t[1])^(l+1) * sol[1][1]
    else
        return sol(ρ)[1]
    end
end

# return scattering wave function Ψsc_xy(x,y)
function ψsc_xy(k, η, lmax, sols)
    function ψsc(x,y)::ComplexF64
        r = hypot(x,y)
        c = x/(r+eps(Float64))
        χl = [(2l+1).*(1im.^l).*(ul(sols[l+1],k*r,l) - coulomb_F(k*r,l,η)  ) for l in 0:lmax]
        P = PLv(c,lmax)
        return dot(χl,P)
    end
    return ψsc
end

function ψi_xy(k, η)
    Amp=gamma(complex(1,η))*exp(-π*η/2)
    function ψi(x,y)::ComplexF64
        return Amp * exp(1im*k*x) * pFq([-1im*η,],[1,],1im*k*(hypot(x,y)-x))
    end
    return ψi
end

# # Radial part in partial-wave expansion, not used
# function ψsc_Rl(k, η, lm, sols, rgrid)
#     ρgrid = rgrid.*k
#     usc=zeros(ComplexF64,length(ρgrid),lm+1)
#     for l=0:lm
#         @inbounds usc[:,l+1] =
#             (2l+1).*(1im.^l).*(sol(ρgrid)[1,:].-CFarray(ρgrid,l,η))
#     end
#     return usc
# end

# # Angular part in partial-wave expansion, not used
# function ψsc_Θl(lm, θgrid)
#     Psc=zeros(lm+1,length(θgrid))
#     for ii=1:length(θgrid)
#         @inbounds Psc[:,ii]=PLv(cos(θgrid[ii]),lm)
#     end
#     return Psc
# end
