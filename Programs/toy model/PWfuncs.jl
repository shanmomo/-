using LinearAlgebra
using DifferentialEquations.OrdinaryDiffEq
using Base.Threads

## spherical bessel functions
function contfrac_spbesselj(l::Int, x, maxiter)
    α = f = eps(Float64)
    bn = β = Δ = 0.
    @inbounds for n = 1:maxiter
        bn = (2l+1+2n)/x
        α = bn - 1/α #(α==0 ? α=eps(Float64) : );
        β = 1/(bn - β) #(α==0 ? α=eps(Float64) : );
        #β ^=(-1);
        Δ = α*β ; f*=Δ ;
        if abs(Δ-1)<eps(Float64)
            break
        end
    end
    return -f
end

function sphericalbessel_jy(l::Int, x; maxiter=10000)
    yo=y=-cos(x)/x; yp=(x*sin(x)+cos(x))/(x^2);
    jo=j=sin(x)/x; jp=(x*cos(x)-sin(x))/(x^2);
    @inbounds for n=0:l-1
        ##################y
        yo = y
        y = yo*(n/x) - yp
        yp = yo - y*(n+2)/x
        ##################j
        jo = j
        if n<x
            j = jo*(n/x) - jp
        else
            j *= contfrac_spbesselj(n,x,maxiter)
        end
        jp = jo - j*(n+2)/x
    end
    return j,y,jp,yp
end

## Lagrange Polynomials
@inline function PLv(l::Int, x::T) where T
    P = zeros(T,l+1)
    P[1]=1. ; P[2]=x
    for n=3:l+1
        @inbounds P[n]=((2n-3)*x*P[n-1]-(n-2)*P[n-2])/(n-1)
    end
    return P
end

@inline function QLv(l::Int, x::T) where T
    Q = zeros(T,l+1)
    Q[1]=log((1+x)/(1-x))/2; Q[2]=x*Q[1]-1
    for n=3:l+1
        @inbounds Q[n]=((2n-3)*x*Q[n-1]-(n-2)*Q[n-2])/(n-1)
    end
    return Q
end

## Initial Condition for Radial Function
function ul_ini!(l, ρmax, ρ0tolp1)
    u0 = ComplexF64[ρ0tolp1, 1.]
    ρspan = (ρ0tolp1*(l+1), ρmax)
    return u0, ρspan
end

## The Radial Differential Equation of ul
function ul_prob!(du, u, paras, ρ)
    l, k, Vn = paras
    @inbounds du[1] = u[2]
    @inbounds du[2] = ((l^2 + l) / (ρ^2) + Vn(ρ/k) - 1.0) * u[1]
end

## external解析解
function Hl(l, ρ)
    j,y,jp,yp = sphericalbessel_jy(l,ρ)
    H = ρ * complex(-y,j)    #
    Hp = complex(-y,j) + ρ * complex(-yp,jp)
    return H, conj(H), Hp, conj(Hp)
end
function ul_ext(l,ρ,Sl)
    H = Hl(l,ρ)
    u = 1im/2 *(H[2]-Sl*H[1])
    # up = 1im/2 *(H[4]-Sl*H[3])
    return u
end

## Match the Boundary Condition
function ul_match!(l, sol)   # Boundary Condition & Sl & Modify ul
    # mc = length(sol.t)
    # while (abs(sol[mc][1] / sol[mc][2]) |> (x -> (x + x^-1))) > 4.0
    #     mc = mc-1
    # end       # 避免匹配边界为节点的情况
    H = Hl(l, sol.t[end])
    u, up = sol.u[end]
    # equation: A*u = ρ( h2 + Sl*h1 )
    B, Sl = [-1im/2 *(H[2]*H[3]-H[1]*H[4]) ; H[2]*up-H[4]*u ] / (H[1]*up-H[3]*u)
    sol.u .*= B
    sol.k .*= B
    return Sl
end

## Core of The Program
function scatt(lmax, ρmax, k, Vn;
                Sltol=1e-7, ρstep=0.4, ρ0tolp1=0.3, l_adaptive=true)
    ls = []
    Sls = []
    sols = []
    for l = 0:lmax
        push!(ls, l)
        u0, ρspan = ul_ini!(l, ρmax, ρ0tolp1)
        ul_sol = solve(
            ODEProblem(ul_prob!, u0, ρspan, [l, k, Vn]),
            TsitPap8();
            dt = ρstep,
            adaptive = false
        )
        Sl = ul_match!(l, ul_sol)
        push!(Sls, Sl)
        push!(sols, ul_sol)
        if (abs(1-Sl) < Sltol) & (l_adaptive)
            break
        end
    end
    return Int64.(ls), ComplexF64.(Sls), sols
end


## Scattering Amplitudes and Cross Sections

function σ(θ, k, ls, Sls, sols)
    lm = ls[end]
    χl = [(2l+1) .*( Sls[l+1] - 1. ) for l in 0:lm ]
    P = PLv(lm,cos(θ))
    Q = QLv(lm,cos(θ))
    Q̃ = @.(P/2-1im/π*Q)
    f = sum(χl .* P) / (2im *k)
    fF = sum(χl .* Q̃) / (2im *k)
    fN = sum(χl .* conj.(Q̃)) / (2im *k)
    # # σR = (η/(2*k*sin(θ/2)^2))^2
    # fl = sum(χl[1:20].* P[1:20]) / (2im *k)
    # fg = sum(χl[21:end].* P[21:end]) / (2im *k)
    return abs(f)^2, abs(fF)^2, abs(fN)^2
    # return  abs(f)^2, abs(fl)^2, abs(fg)^2
end

## Interpolation Functions

# Radial Partial Wave Function 任意r位置
function ul(l, ρ, Sl, sol)
    if ρ < sol.t[1]
        return (ρ/sol.t[1])^(l+1) * sol[1][1]
    elseif ρ < ρmax
        return sol(ρ)[1]
    else
        return ul_ext(l, ρ, Sl)
    end
end

# return scattering wave function Ψsc(r,θ)*kr
function σ_r(θ, r, k, Sls, sols, lmax)
        χl = [(2l+1) .*(1im.^l) .*( ul( l, k*r, Sls[l+1], sols[l+1]) - k*r* sphericalbessel_jy(l,k*r)[1] ) for l in 0:lmax ]
        P = PLv(lmax ,cos(θ))
        return abs.(sum(χl.*P)/k).^2
end
