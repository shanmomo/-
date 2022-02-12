using LinearAlgebra
using DifferentialEquations.OrdinaryDiffEq
using SpecialFunctions
using HypergeometricFunctions
using GSL

## Special Functions

function σc_cal(l,η)
    return imag(loggamma(complex(1+l,η)))
end

function coulomb_H(ρ::Float64,l::Int,η)
    FG = map(x->x.val,sf_coulomb_wave_FG_e(η,ρ,l,0,0.,0.))
    return complex(FG[3],FG[1]), complex(FG[4],FG[2])
end

function coulomb_F(ρ::Float64,l::Int,η)::Float64
    if (ρ<(l*0.3)) & (l>20) #avoid over/underflow
        return 0.0
    else
        F=sf_coulomb_wave_FG_e(η,ρ,l,0,0.,0.)[1].val
        return Float64(F*(~isnan(F)))
    end
end

function CFarray(ρgrid,l::Int,η)
    return [sf_coulomb_wave_FG_e(η,ρ,l,0,0.,0.)[1].val for ρ in ρgrid]
end

function coulomb_G(ρ::Float64,l::Int,η)
    return sf_coulomb_wave_FG_e(η,ρ,l,0,0.,0.)[3].val
end

@inline function PLv(x::T, l::Int) where T
    P = zeros(T,l+1)
    P[1]=1. ; P[2]=x
    for n=3:l+1
        @inbounds P[n]=((2n-3)*x*P[n-1]-(n-2)*P[n-2])/(n-1)
    end
    return P
end

@inline function QLv(x::T, l::Int) where T
    Q = zeros(T,l+1)
    Q[1]=log((1+x)/(1-x))/2; Q[2]=x*Q[1]-1
    for n=3:l+1
        @inbounds Q[n]=((2n-3)*x*Q[n-1]-(n-2)*Q[n-2])/(n-1)
    end
    return Q
end
