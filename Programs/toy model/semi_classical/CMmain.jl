
## 常数定义和参数输入
AMU = 931.5     # MeV
ħc = 197.32     # MeV⋅fm

## 调用
include("CMinput.jl")
include("CMfuncs.jl")

## 运动轨迹

p = [μ, k0, Ecm]
tspan = (0. , t0)
rilist = [ hypot(x0,b) for b in b]
θilist = [ atan(b/x0)+ π  for b in b]

# u0list = [ [ rilist[ii], θilist[ii], v0*cos(θilist[ii]),
#          -v0*sin(θilist[ii])/rilist[ii], 0, 0, 0] for ii in 1:n ]
u0list = [ [ rilist[ii], θilist[ii], v0*cos(θilist[ii]),
         -v0*sin(θilist[ii])/rilist[ii], 0, 0] for ii in 1:n ]

sols = [solveCMproblem(u0,tspan,p) for u0 in u0list]
