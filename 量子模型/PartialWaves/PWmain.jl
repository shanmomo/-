include("PWspfuncs.jl")
include("PWfuncs.jl") # 过程函数调用

## 常数定义和参数输入
AMU = 931.5 # MeV
e²4πϵ = 1.44 # MeV⋅fm
ħc = 197.32 # MeV⋅fm
include("PWinput.jl")

## 分波法求解
ls, σcs, Sls, sols =
    scatt(ρmax, lmax, η, k, V_Coulomb, V_Nuclear;
        relρ0=relρ0, Sltol=Sltol, tstep=ρstep, l_adaptive=l_adaptive)

## 散射截面
fn,fnF,fnN,fc,fcF,fcN,ft,σR = fcal(k,Sls,σcs,η)
