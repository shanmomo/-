include("PWfuncs.jl") # 过程函数调用

## 常数定义和参数输入
AMU = 931.5 # MeV
ħc = 197.32 # MeV⋅fm

include("PWinput.jl")

## 分波法求解
ls, Sls, sols =
    scatt(lmax, ρmax, k, V_Nuclear;
        Sltol=1e-7, ρstep=0.4, ρ0tolp1=0.3, l_adaptive=true)
