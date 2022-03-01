include("PWfuncs.jl") # 过程函数调用

## 常数定义和参数输入
AMU = 931.5 # MeV
e²4πϵ = 1.44 # MeV⋅fm
ħc = 197.32 # MeV⋅fm
# include("O_C_input.jl")
include("He_Zr_input.jl")

## 分波法求解
ls, σcs, Sls, sols =
    scatt(lmax, η, ρmax, k, V_Coulomb, V_Nuclear;
            Sltol=Sltol, ρstep=ρstep, ρ0tolp1=ρ0tolp1, l_adaptive=l_adaptive)

## 散射截面
fn,fnF,fnN,fc,fcF,fcN,ft,σR = fcal(η, k, ls, σcs, Sls)

## 与FRESCO计算标准结果对比
# include("O_C_fresco.jl")
include("He_Zr_fresco.jl")
