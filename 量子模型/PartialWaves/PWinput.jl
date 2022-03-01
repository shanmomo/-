## 分波法参数输入

## --------PHYSICAL_REVIEW C, VOLUME 62, 044601 ———————
# 重离子类型和能量
# Mp,Zp = 16.0,8.0
# Mt,Zt = 12.0,6.0
# Elab = 132.0 # MeV

Mp,Zp = 4.0,2.0
Mt,Zt = 90.0,40.0
Elab = 79.5 # MeV

# 求解范围
lmax = 200 #最高分波阶数
ρmax = 400. #最大无量纲半径，数值不应小于最高分波阶数

# 精度&步长
Sltol = 1e-8 # break, when |Sl-1|<Sltol
ρstep = 0.3 # unitless radius ρ step.
relρ0 = 0.3 # relative parameter for ρ_start
l_adaptive = true # if false, only break at l=lmax

## 参数计算
Ecm = Elab * Mt/(Mp+Mt)
μ = Mp*Mt/(Mp+Mt)*AMU
k = sqrt(2*Ecm*μ)/ħc
η = Zp*Zt*e²4πϵ*μ/(ħc^2*k)
unitlessR = cbrt(Mp) + cbrt(Mt)
rmax = ρmax/k

## 势函数参数和定义
# 库伦势, 无量纲半径为参数。也可以使用double-fold类型，注意连续性
# Rc_scale = 1.2 # 库伦势半径因子: fm
Rc_scale = 1.3 # 库伦势半径因子: fm
Rc = Rc_scale * unitlessR
ρVc = Rc * k

function Vc_1(ρ::Float64)::Float64
    if ρ > ρVc
        return 2*η/ρ #渐进行为
    else
        return 2*η/ρVc*(3/2-ρ^2/(2*ρVc^2))
    end
end

# 核势
# Wood-Saxon: [depth(MeV), width_factor(fm), skewness(fm)]
# WS1_V = [282.2, 0.586, 0.978]
# WS1_W = [13.86, 1.183, 0.656]

WS1_V = [141.2, 1.225, 0.821]
WS1_W = [18.49, 1.575, 0.565]

dVn = WS1_V[1] / Ecm # MeV/Ecm
rVn = WS1_V[2] * unitlessR #fm
aVn = WS1_V[3] #fm

dWn = WS1_W[1] / Ecm # MeV/Ecm
rWn = WS1_W[2] * unitlessR #fm
aWn = WS1_W[3] #fm

# dWD = WS1_V[1] / Ecm # MeV/Ecm
# rVD = WS1_V[2] * unitlessR #fm
# aVD = WS1_V[3] #fm

normal_WS1(r::Float64)::ComplexF64 = Complex(-dVn/(exp((r-rVn)/aVn)+1.), -dWn/(exp((r-rWn)/aWn)+1.))

## 主程序中具体选取的势能，只能有两项V_Coulomb,V_Nuclear
V_Coulomb(ρ) = Vc_1(ρ)
V_Nuclear(r) = normal_WS1(r)
