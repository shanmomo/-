## ------- Nuclear Physics A291 (1977) 93-125 ———————

# 重离子类型和能量
Mp,Zp = 4.0,2.0
Mt,Zt = 90.0,40.0
Elab = 79.5 # MeV

# 势参数
Rc_scale = 1.3  # 库伦势半径因子: fm
WSs = [
        [141.2, 1.225, 0.821, 18.49, 1.575, 0.565] # 79.5 MeV
        # [133.3, 1.235, 0.805, 19.63, 1.571, 0.562], # 99.5 MeV
        # [130.0, 1.231, 0.821, 20.03, 1.572, 0.568]  # 118 MeV
        ]   # 标准WS模型

# 求解范围
lmax = 200 #最高分波阶数
ρmax = 300. #最大无量纲半径，数值不应小于最高分波阶数

# 精度&步长
Sltol = 1e-8 # break, when |Sl-1|<Sltol
ρstep = 0.3 # unitless radius ρ step.
ρ0tolp1 = 0.3 # relative parameter for ρ_start
l_adaptive = true # if false, only break at l=lmax

## 参数计算
Ecm = Elab * Mt/(Mp+Mt)
μ = Mp*Mt/(Mp+Mt)*AMU
k = sqrt(2*Ecm*μ)/ħc
η = Zp*Zt*e²4πϵ*μ/(ħc^2*k)
unitlessR = cbrt(Mt)
rmax = ρmax/k

## 势函数计算

# 库伦势, 无量纲半径为参数。也可以使用double-fold类型，注意连续性
Rc = Rc_scale * unitlessR
ρVc = Rc * k
function Vc_1(ρ::Float64)::Float64
    if ρ > ρVc
        return 2*η/ρ #渐进行为
    else
        return 2*η/ρVc*(3/2-ρ^2/(2*ρVc^2))
    end
end

# 核势 Wood-Saxon: [depth(MeV), width_factor(fm), skewness(fm)]
struct WSparams
    dV::Float64
    rV::Float64
    aV::Float64
    dW::Float64
    rW::Float64
    aW::Float64
end

function parseWSparams(WS)  # 约化参数并写成WSparams类型
    return WSparams((WS.*repeat([Ecm^-1, unitlessR, 1],2))...)  # 放入结构体需解包
end

function genericWS(r::Float64, P::WSparams)::ComplexF64
    return ComplexF64(
    -P.dV/(exp((r-P.rV)/P.aV)+1.),
    -P.dW/(exp((r-P.rW)/P.aW)+1.)
    )
end

WSp = parseWSparams.(WSs)

## 主程序中具体选取的势能，只能有两项V_Coulomb,V_Nuclear
V_Coulomb(ρ) = Vc_1(ρ)
V_Nuclear(r) = sum(genericWS.(r,WSp))
