
## --------PHYSICAL_REVIEW C, VOLUME 62, 044601 ———————

# 重离子类型和能量
Mp,Zp = 16.0,8.0
Mt,Zt = 12.0,6.0
Elab = 132.0 # MeV

## 参数计算

Ecm = Elab*Mt/(Mp+Mt)       # 质心系能量 MeV
μ = Mp*Mt/(Mp+Mt)           # 折合质量
unitlessR = cbrt(Mp)+cbrt(Mt)
k0 = sqrt(2*μ*AMU*Ecm)/ħc   # 初波矢 /fm
v0 = sqrt(2*Ecm/μ)          # 初始速度  √MeV

# b = sqrt(40*81*ħc^2/Ecm/(μ*AMU))    # fm


## 势函数参数和定义

# 核势    Wood-Saxon: [depth(MeV), width_factor(fm), skewness(fm)]
WS1_V = [282.2, 0.586, 0.978]
WS1_W = [13.86, 1.183, 0.656]

dVn = WS1_V[1]              # Mev
rVn = WS1_V[2] * unitlessR  # fm
aVn = WS1_V[3]              # fm

dWn = WS1_W[1]              # MeV
rWn = WS1_W[2] * unitlessR  # fm
aWn = WS1_W[3]              # fm

# 中心势 MeV
cenV(r) = Complex(-dVn/(exp((r-rVn)/aVn)+1.), -dWn/(exp((r-rWn)/aWn)+1.))
realcenV(r) = -dVn/(exp((r-rVn)/aVn)+1.)

using FiniteDifferences
dcenV(r) = central_fdm(5, 1)(realcenV,r)

## 入射粒子分布

x0 = -100.0      # 初始位置
b0 = 15.0       # 入射范围 fm
n = 1000          # 粒子数
blist = collect(LinRange(0,b0,n))
Δb = exp.(-2 *(blist/b0 .- 1.).^2)
Δb = Δb / sum(Δb) * b0
b = zeros(n); b[1] = Δb[1]
for i in 2:n
    b[i] = b[i-1] + Δb[i]
end
wb0 = zeros(n)
for i in 1:n
    wb0[i] = b[i] * Δb[i]
end
wb0 = wb0 / sum(wb0)    # 权重分布

L2 = @.( 2*Ecm*μ*AMU / ħc^2 *b^2 )
l = sqrt.(L2)

## 运动时间
t0 = 200. /v0    # fm/√MeV
