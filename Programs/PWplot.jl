using Base.Threads # 多线程支持
using Plots   # 绘图
using LaTeXStrings  # 绘图LaTeX支持
gr() # 绘图Backend
theme(:default) # 绘图风格

## 绘制不同阶波函数
p0=plot(size = (800,600), minorgrid = true)
title!("Radial wave function real parts")
function plotureal(l)
    plot!(sols[l+1].t ./ k,real.(sols[l+1][1,:]),
        labels="l=$l")
end
plotureal(10)
plotureal(30)
plotureal(50)
xaxis!(L"r(\mathrm{fm})",(0,25))

# #coulomb_F(ρ::Float64,l::Int,η)::Float64
# ρ50 = sols[11].t
# CF50 = map( ρ -> coulomb_F(ρ,10,η) , ρ50)
# plot!(ρ50./k,CF50,labels="coulomb function",lw=1)

## 散射振幅和散射截面 Ref:PWfuns
begin
    θgrid   = collect( LinRange(0.001, pi-0.001, 1000) )
    σ_total = abs.(ft.(θgrid) ).^2
    σ_far   = abs.(fnF.(θgrid)).^2
    σ_near  = abs.(fcN.(θgrid) .+ fnN.(θgrid)).^2
    σ_Ru = σR.(θgrid)
end

begin
    p1=plot(size = (800,600), minorgrid = true)
    title!(L"^{16}\mathrm{O}+^{12}\mathrm{C}~~\mathrm{at}~~ E_{\mathrm{lab}}=132~\mathrm{MeV}")
    xaxis!(L"\theta / \mathrm{deg}",(0,180))
    yaxis!(L"\sigma/\sigma_{R}", (2e-6,2),:log)
    plot!(rad2deg.(θgrid), σ_total./σ_Ru,
        labels="Total",lw=2)
    plot!(rad2deg.(θgrid), σ_far./σ_Ru,
        labels="Far", ls=:dash)
    plot!(rad2deg.(θgrid), σ_near./σ_Ru,
        labels="Near",ls=:dash)
end
savefig(p1,"Fig\\O_C_scatter.svg")

plt_l = 80
ψsc = ψsc_xy(k, η, plt_l, sols)
ψi = ψi_xy(k, η)

begin
    xgrid = collect(-15:0.05:45)
    ygrid = collect(0.01:0.05:20) #y只用上半平面就可以
    xymesh = [[x,y] for y in ygrid, x in xgrid]
    rmesh = map(p->hypot(p...),xymesh)
end

ψsc_mesh = zeros(ComplexF64,size(xymesh))
@threads for ii = 1:length(xymesh)
    x,y = xymesh[ii]
    @inbounds ψsc_mesh[ii] = ψsc(x,y)
end

ψi_mesh = zeros(ComplexF64,size(xymesh))
@threads for ii = 1:length(xymesh)
    x,y = xymesh[ii]
    @inbounds ψi_mesh[ii] = ψi(x,y)
end

ψ_mesh = ψi_mesh .+ (ψsc_mesh ./ rmesh ./k )
#logscaling(x) = abs.(x).^2  |> x->log.(x .+ 1e-1)

symmfill(M) = vcat(reverse(M,dims=1),M)
yplt=vcat(-reverse(ygrid;dims=1),ygrid)
xplt=xgrid

theme(:lime)
plt_zdata = real.(ψ_mesh) |> symmfill
p2 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x", (xplt[1], xplt[end])),
    yaxis = (L"y", (yplt[1], yplt[end])),
    clim  = (-2, 2),  #default clim
    c = :RdBu_9,
    colorbar_title = (L"\mathrm{Re}(\psi)"),
    aspect_ratio = :equal
)
savefig(p2,"Fig\\O_C_scatter_RePsi.svg")

plt_zdata = abs.(ψ_mesh).^2 |> symmfill
p2 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x", (xplt[1], xplt[end])),
    yaxis = (L"y", (yplt[1], yplt[end])),
    clim  = (-0, 3),  #default clim
    c = :thermal,
    colorbar_title = (L"|\psi|^2"),
    aspect_ratio = :equal
)
savefig(p2,"Fig\\O_C_scatter_abs2Psi.svg")

plt_zdata = real.(ψsc_mesh./rmesh./k) |> symmfill
p2 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x", (xplt[1], xplt[end])),
    yaxis = (L"y", (yplt[1], yplt[end])),
    clim  = (-2, 2),  #default clim
    c = :RdBu,
    colorbar_title = (L"\mathrm{Re}(\psi_{\mathrm{sc}})"),
    aspect_ratio = :equal
)
savefig(p2,"Fig\\O_C_scatter_RePsisc.svg")

plt_zdata = abs.(ψsc_mesh./rmesh./k).^2 |> symmfill
p2 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x", (xplt[1], xplt[end])),
    yaxis = (L"y", (yplt[1], yplt[end])),
    clim  = (-0, 2),  #default clim
    c = :thermal,
    colorbar_title = (L"|\psi_{\mathrm{sc}}|^2"),
    aspect_ratio = :equal
)
savefig(p2,"Fig\\O_C_scatter_abs2Psisc.svg")



#
#
# # julia画图示例
# 1
#
# using Plots
# Plots.default()
# gr() # 设置backend
#
# # 设置不同的图窗主题，可以选择.运行其中一行就行
# theme(:default)
# theme(:ggplot2)
# theme(:sand)
# theme(:dao)
# theme(:juno)
# theme(:vibrant)
#
#
# θ = LinRange(0,10π,1000)
# x = cos.(θ)
# y = sin.(θ)
# z = 0.1 .*θ
# plot(x,y,z,line_z = z)
#
#
#
#
# x = LinRange(0,2π,200)
# y = LinRange(0,2π,200)
# z = [sin(x)*sin(y) for x in x, y in y]
#
#
# surface(x,y,z,camera = (15,65))
# contour(x,y,z)
#
# Plots.
# heatmap(x,y,z,
#     aspect_ratio = 1, #x,y轴放缩比例，默认是自动
#     # c = :RdBu, # colormap的选项，可以查手册。默认是根据主题变的
#     xlims = (0,2π), # 坐标轴范围
#     ylims = (0,2π),
#     )
#
