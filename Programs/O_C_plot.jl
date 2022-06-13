using Base.Threads # 多线程支持
using Plots   # 绘图
using LaTeXStrings  # 绘图LaTeX支持
gr() # 绘图Backend
theme(:default) # 绘图风格

## 绘制不同阶波函数
# p0=plot(size = (800,600), minorgrid = true)
# title!("Radial wave function real parts")
# function plotureal(l)
#     plot!(sols[l+1].t ./ k,real.(sols[l+1][1,:]),
#         labels="l=$l")
# end
# plotureal(10)
# plotureal(30)
# plotureal(50)
# xaxis!(L"r(\mathrm{fm})",(0,25))

# #coulomb_F(ρ::Float64,l::Int,η)::Float64
# ρ50 = sols[11].t
# CF50 = map( ρ -> coulomb_F(ρ,10,η) , ρ50)
# plot!(ρ50./k,CF50,labels="coulomb function",lw=1)

## 散射振幅和散射截面 Ref:PWfuns
begin
    θgrid   = collect( LinRange(0.001, pi-0.001, 1000) )
    σ_total = abs.(ft.(θgrid) ).^2
    σ_far   = abs.(fcF.(θgrid) .+ fnF.(θgrid)).^2
    σ_near  = abs.(fcN.(θgrid) .+ fnN.(θgrid)).^2
    σ_c_far = abs.(fcF.(θgrid) ).^2
    σ_c_near= abs.(fcN.(θgrid) ).^2
    σ_n_far = abs.(fnF.(θgrid) ).^2
    σ_n_near= abs.(fnN.(θgrid) ).^2
    σ_n = abs.(fn.(θgrid) ).^2
    σ_Ru = σR.(θgrid)
end

# 与FRESCO计算散射截面对比
begin
    p0=plot(size = (800,600), minorgrid = true,
            legendfontsize = 14, tickfontsize = 12, guidefontsize = 18,
            framestyle = :box, margins = 3mm)
    xaxis!(L"\Theta_{\mathrm{c.m.}} (\mathrm{deg})",(0,180))
    yaxis!(L"\sigma/\sigma_{R}", (2e-5,2),:log)
    plot!(rad2deg.(θgrid), σ_total./σ_Ru,
        labels="This program",lw=3)
    plot!(θf,σf,labels="FRESCO",ls=:solid,lw=2,linealpha=0.8)
    annotate!(60,1e0,text(L"^{16}\mathrm{O}+^{12}\mathrm{C}",24))
    annotate!(70,0.2e0,text(L"E_{\mathrm{lab}}=132.0~\mathrm{MeV}",24))
end
savefig(p0,"Fig\\O_C_Comparison.svg")

#绘制近远端散射截面
begin
    p1=plot(size = (800,600), minorgrid = true,
            legendfontsize = 12, tickfontsize = 12, guidefontsize = 18,
            framestyle = :box, margins = 3mm)
    xaxis!(L"\Theta_{\mathrm{c.m.}} (\mathrm{deg})",(0,180))
    yaxis!(L"\sigma/\sigma_{R}", (2e-6,2),:log)
    plot!(rad2deg.(θgrid), σ_total./σ_Ru,
        labels="Total",lw=2)
    plot!(rad2deg.(θgrid), σ_far./σ_Ru,
        labels="Far", ls=:dash)
    plot!(rad2deg.(θgrid), σ_near./σ_Ru,
        labels="Near",ls=:dash)
    annotate!(60,1e0,text(L"^{16}\mathrm{O}+^{12}\mathrm{C}",24))
    annotate!(70,0.2e0,text(L"E_{\mathrm{lab}}=132.0~\mathrm{MeV}",24))
end
savefig(p1,"Fig\\O_C_Scatter.svg")


# begin
#     p1=plot(size = (800,600), minorgrid = true)
#     title!(L"^{16}\mathrm{O}+^{12}\mathrm{C}~~\mathrm{at}~~ E_{\mathrm{lab}}=132.0~\mathrm{MeV}")
#     xaxis!(L"\theta / \mathrm{deg}",(0,180))
#     yaxis!(L"\sigma", (2e-6,1e4),:log)
#     plot!(rad2deg.(θgrid), σ_n,
#         labels="Nuclear",lw=2)
#     plot!(rad2deg.(θgrid), σ_n_far./σ_Ru,
#             labels="Nuclear Far", ls=:dash)
#     plot!(rad2deg.(θgrid), σ_n_near./σ_Ru,
#             labels="Nuclear Near", ls=:dash)
# end


## 散射过程
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

theme(:default)
# plt_zdata = real.(ψ_mesh) |> symmfill
# p2 = heatmap(size = (600,360),
#     xplt,  # x轴
#     yplt,  # y
#     plt_zdata,  # z值
#     xaxis = (L"x", (xplt[1], xplt[end])),
#     yaxis = (L"y", (yplt[1], yplt[end])),
#     clim  = (-2, 2),  #default clim
#     c = :RdBu_9,
#     colorbar_title = (L"\mathrm{Re}(\psi)"),
#     aspect_ratio = :equal
# )
# savefig(p2,"Fig\\O_C_scatter_RePsi.svg")
#
# plt_zdata = abs.(ψ_mesh).^2 |> symmfill
# p2 = heatmap(size = (600,360),
#     xplt,  # x轴
#     yplt,  # y
#     plt_zdata,  # z值
#     xaxis = (L"x", (xplt[1], xplt[end])),
#     yaxis = (L"y", (yplt[1], yplt[end])),
#     clim  = (-0, 3),  #default clim
#     c = :thermal,
#     colorbar_title = (L"|\psi|^2"),
#     aspect_ratio = :equal
#     )
# savefig(p2,"Fig\\O_C_scatter_abs2Psi.svg")

# plt_zdata = real.(ψsc_mesh./rmesh./k) |> symmfill
# p2 = heatmap(size = (600,360),
#     xplt,  # x轴
#     yplt,  # y
#     plt_zdata,  # z值
#     xaxis = (L"x", (xplt[1], xplt[end])),
#     yaxis = (L"y", (yplt[1], yplt[end])),
#     clim  = (-2, 2),  #default clim
#     c = :RdBu,
#     colorbar_title = (L"\mathrm{Re}(\psi_{\mathrm{sc}})"),
#     aspect_ratio = :equal
# )
# savefig(p2,"Fig\\O_C_scatter_RePsisc.svg")

plt_zdata = abs.(ψsc_mesh./rmesh./k).^2 |> symmfill
# p2 = heatmap(size = (600,360),
#     xplt,  # x轴
#     yplt,  # y
#     plt_zdata,  # z值
#     xaxis = (L"x", (xplt[1], xplt[end])),
#     yaxis = (L"y", (yplt[1], yplt[end])),
#     clim  = (-0, 2),  #default clim
#     c = :thermal,
#     colorbar_title = (L"|\psi_{\mathrm{sc}}|^2"),
#     aspect_ratio = :equal
# )
# savefig(p2,"Fig\\O_C_scatter_abs2Psisc.svg")

ϕ = LinRange(0,2π,1000)
L = argmin(abs.(real.(Sls).-0.5))-1
R12 = (sqrt(L*(L+1)+η^2)+η) /k
Rc

p3 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"z(\mathrm{fm})", (xplt[1], xplt[end])),
    yaxis = (L"R(\mathrm{fm})", (yplt[1], yplt[end])),
    clim  = (-0, 2),  #default clim
    c = :thermal,
    colorbar_title = (L"|\psi_{\mathrm{N}}|^2"),
    aspect_ratio = :equal
    )
plot!( Rc*cos.(ϕ), Rc*sin.(ϕ), lw=2, label="strong absorption radius" )
savefig(p3,"Fig\\O_C_radius.svg")
