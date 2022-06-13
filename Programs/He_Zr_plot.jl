using Base.Threads # 多线程支持
using Plots   # 绘图
using LaTeXStrings  # 绘图LaTeX支持
gr() # 绘图Backend
theme(:default) # 绘图风格

## 散射振幅和散射截面 Ref:PWfuns
begin
    θgrid   = collect( LinRange(0.001, pi-0.001, 1000) )
    σ_total = abs.(ft.(θgrid) ).^2
    σ_far   = abs.(fnF.(θgrid)).^2
    σ_near  = abs.(fcN.(θgrid) .+ fnN.(θgrid)).^2
    σ_Ru = σR.(θgrid)
end

# 与FRESCO计算散射截面对比
begin
    p0=plot(size = (800,600), minorgrid = true,
            legendfontsize = 14, tickfontsize = 12, guidefontsize = 18,
            framestyle = :box, margins = 3mm)
    xaxis!(L"\Theta_{\mathrm{c.m.}} (\mathrm{deg})",(0,180))
    yaxis!(L"\sigma/\sigma_{R}", (2e-6,2),:log)
    plot!(rad2deg.(θgrid), σ_total./σ_Ru,
        labels="This program",lw=3)
    plot!(θf,σf,labels="FRESCO",ls=:solid,lw=2,linealpha=0.8)
    annotate!(60,1e0,text(L"\alpha+^{90}\mathrm{Zr}",24))
    annotate!(70,0.2e0,text(L"E_{\mathrm{lab}}=79.5~\mathrm{MeV}",24))
end
savefig(p0,"Fig\\α_Zr_Comparison.svg")

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
    annotate!(60,1e0,text(L"\alpha+^{90}\mathrm{Zr}",24))
    annotate!(70,0.2e0,text(L"E_{\mathrm{lab}}=79.5~\mathrm{MeV}",24))
end
savefig(p1,"Fig\\α_Zr_Scatter.svg")

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

theme(:lime)
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
# # savefig(p2,"Fig\\O_C_scatter_RePsi.svg")
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
# )
# # savefig(p2,"Fig\\O_C_scatter_abs2Psi.svg")
#
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
# # savefig(p2,"Fig\\O_C_scatter_RePsisc.svg")

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
# savefig(p2,"Fig\\α_Zr_scatter_abs2Psisc.svg")

ϕ = LinRange(0,2π,1000)
p3 = heatmap(size = (600,360),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x", (xplt[1], xplt[end])),
    yaxis = (L"y", (yplt[1], yplt[end])),
    clim  = (-0, 2),  #default clim
    c = :thermal,
    # colorbar_title = (L"|\psi_{\mathrm{sc}}|^2"),
    aspect_ratio = :equal
)
plot!( Rc*cos.(ϕ), Rc*sin.(ϕ), lw=1 , label="Nuclear Radius" )
savefig(p3,"Fig\\α_Zr_radius.svg")
