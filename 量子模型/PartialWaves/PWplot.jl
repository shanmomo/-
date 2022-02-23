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
# savefig(p0,"不同阶波函数.svg")

# #coulomb_F(ρ::Float64,l::Int,η)::Float64
# ρ50 = sols[11].t
# CF50 = map( ρ -> coulomb_F(ρ,10,η) , ρ50)
# plot!(ρ50./k,CF50,labels="coulomb function",lw=1)

## 散射振幅和散射截面 Ref:PWfuns
θgrid   = collect( LinRange(0.001, pi-0.001, 1000) )
σ_total = abs.(ft.(θgrid) ).^2
σ_far   = abs.(fnF.(θgrid)).^2
σ_near  = abs.(fcN.(θgrid) .+ fnN.(θgrid)).^2
σ_Ru = σR.(θgrid)
# begin
#     p1=plot(size = (800,600), minorgrid = true,
#     titlefontsize=18,
#     legendfontsize=12,
#     guidefontsize=14,
#     tickfontsize=10)
#     title!(L"^{16}\mathrm{O}+^{12}\mathrm{C}~~\mathrm{at}~~ E_{\mathrm{lab}}=132~\mathrm{MeV}")
#     xaxis!(L"\theta / \mathrm{deg}",(0,180))
#     yaxis!(L"\sigma/\sigma_{R}", (2e-6,2),:log)
#     plot!(rad2deg.(θgrid), σ_total./σ_Ru,
#         labels="Total",lw=2)
#     plot!(rad2deg.(θgrid), σ_far./σ_Ru,
#         labels="Far", ls=:dash)
#     plot!(rad2deg.(θgrid), σ_near./σ_Ru,
#         labels="Near",ls=:dash)
# end

begin
    p1=plot(size = (800,600), minorgrid = true,
    titlefontsize=18,
    legendfontsize=12,
    guidefontsize=14,
    tickfontsize=10)
    title!(L"\alpha+^{90}\mathrm{Zr}~~\mathrm{at}~~ E_{\mathrm{lab}}=79.5~\mathrm{MeV}")
    xaxis!(L"\theta / \mathrm{deg}",(0,180))
    yaxis!(L"\sigma/\sigma_{R}", (2e-8,2),:log)
    plot!(rad2deg.(θgrid), σ_total./σ_Ru,
        labels="Total",lw=2)
    plot!(rad2deg.(θgrid), σ_far./σ_Ru,
        labels="Far", ls=:dash)
    plot!(rad2deg.(θgrid), σ_near./σ_Ru,
        labels="Near",ls=:dash)
end

savefig(p1,"α_Zr_scatter.svg")

plt_l = 60
ψsc = ψsc_xy(k, η, plt_l, sols)
ψi = ψi_xy(k, η)

xgrid = collect(-10:0.05:25)
ygrid = collect(0.01:0.05:20) #y只用上半平面就可以
xymesh = [[x,y] for y in ygrid, x in xgrid]
rmesh = map(p->hypot(p...),xymesh)
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

#logscaling(x) = abs.(x).^2  |> x->log.(x .+ 1e-1)

symmfill(M) = vcat(reverse(M,dims=1),M)
yplt=vcat(-reverse(ygrid;dims=1),ygrid)
xplt=xgrid
symmfill(ψsc_mesh)

plt_zdata = abs.(ψsc_mesh./ rmesh).^2 |> symmfill
p2 = heatmap(size = (500,500),
    xplt,  # x轴
    yplt,  # y
    plt_zdata,  # z值
    xaxis = (L"x(\mathrm{fm})", (xplt[1], xplt[end])),
    yaxis = (L"y(\mathrm{fm})", (yplt[1], yplt[end])),
    #caxis = (:log),
    #clim  = (0, 20),  #default clim
    c = :thermal,
    colorbar_title = (L"|\psi_{N}|^2"),
    aspect_ratio = :equal,
    title = ("Wave Function Visualization")
)

savefig(p2,"O_C_scatter_psi_sc.svg")
