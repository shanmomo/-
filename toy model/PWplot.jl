## 绘制散射过程截面图
using Base.Threads  # 多线程支持
using Plots         # 绘图
using LaTeXStrings  # 绘图LaTeX支持
gr()            # 绘图Backend
theme(:default) # 绘图风格

## 散射截面
θlist = collect( LinRange(0.001, pi-0.001, 1000) )
σlist = [σ(θ, k, ls, Sls, sols) for θ in θlist]

begin
    p0=plot(size = (800,600), minorgrid = true)
    title!(L"^{16}\mathrm{O}+^{12}\mathrm{C}~~\mathrm{at}~~ E_{\mathrm{lab}}=132.0~\mathrm{MeV}")
    xaxis!(L"\theta / \mathrm{deg}",(0,180))
    yaxis!(L"\sigma", (2e-5,2e4),:log)
    plot!(rad2deg.(θlist), σlist, labels="Nuclear",lw=2)
end
savefig(p0,"Fig\\sigma.svg")

# heatmap( rad2deg.(θlist), [1], reshape(log10.(σlist),1,length(θlist)),
            # clims = (-4,2) )

## 散射截面变化
plt_l = 100
θgrid = collect( LinRange(0.001, pi-0.001, 400) )
rgrid = [10^n for n in LinRange(0, 3, 200)]
σmatrix = [σ_r(θ, r, k, Sls, sols, plt_l) for r in rgrid, θ in θgrid]
cmatrix = log10.(σmatrix)

p1 = heatmap(size = (800,600), legend = :bottomright,
    rad2deg.(θgrid), rgrid, cmatrix,
    xaxis = (L"\theta(\degree)", (0,180)),
    yaxis = (L"r(\mathrm{fm})", (rgrid[1], rgrid[end]), :log),
    clim  = (-4, 2),  #default clim
    c = :thermal,
    colorbar_title = (L"\lg(\sigma)")
)
radius = WSp[1].rV * unitlessR .*ones(size(θgrid))
plot!(rad2deg.(θgrid), radius, lw=2, labels="Radius")
savefig(p1,"Fig\\sigma_r.svg")
