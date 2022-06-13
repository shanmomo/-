## 绘制散射过程截面图
using Base.Threads  # 多线程支持
# using Plots         # 绘图
# using LaTeXStrings  # 绘图LaTeX支持
# gr()            # 绘图Backend
# theme(:default) # 绘图风格

## 散射截面
θlist = collect( LinRange(0.001, pi-0.001, 1000) )
σlist = [σ(θ, k, ls, Sls, sols)[1] for θ in θlist]
σFar = [σ(θ, k, ls, Sls, sols)[2] for θ in θlist]
σNear = [σ(θ, k, ls, Sls, sols)[3] for θ in θlist]

# begin
#     p0=plot(size = (800,600), minorgrid = true)
#     title!(L"^{16}\mathrm{O}+^{12}\mathrm{C}~~\mathrm{at}~~ E_{\mathrm{lab}}=132.0~\mathrm{MeV}")
#     xaxis!(L"\theta / \mathrm{deg}",(0,180))
#     yaxis!(L"\sigma", (2e-10,2e4),:log)
#     plot!(rad2deg.(θlist), σlist, labels="Nuclear",lw=2)
# end
# savefig(p0,"Fig\\sigma.svg")


## 散射截面变化
plt_l = 100
θgrid = collect( LinRange(0.001, pi-0.001, 400) )
rgrid = [10^n for n in LinRange(0, 3, 200)]
# σmatrix = [σ_r(θ, r, k, Sls, sols, plt_l) for r in rgrid, θ in θgrid]
σmatrix = [σ_r(θ, r, k, Sls, sols, plt_l) for θ in θgrid, r in rgrid]
cmatrix = log10.(σmatrix)

L = argmin(abs.(real.(Sls).-0.5))-1
R12 = sqrt(L*(L+1))/k
radius = R12 .* ones(size(θgrid))

# radius = WSp[1].rV .*ones(size(θgrid))

# p1 = heatmap(size = (800,600), legend = :bottomright,
#     rad2deg.(θgrid), rgrid, cmatrix,
#     xaxis = (L"\theta(\degree)", (0,180)),
#     yaxis = (L"r(\mathrm{fm})", (rgrid[1], rgrid[end]), :log),
#     clim  = (-4, 2),  #default clim
#     c = :thermal,
#     colorbar_title = (L"\lg(\sigma)")
# )
# plot!(rad2deg.(θgrid), radius, lw=2, linecolor = "red",
#         labels=L"\mathrm{Radius}~R_V")
# savefig(p1,"Fig\\sigma_r.svg")

## 拼接
using CairoMakie

begin
	f = Figure(resolution = (800,800), # 大小
	fontsize = 18) #字号
	# 定义2x2网格，行高[2,3]，列宽[9,1]
	rowgrids = [1:2, 3:5] ; colgrids = [1:9, 10]
	ax1 = Axis(f[rowgrids[1], colgrids[1]], ylabel = L"𝜎(𝜃)")
	ax2 = Axis(f[rowgrids[2], colgrids[1]], ylabel = L"r(\mathrm{fm})", xlabel = L"𝛩_{\mathrm{c.m.}}(\deg)")
	ax1.title = "Nuclear Rainbow"
	# 在ax1画线图
	L1 = lines!(ax1, rad2deg.(θlist), σlist, linewidth=2)
	L11 = lines!(ax1, rad2deg.(θlist), σFar, linewidth=1.5, linestyle=:dash)
	L12 = lines!(ax1, rad2deg.(θlist), σNear, linewidth=1.5, linestyle=:dash)
	# Legend(f[rowgrids[1],colgrids[2]], ax1)
	# 在ax2画色图
	hm2 = heatmap!(ax2, rad2deg.(θgrid), rgrid, cmatrix,
		colormap = :thermal,#:Spectral,
		colorrange = (-4, 2) )
	ax3 = Colorbar(f[rowgrids[2], colgrids[2]], hm2, label = L"\lg(𝜎(r,𝜃))")
	L2 = lines!(ax2, rad2deg.(θgrid), radius, linewidth=2, color = :red)
	Legend(f[rowgrids[1],colgrids[2]], [L1,L11,L12,L2], [L"\mathrm{Total}",L"\mathrm{Far}",L"\mathrm{Near}", L"Radius$~R_{1/2}$" ])
	# Legend(f[rowgrids[1],colgrids[2]], [L1,L11,L12,L2], [L"\mathrm{Total}",L"\mathrm{l_<}",L"\mathrm{l_>}", L"Radius$~R_{1/2}$" ])

	# 关联x轴
	linkxaxes!(ax1,ax2)
	xlims!(ax1, low=0, high=180)
	ylims!(ax1, low=2e-6, high=2e4)
	# ylims!(ax1, low=2e-10, high=2e4)

	# 坐标轴格式
	ax1.yscale=log10
	ax1.xaxisposition = :top
	ax2.yscale=log10
	ax1.xminorticksvisible = true
	ax1.xminorgridvisible = true
	ax1.xminorticks = IntervalsBetween(10)
	ax2.xminorticksvisible = true
	ax2.xminorgridvisible = true
	ax2.xminorticks = IntervalsBetween(10)
	# 隐藏一些线
	hidespines!(ax1, :b) #; hidespines!(ax2, :t)
	ax1.xticklabelsvisible = false
	#设置间隙
	rowgap!(f.layout, 0)
	colgap!(f.layout, 10)
	# ax1.yticklabelfont = cmr
	f
end
# save("Fig\\toy_quantum_new.pdf", f, pt_per_unit = 1)

# using Dierckx
# δl = angle.(Sls)./2
# llist = LinRange(0,100,1000)
# spl = Spline1D(ls, δl)
# Θl = 2 .* derivative(spl, llist; nu=1)
#
# using Plots
# plot(llist,Θl)
# plot(ls,δl)

# using CSV
# using DataFrames
#
# Sl = DataFrame( l = ls, Sl = Sls)
# CSV.write("C:\\Users\\Lenovo\\Desktop\\Sl.csv",Sl)
