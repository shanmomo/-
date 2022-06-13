using Plots
using LaTeXStrings
gr()            # 绘图Backend
theme(:default) # 绘图风格

##

# Θ = [sol[end][2] for sol in sols]
# # Θ = [mod2pi(sol[end][2]) for sol in sols]
# # for i = 1:length(Θ)
# #     if Θ[i] > π
# #         Θ[i] += -2π
# #     end
# # end
#
# p3=plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         framestyle = :origin,legend=false,
#         guide_position = :top)
#         title!(L"\Theta -l")
#         plot!(l,rad2deg.(Θ),lw=2)
#         xaxis!(L"l")
#         yaxis!(L"\Theta(\degree)")

## θ-b
# θlist = [rad2deg.(sol[end][2]) for sol in sols]
# p3=plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         framestyle = :origin,legend=false,
#         guide_position = :top)
# title!(L"\Theta -b")
# plot!(b,θlist,lw=2)
# xaxis!(L"b(\mathrm{fm})")
# yaxis!(L"\Theta(\degree)")
# savefig(p3,"Thetab.svg")

## θ-b

θ = [rad2deg.(π .- abs.(mod.(theta(sol),2π) .-π)) for sol in sols]

# begin
#     p1 = plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         framestyle = :origin,legend=false)
#         title!(L"\theta -b")
#         plot!(b,θ,lw=2)
#         xaxis!(L"b(\mathrm{fm})")
#         yaxis!(L"\theta(\degree)")
# end
# # savefig(p1,"θb.svg")

## ϕ-b

# ϕ = [ (phase(sol,k0)-π)/π for sol in sols]
ϕ = [ phase(sol,k0) for sol in sols]

# begin
#     p2 = plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         framestyle = :origin,legend=false)
#         title!(L"\phi -b")
#         plot!(b, ϕ, lw=2 )
#         plot!(b, zeros(size(b)), ls = :dash, color = :red )
#         xaxis!(L"b(\mathrm{fm})")
#         yaxis!(L"\phi(\pi)")
# end
# savefig(p2,"ϕb.svg")

## A-b

A = [exp.(-sol[end][6]) for sol in sols]

# begin
#     p3 = plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         framestyle = :origin,legend=false)
#         title!(L"A -b")
#         plot!(b, A, lw=2 )
#         xaxis!(L"b(\mathrm{fm})")
#         yaxis!(L"A")
# end

## 合并
using CairoMakie

begin
	f = Figure(resolution = (900,900), # 大小
	fontsize = 18) #字号
	# 定义3x1网格，行高[3,3,3]，列宽[9]
	# rowgrids = [1:2, 3:4, 5:6] ; colgrids = [1:9]
	rowgrids = [1:2, 3:6] ; colgrids = [1:6]
	ax1 = Axis(f[rowgrids[1], colgrids[1]], ylabel = L"𝜃(\degree)", title = "Semi-classical Method")
	# ax2 = Axis(f[rowgrids[2], colgrids[1]], ylabel = L"\cos(𝜙)" )
	# ax3 = Axis(f[rowgrids[3], colgrids[1]], ylabel = L"𝐴", xlabel = L"𝑙")
	ax2 = Axis(f[rowgrids[2], colgrids[1]], ylabel = L"𝐴\cos(𝜙)", xlabel = L"𝑙" )

	L1_0 = lines!(ax1, l, 78*ones(size(l)), linewidth=2, linestyle=:dash, color=:red)
	L1_1 = lines!(ax1, [14.6,14.6], [0,78], linewidth=2, linestyle=:dash, color=:red)
	L1_2 = lines!(ax1, [26.65,26.65], [0,78], linewidth=2, linestyle=:dash, color=:red)
	L1 = lines!(ax1, l, θ, linewidth=2)
	axislegend(ax1, [L1], ["散射角"])
	L2_0 = lines!(ax2, l, zeros(size(l)), linestyle=:dash, color=:black)
	L2_1 = lines!(ax2, [14.6,14.6], [0,0.5], linewidth=2, linestyle=:dash, color=:red)
	L2_2 = lines!(ax2, [26.65,26.65], [0,0.5], linewidth=2, linestyle=:dash, color=:red)
	L2 = lines!(ax2, l, A.*cos.(ϕ), linewidth=2)
	axislegend(ax2, [L2], ["干涉信息"])
	# L3 = lines!(ax3, l, A, lw=2)
	# axislegend(ax3, [L3], ["幅值"])

	# 关联x轴
	linkxaxes!(ax1,ax2,ax3)
	xlims!(ax1, low=10, high=30)
	ylims!(ax1, low=0, high=180)
	ylims!(ax2, low=-0.1, high=0.1)
	# ylims!(ax3, low=0, high=1)

	# 坐标轴格式
	ax1.xaxisposition = :top
	ax1.xminorticksvisible = true
	ax1.xminorgridvisible = true
	ax1.xminorticks = IntervalsBetween(10)
	ax2.xminorticksvisible = true
	ax2.xminorgridvisible = true
	ax2.xminorticks = IntervalsBetween(10)
	ax2.yminorticksvisible = true
	ax2.yminorgridvisible = true
	ax2.yminorticks = IntervalsBetween(10)
	# ax3.xminorticksvisible = true
	# ax3.xminorgridvisible = true
	# ax3.xminorticks = IntervalsBetween(10)

	# 隐藏一些线
	hidespines!(ax1, :b)
	# hidexdecorations!(ax2, minorgrid = false)
	ax1.xticklabelsvisible = false
	#设置间隙
	rowgap!(f.layout, 0)
	# lengend
	f
end
# save("Fig\\toy model\\toy_classical.pdf", f, pt_per_unit = 1)
