using Plots
using LaTeXStrings
gr()            # 绘图Backend
theme(:default) # 绘图风格

# ## 势函数
# rlist = collect(LinRange(0,15. ,500))
# Vclist = [Vc(r) for r in rlist]
# WSreallist = [real(normal_WS1(r)) for r in rlist]
# cenVlist = [realcenV(r) for r in rlist]
#
# p0=plot(size = (800,500), minorgrid = true,
#         titlefontsize=18,legendfontsize=12,
#         guidefontsize=14,tickfontsize=10,
#         legend=:bottomright)
#     title!("Real potential")
#     plot!(rlist,Vclist,lw=1.5,labels="Coulomb potential")
#     plot!(rlist,WSreallist,lw=1.5,labels="Nuclear potential")
#     plot!(rlist,cenVlist,lw=3,labels="Total potential")
#     xaxis!(L"r(\mathrm{fm})")
#     yaxis!(L"V(\mathrm{MeV})")
#
# savefig(p0,"RealPotential.svg")
#
# ## 运动轨迹 (直角坐标系)
# getux(sol) = [cos(u[2])*u[1] for u in sol.u]
# getuy(sol) = [sin(u[2])*u[1] for u in sol.u]
# getcolor(sol) = [cos.(u[5]) for u in sol.u]
# getalpha(sol) = [0.2*exp.(-u[6]) for u in sol.u]
#
# xdatas = [getux(sol) for sol in sols]
# ydatas = [getuy(sol) for sol in sols]
# colordatas = [getcolor(sol) for sol in sols]
# alphadatas = [getalpha(sol) for sol in sols]
#
# p1=plot(size = (900,600), minorgrid = true,
#         legend = false, axis = nothing, titlefontsize=20,
#         xlim = (-15,45), ylim = (-20,20))
# for ii in 1:n
#      plot!(xdatas[ii],ydatas[ii],
#      lw = 150*Δb[ii],
#      line_z = colordatas[ii],
#      seriesalpha = alphadatas[ii])
# end
# for ii in 1:n
#      plot!(xdatas[ii],-ydatas[ii],
#      lw = 150*Δb[ii],
#      line_z = colordatas[ii],
#      seriesalpha = alphadatas[ii])
# end
# title!("Trajectory with absorption")
# savefig(p1,"TrajWizAb.png")
#
# # 无吸收轨迹
# getcolor2(sol) = [cos.(u[7]) for u in sol.u]
# colordatas2 = [getcolor2(sol) for sol in sols]
#
# p2=plot(size = (900,600), minorgrid = true,
#         legend = false, axis = nothing, titlefontsize=20,
#         xlim = (-15,45), ylim = (-20,20))
# for ii in 1:n
#      plot!(xdatas[ii],ydatas[ii],
#      lw = 150*Δb[ii],
#      line_z = colordatas2[ii],
#      seriesalpha = 0.2)
# end
# for ii in 1:n
#      plot!(xdatas[ii],-ydatas[ii],
#      lw = 150*Δb[ii],
#      line_z = colordatas2[ii],
#      seriesalpha = 0.2)
# end
# title!("Trajectory without absorption")
# savefig(p2,"TrajWizoutAb.png")

## θ-b
θlist = [rad2deg.(sol[end][2]) for sol in sols]
p3=plot(size = (800,500), minorgrid = true,
        titlefontsize=18,legendfontsize=12,
        guidefontsize=14,tickfontsize=10,
        framestyle = :origin,legend=false,
        guide_position = :top)
title!(L"\Theta -b")
plot!(b,θlist,lw=2)
xaxis!(L"b(\mathrm{fm})")
yaxis!(L"\Theta(\degree)")
savefig(p3,"Thetab.svg")

θ = [rad2deg.(π .- abs.(mod.(sol[end][2],2π) .-π)) for sol in sols]
p4=plot(size = (800,500), minorgrid = true,
        titlefontsize=18,legendfontsize=12,
        guidefontsize=14,tickfontsize=10,
        framestyle = :origin,legend=false)
title!(L"\theta -b")
plot!(b,θ,lw=2)
xaxis!(L"b(\mathrm{fm})")
yaxis!(L"\theta(\degree)")
savefig(p4,"θb.svg")

# ## 运动轨迹 (极坐标系)
# rdatas = [sol[1,:] for sol in sols]
# θdatas = [sol[2,:] for sol in sols]
#
# p3=plot(size = (800,800),legend=false)
# for ii in 1:n
#      plot!(θdatas[ii],rdatas[ii], proj = :polar,
#      lw = 30*Δb[ii],
#      line_z = colordatas[ii],
#      seriesalpha = alphadatas[ii])
# end
# for ii in 1:n
#      plot!(2π .-θdatas[ii],rdatas[ii], proj = :polar,
#      lw = 30*Δb[ii],
#      line_z = colordatas[ii],
#      seriesalpha = alphadatas[ii])
# end
