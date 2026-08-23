
using Plots
using Plots.PlotMeasures
#plots for singlet simulation
function plots_run_singlet(nΔt,vt,#ezrinbt,ezrinut,
  ract,rhot,pPNG,
   partition,L,Δt,T,xplot,σₐt,χt)
 # xplot = range(0, L, length=partition)
 
  tplot = range(0, T, length=nΔt)
  ract_mirror = zeros(nΔt,partition*2)
  MCAt_mirror = zeros(nΔt,partition*2)
  rhot_mirror = zeros(nΔt,partition*2)
  x2plot = range(-L, L, length=partition*2) 

  for i in 1:1:nΔt
    for j in 1:1:(partition)
      ract_mirror[i,j]=ract[i,partition-j+1]
      ract_mirror[i,2*partition+1-j]=ract[i,partition-j+1]
      rhot_mirror[i,j]=rhot[i,partition-j+1]
      rhot_mirror[i,2*partition+1-j]=rhot[i,partition-j+1]
     # MCAt_mirror[i,j]=ezrinbt[i,partition-j+1]
     # MCAt_mirror[i,2*partition+1-j]=ezrinbt[i,partition-j+1]
    end 
  end
  heatmap(tplot,x2plot,transpose(MCAt_mirror),  size=(250, 220),margin=15px, plot_title="MCA protein",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_MCA_mi.pdf")
  
  # mi = minimum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # ma = maximum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # heatmap(tplot,x2plot,transpose(ract_mirror./ract_mirror[30,50]), clim=(1.0*mi ,1.3 ),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(ract_mirror),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac_mi.pdf")

  heatmap(tplot,xplot,transpose(ract),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac.pdf")

  # mi =  minimum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])
  # ma = maximum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  #heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(rhot_mirror),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho_mirror.pdf")

  heatmap(tplot,xplot,transpose(rhot),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho.pdf")

  n=18
  plot_lines =  trunc(Int,nΔt/(n-8))

  vplot=vt[1:plot_lines:end,:]
  p5 = plot(xplot,transpose(vplot), legend=false,color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  # ρtplot=ezrinbt[5:plot_lines:end,:]
  # p2 = plot(xplot,transpose(ρtplot), legend=false, color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  # #ylims!(0.045, 0.12)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA protein")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot.pdf")

  ractplot=ract[:2:plot_lines:end,:]
  p2 = plot(xplot,transpose(ractplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rac R")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[:2:plot_lines:end,:]
  p1 =plot(xplot,transpose(rhotplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rho ρ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf") 
  
   σₐtplot=σₐt[:2:plot_lines:end,:]
  p3 = plot(xplot,transpose(σₐtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("σₐt")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"sigmat.pdf")

  χtplot=χt[:2:plot_lines:end,:]
  p4 =plot(xplot,transpose(χtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("χ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"chit.pdf") 
  
  plot(p1, p2,p3,p4,p5, layout = 5, size=(800, 400))
  savefig(pPNG*"all.pdf") 
   
  ractplot=ract[:,1]./ract[1,1]
  ractplot2=ract[:,partition]./ract[1,1]
  plot(tplot,[ractplot ractplot2], label=["Back" "Front"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 2.75)
  xlabel!("time (sec)")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rac_2.pdf")

  rhotplot=rhot[:,1]./rhot[1,1]
  rhotplot2=rhot[:,partition]./rhot[1,1]
  plot(tplot,[rhotplot rhotplot2], label=["Back" "Front"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rho_2.pdf")
  

  ρtplot=ezrinbt[:,1]./ezrinbt[1,1]
  ρtplot2=ezrinbt[:,partition]./ezrinbt[1,1]
  plot(tplot,[ρtplot ρtplot2], label=["Back" "Front"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("MCA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"MCA_2.pdf")
  

  # ρtplot=ezrinut[:1:plot_lines:end,:]
  # plot(xplot,transpose(ρtplot), legend=false,palette = :hawaii25,size=(300, 220), margin = 5px,linealpha=0.3,linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("ρ_u")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rho0t.pdf")
  
  # plot(xplot, [ezrinut[end,:], ezrinbt[end,:]], label=["ρ_u(ξ,t=$T)" "ρ_b(ξ,t=$T)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rho.pdf")
  
  # plot(xplot, [ezrinut[1,:], ezrinbt[1,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhoINITIAL.pdf")

  # plot(xplot, [ezrinut[2,:], ezrinbt[2,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot1.pdf")

  
  vplot=vt[end,:]
  p1 = plot(xplot,vplot, legend=false,size=(300, 220), margin = 5px, alpha = 0.9)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  # ρtplot=ezrinbt[end,:]
  # p2 = plot(xplot,ρtplot, legend=false,size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("ρ_b")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot.pdf")

  ractplot=ract[end,:]
  p5 = plot(xplot,ractplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rac")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[end,:]
  p6 =plot(xplot,rhotplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rhoA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")
  

  plot(p1, p5,p6, layout = 3, plot_titlefontsize=10,size=(800, 400))
  savefig(pPNG*"all_t=$T.pdf")
end

function plots_run_singlet(nΔt,vt,vnt,xt,ezrinbt,#ezrinut,
  ract,rhot,pPNG,
   partition,L,Δt,T,xplot,σₐt,χt,tent)
 # xplot = range(0, L, length=partition)
 
  tplot = range(0, T, length=nΔt)
  ract_mirror = zeros(nΔt,partition*2)
  MCAt_mirror = zeros(nΔt,partition*2)
  rhot_mirror = zeros(nΔt,partition*2)
  x2plot = range(-L, L, length=partition*2) 

  for i in 1:1:nΔt
    for j in 1:1:(partition)
      ract_mirror[i,j]=ract[i,partition-j+1]
      ract_mirror[i,2*partition+1-j]=ract[i,partition-j+1]
      rhot_mirror[i,j]=rhot[i,partition-j+1]
      rhot_mirror[i,2*partition+1-j]=rhot[i,partition-j+1]
     MCAt_mirror[i,j]=ezrinbt[i,partition-j+1]
     MCAt_mirror[i,2*partition+1-j]=ezrinbt[i,partition-j+1]
    end 
  end
  heatmap(tplot,x2plot,transpose(MCAt_mirror),  size=(250, 220),margin=15px, plot_title="MCA protein",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_MCA_mi.pdf")
  
  # mi = minimum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # ma = maximum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # heatmap(tplot,x2plot,transpose(ract_mirror./ract_mirror[30,50]), clim=(1.0*mi ,1.3 ),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(ract_mirror),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac_mi.pdf")

  heatmap(tplot,xplot,transpose(ract),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac.pdf")

  # mi =  minimum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])
  # ma = maximum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  #heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(rhot_mirror),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho_mirror.pdf")

  heatmap(tplot,xplot,transpose(rhot),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho.pdf")

  n=18
  plot_lines =  trunc(Int,nΔt/(n-8))

  vplot=vt[1:plot_lines:end,:]
  p5 = plot(xplot,transpose(vplot), legend=false,color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (radians)")
  ylabel!("v tan")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  vplot=vnt[1:plot_lines:end,:]
  plot(xplot,transpose(vplot), legend=false,color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (radians)")
  ylabel!("v normal")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vnt.pdf")


  xtplot=xt[1:plot_lines:end,:]
  p6 = plot(xplot,transpose(xtplot), legend=false,color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (radians)")
  ylabel!("x tan")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"xt.pdf")

  ρtplot=ezrinbt[5:plot_lines:end,:] 
  p4 = plot(xplot,transpose(ρtplot), legend=false, color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  #ylims!(0.045, 0.12)
  ylims!(0.5*ezrinbt[1,1], 1.7*ezrinbt[1,1])
  xlabel!("Cell perimeter (μm)")
  ylabel!("MCA protein")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhot.pdf")

  for i in 1:1:nΔt
    for j in 1:1:(partition) 
     MCAt_mirror[i,j]=tent[i,partition-j+1]
     MCAt_mirror[i,2*partition+1-j]=tent[i,partition-j+1]
    end 
  end
  heatmap(tplot,x2plot,transpose(MCAt_mirror),  size=(250, 220),margin=15px, plot_title="tension",plot_titlefontsize=10, framestyle = :box)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  savefig(pPNG*"ten_mi.pdf")
  
  tenplot=tent[5:plot_lines:end,:]
  p3 = plot(xplot,transpose(tenplot), legend=false, color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  #ylims!(0.045, 0.12)
  xlabel!("Cell perimeter (radians)")
  ylabel!("Tension")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ten.pdf")

  ractplot=ract[:2:plot_lines:end,:]
  p2 = plot(xplot,transpose(ractplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # # ylims!(0, Ntot/L)
  #    ylims!(0.5, 1.7)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rac R")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

 
  rhotplot=rhot[:2:plot_lines:end,:]
  p1 =plot(xplot,transpose(rhotplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  #  ylims!(0.5, 1.7)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rho ρ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf") 
  
   σₐtplot=σₐt[:2:plot_lines:end,:]
  plot(xplot,transpose(σₐtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("σₐt")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"sigmat.pdf")

  χtplot=χt[:2:plot_lines:end,:]
  plot(xplot,transpose(χtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("χ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"chit.pdf") 
  
  plot(p1, p2,p3,p4,p5, p6, layout = 6, size=(800, 400))
  savefig(pPNG*"all.pdf") 
   
  ractplot=ract[:,1]./ract[1,1]
  ractplot2=ract[:,partition]./ract[1,1]
  plot(tplot,[ractplot ractplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 2.75)
  xlabel!("time (sec)")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rac_2.pdf")

  rhotplot=rhot[:,1]./rhot[1,1]
  rhotplot2=rhot[:,partition]./rhot[1,1]
  plot(tplot,[rhotplot rhotplot2], label=["Back" "Front"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rho_2.pdf")
 


  ρtplot=ezrinbt[:,1]./ezrinbt[1,1]
  ρtplot2=ezrinbt[:,partition]./ezrinbt[1,1]
  plot(tplot,[ρtplot ρtplot2], label=["Back" "Front"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("MCA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"MCA_2.pdf")
  
  vplot=vt[end,:]
  p1 = plot(xplot,vplot, legend=false,size=(300, 220), margin = 5px, alpha = 0.9)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  ρtplot=ezrinbt[end,:]
  p2 = plot(xplot,ρtplot, legend=false,size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("ρ_b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhot.pdf")

  ractplot=ract[end,:]
  p5 = plot(xplot,ractplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rac")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[end,:]
  p6 =plot(xplot,rhotplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rhoA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")
  

  plot(p1, p5,p6, layout = 3, plot_titlefontsize=10,size=(800, 400))
  savefig(pPNG*"all_t=$T.pdf")
end


function plots_run_singlet_turnover(nΔt,vt,#ezrinbt,ezrinut,
  ract,rhot,pPNG,
   partition,L,Δt,T,xplot,ett,σₐt,χt)
 # xplot = range(0, L, length=partition)
 
  tplot = range(0, T, length=nΔt)
  ract_mirror = zeros(nΔt,partition*2)
  MCAt_mirror = zeros(nΔt,partition*2)
  rhot_mirror = zeros(nΔt,partition*2)
  x2plot = range(-L, L, length=partition*2) 

  for i in 1:1:nΔt
    for j in 1:1:(partition)
      ract_mirror[i,j]=ract[i,partition-j+1]
      ract_mirror[i,2*partition+1-j]=ract[i,partition-j+1]
      rhot_mirror[i,j]=rhot[i,partition-j+1]
      rhot_mirror[i,2*partition+1-j]=rhot[i,partition-j+1]
     # MCAt_mirror[i,j]=ezrinbt[i,partition-j+1]
     # MCAt_mirror[i,2*partition+1-j]=ezrinbt[i,partition-j+1]
    end 
  end
  heatmap(tplot,x2plot,transpose(MCAt_mirror),  size=(250, 220),margin=15px, plot_title="MCA protein",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_MCA_mi.pdf")
  
  # mi = minimum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # ma = maximum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  # heatmap(tplot,x2plot,transpose(ract_mirror./ract_mirror[30,50]), clim=(1.0*mi ,1.3 ),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(ract_mirror),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac_mi.pdf")

  heatmap(tplot,xplot,transpose(ract),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac.pdf")

  # mi =  minimum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])
  # ma = maximum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  #heatmap(tplot,x2plot,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  heatmap(tplot,x2plot,transpose(rhot_mirror),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho_mirror.pdf")

  heatmap(tplot,xplot,transpose(rhot),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("time (sec)")
  ylabel!("Cell perimeter (μm)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho.pdf")

  n=18
  plot_lines =  trunc(Int,nΔt/(n-8))

  eplot=ett[1:plot_lines:end,:]
  p3 = plot(xplot,transpose(eplot), legend=false,color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("e")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  # ρtplot=ezrinbt[5:plot_lines:end,:]
  # p2 = plot(xplot,transpose(ρtplot), legend=false, color = :matter, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  # #ylims!(0.045, 0.12)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA protein")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot.pdf")

  ractplot=ract[:2:plot_lines:end,:]
  p2 = plot(xplot,transpose(ractplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[:2:plot_lines:end,:]
  p1 =plot(xplot,transpose(rhotplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")  

  σₐtplot=σₐt[:2:plot_lines:end,:]
  p4 = plot(xplot,transpose(σₐtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("σₐt")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"sigmat.pdf")

  χtplot=χt[:2:plot_lines:end,:]
  p5 =plot(xplot,transpose(χtplot), legend=false,color =cgrad(:matter, rev=false), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("χ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"chit.pdf") 
  
  plot(p1, p2,p3,p4,p5, layout = 5, size=(800, 400))
  savefig(pPNG*"all.pdf") 
   
  ractplot=ract[:,1]./ract[1,1]
  ractplot2=ract[:,partition]./ract[1,1]
  plot(tplot,[ractplot ractplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 2.75)
  xlabel!("time (sec)")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rac_2.pdf")

  rhotplot=rhot[:,1]./rhot[1,1]
  rhotplot2=rhot[:,partition]./rhot[1,1]
  plot(tplot,[rhotplot rhotplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rho_2.pdf")

  ρtplot=ezrinbt[:,1]./ezrinbt[1,1]
  ρtplot2=ezrinbt[:,partition]./ezrinbt[1,1]
  plot(tplot,[ρtplot ρtplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("time (sec)")
  ylabel!("MCA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"MCA_2.pdf")
  

  # ρtplot=ezrinut[:1:plot_lines:end,:]
  # plot(xplot,transpose(ρtplot), legend=false,palette = :hawaii25,size=(300, 220), margin = 5px,linealpha=0.3,linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("ρ_u")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rho0t.pdf")
  
  # plot(xplot, [ezrinut[end,:], ezrinbt[end,:]], label=["ρ_u(ξ,t=$T)" "ρ_b(ξ,t=$T)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rho.pdf")
  
  # plot(xplot, [ezrinut[1,:], ezrinbt[1,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhoINITIAL.pdf")

  # plot(xplot, [ezrinut[2,:], ezrinbt[2,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("MCA (t=$T)")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot1.pdf")

  
  vplot=vt[end,:]
  p1 = plot(xplot,vplot, legend=false,size=(300, 220), margin = 5px, alpha = 0.9)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  # ρtplot=ezrinbt[end,:]
  # p2 = plot(xplot,ρtplot, legend=false,size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1)
  # # ylims!(0, Ntot/L)
  # xlabel!("Cell perimeter (μm)")
  # ylabel!("ρ_b")
  # # plot!(legend=:topright, legendcolumns=3)
  # savefig(pPNG*"rhot.pdf")

  ractplot=ract[end,:]
  p2 = plot(xplot,ractplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rac")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[end,:]
  p3 =plot(xplot,rhotplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("Cell perimeter (μm)")
  ylabel!("rhoA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")
  

  plot(p1, p2,p3, layout = 3, plot_titlefontsize=10,size=(800, 400))
  savefig(pPNG*"all_t=$T.pdf")
end
