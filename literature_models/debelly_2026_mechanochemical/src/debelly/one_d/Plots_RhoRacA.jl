# Kymographs may requite a bit of playing around their colorbar range clim=(min, max) depending on the opto signals
# Otherwise they may look either too big or too small for the introduced signal

using Gridap
using Plots
using Plots.PlotMeasures
using DelimitedFiles

#plots for the simulations that do heatmap
function plots_heat(time,α₀,β₀,nΔt,Δx,vt,xt,tensiont,ρt,ρ0t,ract,rhot,Mt,λ⁻²,pPNG,α, β, dᵃ, dᵇ,αt,βt)
  xplotp1 = range(0, L, length=partition+1)
  xplotm1 = range(0, L, length=partition-1)
  xplotm2 = range(0, L, length=partition-2)  

  n=18
  plot_lines =  trunc(Int,time/(n-8))
  if plot_lines==0
    plot_lines=1
  end
  vplot=vt[1:plot_lines:time,:]
  p1 = plot(xplotm1,transpose(vplot), legend=false,color = :haline, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)

  ρtplot=ρt[5:plot_lines:time,:]
  p2 = plot(xplotp1,transpose(ρtplot), legend=false, color = :haline, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  ylims!(0.045, 0.12)
  xlabel!("ξ[μm]")
  ylabel!("MCA protein")
  # plot!(legend=:topright, legendcolumns=3)
 
  xxplot=xt[1:plot_lines:time,:] 
  p4 = plot(xplotm1,transpose(xxplot), legend=false, color = :haline, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.8, framestyle = :box)
  # ylims!(0, Ntot/L)
  plot!(cbar=true)
  xlabel!("ξ[μm]")
  ylabel!("x[μm]")
  # plot!(legend=:topright, legendcolumns=3)
    
  tensiontplot=tensiont[6:plot_lines:time,:]
  p3 = plot(xplotm2,transpose(tensiontplot), legend=false,color =:haline, line_z = (1:n+1)',size=(300, 220), margin = 5px, framestyle = :box)
  #ylims!(1, 11)
  xlabel!("ξ[μm]")
  ylabel!("σ[pN/μm]")
  # plot!(legend=:topright, legendcolumns=3)

  ractplot=ract[:4:plot_lines:time,:]
  p5 = plot(xplotp1,transpose(ractplot), legend=false,color =:haline, line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)

  rhotplot=rhot[:4:plot_lines:time,:]
  p6 =plot(xplotp1,transpose(rhotplot), legend=false,color =:haline, line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)

  λ=trunc(sqrt(1/λ⁻²),digits=2)
  plot(p1, p2, p5,p4,p3,p6, layout = 6, plot_title="t2s=$time s λ=$λ τₑ=$τₑ τₐ=$τₐ vp=$vCTE   σ₀=$ten0 ρ₀=$rho0 α=$α β=$β dₐ=$dᵃ dᵦ=$dᵇ",plot_titlefontsize=10,size=(1200, 600))
  savefig(pPNG*"alpha=$α₀ beta=$β₀.png") 
end

#standard array of plots for standard simulation
function plots_run(nΔt,Δx,vt,xt,tensiont,ρt,ρ0t,ract,rhot,Mt,λ⁻²,pPNG,pRac, pRho, koffrac, koffrho,αt,βt,αopto,βopto,topto)
  xplotp1 = range(0, L, length=partition+1)
  xplotm1 = range(0, L, length=partition-1)
  xplotm2 = range(0, L, length=partition-2)
  xplot = range(0, L, length=partition)
  tplot = range(0, T, length=(trunc(Int,T/Δt)+1))
  ract_mirror = zeros((trunc(Int,T/Δt)+1),partition*2+2)
  MCAt_mirror = zeros((trunc(Int,T/Δt)+1),partition*2+2)
  rhot_mirror = zeros((trunc(Int,T/Δt)+1),partition*2+2)
  tent_mirror = zeros((trunc(Int,T/Δt)+1),partition*2-4)
  x2plotp1 = range(-L, L, length=partition*2+2)
  x2plotm2 = range(-L, L, length=partition*2-4)

  for i in 1:1:(trunc(Int,T/Δt)+1)
    for j in 1:1:(partition+1)
      ract_mirror[i,j]=ract[i,partition+2-j]
      ract_mirror[i,2*partition+3-j]=ract[i,partition+2-j]
      rhot_mirror[i,j]=rhot[i,partition+2-j]
      rhot_mirror[i,2*partition+3-j]=rhot[i,partition+2-j]
      MCAt_mirror[i,j]=ρt[i,partition+2-j]
      MCAt_mirror[i,2*partition+3-j]=ρt[i,partition+2-j]
    end 
    for j in 1:1:(partition-2) 
      tent_mirror[i,j]=tensiont[i,partition-1-j]
      tent_mirror[i,2*partition-3-j]=tensiont[i,partition-1-j]
      
    end
  end
  heatmap(tplot,x2plotp1,transpose(MCAt_mirror), clim=(0.055 ,0.095),  size=(250, 220),margin=15px, plot_title="MCA protein",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_MCA_mi.pdf")

  heatmap(tplot,x2plotm2,transpose(tent_mirror),  size=(250, 220),margin=15px, plot_title="σ",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_ten_mi.pdf")
  
  mi = minimum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  ma = maximum(ract_mirror[:,:]./ract_mirror[Int32(topto*4/5),50])
  heatmap(tplot,x2plotp1,transpose(ract_mirror./ract_mirror[30,50]), clim=(1.0*mi ,(1/(αopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac_mi.pdf")

  heatmap(tplot,xplotp1,transpose(ract),  size=(250, 220),margin=15px, plot_title="Rac",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rac.pdf")

  mi =  minimum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])
  ma = maximum(rhot_mirror[:,:]./rhot_mirror[Int32(topto*4/5),50])
  heatmap(tplot,x2plotp1,transpose(rhot_mirror./rhot_mirror[30,50]), clim=(1*mi, (1/(βopto/2+1))*ma ),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho_mirror.pdf")

  heatmap(tplot,xplotp1,transpose(rhot),  size=(250, 220),margin=15px, plot_title="Rho",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_rho.pdf")

  n=18
  plot_lines =  trunc(Int,nΔt/(n-8))

  vplot=vt[1:plot_lines:end,:]
  p1 = plot(xplotm1,transpose(vplot), legend=false,color = :berlin, line_z = (1:n)',size=(300, 220), margin = 5px, alpha = 0.9, framestyle = :box)
  #plot!(cbar=true)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  ρtplot=ρt[5:plot_lines:end,:]
  p2 = plot(xplotp1,transpose(ρtplot), legend=false, color = :haline, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1, framestyle = :box)
  ylims!(0.045, 0.12)
  xlabel!("ξ[μm]")
  ylabel!("MCA protein")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhot.pdf")

#  print((maximum(xt[1:end,:]))) 
  xxplot=xt[1:plot_lines:end,:]
 # print("\n")
  #print(xxplot)
  p4 = plot(xplotm1,transpose(xxplot), legend=false, color = :haline, line_z = (1:n)',size=(300, 220), margin = 5px,linealpha=0.8, framestyle = :box)
  # ylims!(0, Ntot/L)
  plot!(cbar=true)
  xlabel!("ξ[μm]")
  ylabel!("x[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"xt.pdf")
    
  tensiontplot=tensiont[6:plot_lines:end,:]
  p3 = plot(xplotm2,transpose(tensiontplot), legend=false,color =:haline, line_z = (1:n+1)',size=(300, 220), margin = 5px, framestyle = :box)
  #ylims!(1, 11)
  xlabel!("ξ[μm]")
  ylabel!("σ[pN/μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensiont.pdf") 

  ractplot=ract[:2:plot_lines:end,:]
  p5 = plot(xplotp1,transpose(ractplot), legend=false,color =cgrad(:matter, rev=true), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[:2:plot_lines:end,:]
  p6 =plot(xplotp1,transpose(rhotplot), legend=false,color =cgrad(:matter, rev=true), line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")
  α2= trunc(α,digits=2)
  β2= trunc(β,digits=2)
  λ=trunc(sqrt(1/λ⁻²),digits=2)
  plot(p1, p2, p5,p4,p3,p6, layout = 6, plot_title="T=$T s λ=$λ τₑ=$τₑ τₐ=$τₐ vp=$vCTE   σ₀=$ten0 ρ₀=$rho0 α=$α2 β=$β2 dₐ=$dᵃ dᵦ=$dᵇ",plot_titlefontsize=10,size=(800, 400))
  savefig(pPNG*"all.pdf") 
   
  ractplot=ract[:,1]./ract[1,1]
  ractplot2=ract[:,partition]./ract[1,1]
  plot(tplot,[ractplot ractplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 2.75)
  xlabel!("t[s]")
  ylabel!("Rac a")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rac_2.pdf")

  rhotplot=rhot[:,1]./rhot[1,1]
  rhotplot2=rhot[:,partition]./rhot[1,1]
  plot(tplot,[rhotplot rhotplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("t[s]")
  ylabel!("Rho b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rho_2.pdf")
  
  myarray=open(readdlm,"a_betavsalpha_diagram.txt")
  x = myarray[:,1]
  y = myarray[:,2]
  plot((βt[:2:plot_lines:end,:,end]),(αt[:2:plot_lines:end,:,1]), seriestype=:scatter ,markerstrokewidth=0,  markersize = 0.6, legend=false,markercolor =cgrad(:matter, rev=true),size=(300, 220), margin = 1px, framestyle = :box)
  plot!(x,y, legend=false, size=(300, 220), seriestype=:scatter ,markerstrokewidth=0,  markercolor = :red,  markersize = 0.4, margin = 1px)
  ylims!(0, maximum(y))
  xlims!(0, maximum(x))
  ylabel!("α") 
  xlabel!("β")
  savefig(pPNG*"alphabeta.pdf")


  ρtplot=ρt[:,1]./ρt[1,1]
  ρtplot2=ρt[:,partition]./ρt[1,1]
  plot(tplot,[ρtplot ρtplot2], label=["Front" "Back"], grid=false,size=(200, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylims!(0.5, 1.7)
  xlabel!("t[s]")
  ylabel!("MCA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"MCA_2.pdf")
  

  heatmap(tplot,xplot,transpose(αt), clim=(0,1.2*(α₀ + α)), size=(320, 220),margin=15px, plot_title="α",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_alpha.pdf")


  heatmap(tplot,xplot,transpose(βt), clim=(0,1.2*(β₀ + β) ),  size=(320, 220),margin=15px, plot_title="β",plot_titlefontsize=10, framestyle = :box)
  # ylims!(0, Ntot/L)
  xlabel!("t[s]")
  ylabel!("ξ[μm]")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"kymo_beta.pdf")

  bplot=αt[:2:plot_lines:end,:]
  plot(xplot,transpose(bplot),  legend=false,color = :berlin25, line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylabel!("α") 
  xlabel!("ξ")
  savefig(pPNG*"alpha.pdf")

  bplot=βt[:2:plot_lines:end,:]
  plot(xplot,transpose(bplot),  legend=false,color = :berlin25, line_z = (1:n)',size=(300, 220), margin = 5px,linewidth=1, framestyle = :box)
  ylabel!("β") 
  xlabel!("ξ")
  savefig(pPNG*"beta.pdf")
 
 

  
  plot(Mt, size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("t(s)")
  ylabel!("M = ∫(ρ_b+ρ_u)dξ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"Mt.pdf")

 
  
  plot(xplotp1, [ρ0t[end,:], ρt[end,:]], label=["ρ_u(ξ,t=$T)" "ρ_b(ξ,t=$T)"],size=(300, 200),  margin = 5px, linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("MCA (t=$T)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rho.pdf")
  
  plot(xplotp1, [ρ0t[1,:], ρt[1,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("MCA (t=$T)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhoINITIAL.pdf")

  plot(xplotp1, [ρ0t[2,:], ρt[2,:]], label=["ρ_u(ξ,t=1)" "ρ_b(ξ,t=1)"],size=(300, 200),  margin = 5px, linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("MCA (t=$T)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhot1.pdf")
  
  tensiontplot=tensiont[:,:]
  plot((-tensiontplot[:,end]),(tensiontplot[:,1]), legend=false,palette = :berlin25,size=(300, 220), margin = 5px)
  # ylims!(0, Ntot/L)
  xlabel!("-σ₂")
  ylabel!("σ₁")
 # plot!(legend=:outerright, legendcolumns=3)
  savefig(pPNG*"sig1vssig2.pdf")

  
  vplot=vt[end,:]
  p1 = plot(xplotm1,vplot, legend=false,size=(300, 220), margin = 5px, alpha = 0.9)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  ρtplot=ρt[end,:]
  p2 = plot(xplotp1,ρtplot, legend=false,size=(300, 220), margin = 5px,linealpha=0.9,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("ρ_b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhot.pdf")

  ractplot=ract[end,:]
  p5 = plot(xplotp1,ractplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("rac")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"ract.pdf")

  rhotplot=rhot[end,:]
  p6 =plot(xplotp1,rhotplot, legend=false,size=(300, 220), margin = 5px,linewidth=1)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("rhoA")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoAt.pdf")
  

  xxplot=xt[end,:]
  p3 = plot(xplotm1,(xxplot), legend=false,size=(300, 220), margin = 5px,linealpha=0.8)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("x")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"xt.pdf")
    
  tensiontplot=tensiont[end,:]
  p4 = plot(xplotm2,(tensiontplot), legend=false,size=(300, 220), margin = 5px)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Tension")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensiont.pdf") 

  plot(p1, p2, p5,p4,p3,p6, layout = 6, plot_title="T=$T λ⁻²=$λ⁻²  σₐ₀=$σₐ₀  k=$k  η=$η",plot_titlefontsize=10,size=(800, 400))
  savefig(pPNG*"all_t=$T.pdf") 
end

#plots for simulation without rac and rho, only mechanics
function plots_run2(nΔt,vt,xt,tensiont,λ⁻²,pPNG)
  xplotp1 = range(0, L, length=partition+1)
  xplotm1 = range(0, L, length=partition-1)
  xplotm2 = range(0, L, length=partition-2)
  xplot = range(0, L, length=partition)
  plot_lines =  trunc(Int,nΔt/20)
  vplot=vt[1:plot_lines:end,:]
  p1 = plot(xplotm1,transpose(vplot), legend=false,palette = :hawaii25,size=(300, 220), margin = 1px)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("v")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"vt.pdf")

  xxplot=xt[:1:plot_lines:end,:]
  p3 = plot(xplotm1,transpose(xxplot), legend=false,palette = :hawaii25,size=(300, 220), margin = 1px,linealpha=0.8)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("x")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"xt.pdf")

  tensiontplot=tensiont[:1:plot_lines:end,:]
  p4 = plot(xplotm2,transpose(tensiontplot), legend=false,palette = :hawaii25,size=(300, 220), margin = 1px)
  # ylims!(0, Ntot/L)
  xlabel!("ξ[μm]")
  ylabel!("Tension")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensiont.pdf") 

  tensiontplot=tensiont[:,:]
  p2 = plot((-tensiontplot[:,end]),(tensiontplot[:,1]), legend=false,palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("-σ₂")
  ylabel!("σ₁") 
  savefig(pPNG*"sig1vssig2.pdf")


  plot(p1, p2, p3,p4, layout = 4, plot_title="λ⁻²=$λ⁻²   k=$k  η=$η",plot_titlefontsize=12)
  savefig(pPNG*"all.pdf") 
end

#plots for an array pf simulations where friction and viscosity change for each simulation
function plots_end()
  xplotm1 = range(0, L, length=partition-1)
  xplotm2 = range(0, L, length=partition-2)

  #FINAL PLOTS ARTICLE
  whicheta = 1 
  tensiontplot=tension[2,whicheta,end,:] 
  tensiontplot3=ρ_b[2,whicheta,end,:]
  plot((tensiontplot),  label="σ", palette = :hawaii25,size=(300, 200), margin = 2px)
  xlabel!("ξ[μm]")
  ylabel!("Tension σ[pN/μm]")
  axis2 = twinx();
  plot!(axis2,(tensiontplot3), label="ρᵇ",legend=:bottomleft, palette = :berlin25, margin = 2px)
  ylabel!(axis2,"ρᵇ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tenMCA.pdf")
 
  whicheta = 1 
  auxplot1=a[2,whicheta,end,:] 
  auxplot2=b[2,whicheta,end,:]
  plot((auxplot1),  label="Rac",   palette = :hawaii25,size=(300, 200), margin = 2px)
  xlabel!("ξ[μm]")
  ylabel!("Rac")
  axis2 = twinx();
  plot!(axis2,(auxplot2), label="Rho",legend=:bottomleft,  palette = :berlin10, margin = 2px)
  ylabel!(axis2,"Rho")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"RhoRac.pdf")


  tensiontplot=tension[:,:,end,end]
  plot(xχ,tensiontplot, label=leg2,xscale=:log10,palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Tension(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionvisco_theta=180.pdf")

  tensiontplot=tension[:,1,end,end]/(maximum(abs.(tension[:,1,end,end])))
  tensiontplot2=tension2[:,1,end,end]/(maximum(abs.(tension2[:,1,end,end])))
  plot(xχ,[tensiontplot,tensiontplot2], title="T=$T s  M0=$M0 η="*string(trunc(Int,η))*"\nk="*string(trunc(k,digits=0))*" σₐ₀="*string(trunc(σₐ₀,digits=1))*" p vel=$vCTE"
  ,label=["Heterogeneous ρ_b" "Homogeneous ρ_b"],
  xscale=:log10, size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Tension σ₂")
  plot!(legend=:topright)
  savefig(pPNG*"tension_theta=180.pdf")

  tensiontplot=tension[:,:,end,trunc(Int,partition*1/2)] 
  plot(xχ,tensiontplot, label=leg2,xscale=:log10,palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Tension(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionvisco_theta=90.pdf")

  tensiontplot=tension[:,1,end,trunc(Int,partition*1/2)]/(maximum(abs.(tension[:,1,end,trunc(Int,partition*1/2)])))
  tensiontplot2=tension2[:,1,end,trunc(Int,partition*1/2)]/(maximum(abs.(tension2[:,1,end,trunc(Int,partition*1/2)])))
  plot(xχ,[tensiontplot,tensiontplot2], title="θ=90° η=$η",label=["Heterogeneous ρ_b" "Homogeneous ρ_b"],xscale=:log10,
  size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"tension_theta=90.pdf")


  tensiontplot=tension[:,:,end,2] 
  plot(xχ,tensiontplot, label=leg2,xscale=:log10,palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Tension(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionvisco_theta=1.pdf")

  tensiontplot=tension[:,1,end,2]/(maximum(abs.(tension[:,1,end,2])))
  tensiontplot2=tension2[:,1,end,2]/(maximum(abs.(tension2[:,1,end,2])))
  plot(xχ,[tensiontplot,tensiontplot2]
  , title="T=$T s  M0=$M0 η="*string(trunc(Int,η))*"\nk="*string(trunc(k,digits=0))*" σₐ₀="*string(trunc(σₐ₀,digits=1))*" p vel=$vCTE"
  ,label=["Heterogeneous ρ_b" "Homogeneous ρ_b"],xscale=:log10,
  size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Tension σ₁")
  plot!(legend=:topleft)
  savefig(pPNG*"tension_theta=1.pdf")

  tensiontplot=tension[:,1,end,2]/(maximum(abs.(tension[:,1,end,2])))
  tensiontplot2=tension[:,1,end,end]/(maximum(abs.(tension[:,1,end,end])))
  tensiontplot3=ρ_b[:,1,end,end]/(maximum(abs.(ρ_b[:,1,end,end])))
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3], 
  title="η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀))
  ,xscale=:log10, size=(300, 220), margin = 1px, label=["σ₁" "σ₂" "ρₘ"],alpha=0.9)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"Nfinal_T.pdf")

  tensiontplot=tension[:,1,end,2]
  tensiontplot2=tension[:,1,end,end]
  tensiontplot3=ρ_b[:,1,end,end]
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3], title="η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀))
  ,xscale=:log10, size=(300, 220), margin = 1px, label=["σ₁" "σ₂" "ρₘ"],alpha=0.9)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"final_T.pdf")

  tensiontplot=tension[:,1,(trunc(Int64,T/(4*Δt))),2]/(maximum(abs.(tension[:,1,(trunc(Int64,T/(4*Δt))),2])))
  tensiontplot2=tension[:,1,(trunc(Int64,T/(4*Δt))),end]/(maximum(abs.(tension[:,1,(trunc(Int64,T/(4*Δt))),end])))
  tensiontplot3=ρ_b[:,1,(trunc(Int64,T/(4*Δt))),end]/(maximum(abs.(ρ_b[:,1,(trunc(Int64,T/(4*Δt))),end])))
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3], title="η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀))
  ,xscale=:log10, size=(300, 220), margin = 1px, label=["σ₁" "σ₂" "ρₘ"],alpha=0.9)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"Nfinal_Tq.pdf")

  tensiontplot=tension[:,1,(trunc(Int64,T/(4*Δt))),2]
  tensiontplot2=tension[:,1,(trunc(Int64,T/(4*Δt))),end]
  tensiontplot3=ρ_b[:,1,(trunc(Int64,T/(4*Δt))),end]
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3], title="η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀))
  ,xscale=:log10, size=(300, 220), margin = 1px, label=["σ₁" "σ₂" "ρₘ"],alpha=0.9)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"final_Tq.pdf")

  tensiontplot=tension[:,1,(trunc(Int64,T/(4*Δt))),2]/(maximum(abs.(tension[:,1,:,2])))
  tensiontplot2= tension[:,1,(trunc(Int64,2*T/(4*Δt))),2]/(maximum(abs.(tension[:,1,:,2])))
  tensiontplot3= tension[:,1,(trunc(Int64,3*T/(4*Δt))),2]/(maximum(abs.(tension[:,1,:,2])))
  tensiontplot4= tension[:,1,end,2]/(maximum(abs.(tension[:,1,:,2])))
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3,tensiontplot4], title="θ=1° η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀)),
  label=["T="*string((trunc(Int64,T/(4)))) "T="*string((trunc(Int64,T/2))) "T="*string((trunc(Int64,3*T/(4)))) "T=$T" ]
  ,xscale=:log10, size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"tension_theta=1_T.pdf")

  tensiontplot=tension[:,1,(trunc(Int64,T/(4*Δt))),end]/(maximum(abs.(tension[:,1,(trunc(Int64,T/(4*Δt))),end])))
  tensiontplot2= tension[:,1,(trunc(Int64,2*T/(4*Δt))),end]/(maximum(abs.(tension[:,1,(trunc(Int64,2*T/(4*Δt))),end])))
  tensiontplot3= tension[:,1,(trunc(Int64,3*T/(4*Δt))),end]/(maximum(abs.(tension[:,1,(trunc(Int64,3*T/(4*Δt))),end])))
  tensiontplot4= tension[:,1,end,end]/(maximum(abs.(tension[:,1,end,end])))
  plot(xχ,[tensiontplot,tensiontplot2,tensiontplot3,tensiontplot4], title="θ=180°  η="*string(trunc(Int,η))*" k="*string(trunc(Int,k))*" σₐ₀="*string(trunc(Int,σₐ₀)),
  label=["T="*string((trunc(Int64,T/(4)))) "T="*string((trunc(Int64,T/2))) "T="*string((trunc(Int64,3*T/(4)))) "T=$T" ]
  ,xscale=:log10, size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("Normalized Tension")
  plot!(legend=:topright)
  savefig(pPNG*"tension_theta=180_T.pdf")

  whicheta = 1 
  tensiontplot=tension[1:1:end,whicheta,end,:] 
  tensiontplot3=ρ_b[1:1:end,whicheta,end,:]
  plot(transpose(tensiontplot),  label=lag2, legend=false, legendcolumns=2, palette = :hawaii25,size=(500, 300), margin = 2px)
  xlabel!("ξ[μm]")
  ylabel!("Tension σ(ξ,λ⁻²)")
  axis2 = twinx();
  plot!(axis2,transpose(tensiontplot3), label=false, legend=false,palette = :hawaii25, margin = 2px)
  ylabel!(axis2,"ρ_b(ξ,λ⁻²)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionfrict_x.pdf")

  whicheta = 1
  tensiontplot=tension[:,whicheta,end,:] 
  tensiontplot3=ρ_b[:,whicheta,end,1:1:end-3]
  plot(transpose(tensiontplot+tensiontplot3),  label=lag2, legend=false, legendcolumns=2, palette = :hawaii25,size=(500, 300), margin = 2px)
  xlabel!("ξ[μm]")
  ylabel!("Tension + ρ_b")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionfrict_SUM_x.pdf")

  whicheta = 1
  tensiontplot=tension[:,whicheta,end,:] 
  tensiontplot3=2.5*ρ_b[:,whicheta,end,1:1:end-3]
  plot(transpose(tensiontplot+tensiontplot3),  label=lag2, legend=false, legendcolumns=2, palette = :hawaii25,size=(500, 300), margin = 2px)
  xlabel!("ξ[μm]")
  ylabel!("Tension + ϵ⋅ρ")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionfrict_SUM2_x.pdf")


  tensiontplot=tension[2,:,end,:]
  tensiontplot2=tension2[2,:,end,:]
  plot(xplotm2,transpose(tensiontplot), label=leg2, legend=:topright, palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("ξ[μm]")
  ylabel!("Tension(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"tensionvisco_x.pdf")

  ρplot=ρ_b[:,:,end,end]
  plot(xχ,ρplot, label=leg2,xscale=:log10, palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("ρ_b(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"rhovisco.pdf")

  vplot=v[:,:,end,end]
  vplot2=v2[:,:,end,end]
  plot(xχ,vplot,label=leg2,xscale=:log10, palette = :berlin25,size=(300, 220), margin = 1px)
  xlabel!("Friction λ⁻² [pN⋅s/um^3]")
  ylabel!("v(η)")
  # plot!(legend=:topright, legendcolumns=3)
  savefig(pPNG*"Vvisco.pdf")
end