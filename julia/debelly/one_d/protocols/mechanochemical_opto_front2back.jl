# Standard code to run a simulation using the mechanochemical model presented in article "INSERT DOI"
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# One can work with pure local inhibition by changing the following mechanical parameters:
# vCTE=0   α=0   β=0   σₐ₀=0
# Written by Andreu F Gallen working in Turlier lab and in collaboration with Orion Weiner's lab

using Gridap
using Plots
using Plots.PlotMeasures
include("Plots_RhoRacA.jl")

function plotting(ylab,po,pPNG,i)
  plot(po)
  xlabel!("ξ[μm]")
  ylabel!(ylab)
  savefig(pPNG*"rho"*i*".png")
end

function conservation(sρ0,sρ,Minitial)
  return 0.1*(Minitial-(sρ0+sρ))
end
 
function threshold(x,x₀,xth)
  return  (0.5 * (tanh.(x/x₀ .- xth/x₀).+1)) 
end

#function to run a single simulation with a few given parameters
function run(χ,λ⁻²,η,T,Δt,part) 

  # Time discretisation parameters
  t₀  = 0.0
  t   = t₀  
  nΔt = trunc(Int,T/Δt)

  domain = (0,L)
  partition = (part)
  model = CartesianDiscreteModel(domain,partition) 
  Δx = L/partition

 # Lets write down some parameters in a txt just in case 
  pVTU="./VTU/"*simulation
  write_outputs && mkpath(pVTU)
  pPNG="./PNG/"*simulation
  write_outputs && mkpath(pPNG)
  write_outputs && mkpath(pPNG*"ezrin_time/") 
  write_outputs && mkpath(pPNG*"ezrin_unbound_time/") 
  write_outputs && mkpath(pPNG*"Rac_time/") 
  write_outputs && mkpath(pPNG*"Rho_time/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  write_outputs && cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)


  order = 1
  #Starting Boundary conditions
  x₀(x) = 0
  xₗ(x) = 0

  v₀(x) = 0
  vₗ(x) = 0

  # Creating FE space
  reffe = ReferenceFE(lagrangian,Float64,order)
  W0 = TestFESpace(model,reffe;dirichlet_tags=["tag_1"])
  WD0 = TestFESpace(model,reffe;dirichlet_tags=["tag_1","tag_2"])
  X = TrialFESpace(WD0,[x₀,xₗ])
  V = TrialFESpace(WD0,[v₀,vₗ])

  degree = 2
  #triangulate FE space
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  #Neumann and defining trial spaces and test space for ρ_b 
  Q0 = TestFESpace(model,reffe)
  RHO =  TrialFESpace(Q0)
  RHO0 =  TrialFESpace(Q0)
  Rac =  TrialFESpace(Q0)
  Rho =  TrialFESpace(Q0)

  #Building the vectors used for introducing opto influence as an increase in α and β
  α₀o(x)= αopto*exp(-0.5*(x[1])^2/((wrac)^2)) 
  β₀o(x)= βopto*exp(-0.5*(x[1]-L)^2/((wrac)^2))#2.25*β₀*exp(-0.5*(x[1]-L)^2/((L/4)^2))
  α₀opto = zeros(partition)
  β₀opto = zeros(partition)
  for j in 1:1:(partition)
    α₀opto[j] = α₀o(j*Δx)
    β₀opto[j] = β₀o(j*Δx)
  end
  
  #defining Starting conditions for some variables and dummy variables to be able to build the equations
  uh_ρ = interpolate_everywhere(kon*M0/partition/(koff+kon),RHO)
  uh_ρ_old = uh_ρ
  uh_ρ0 = interpolate_everywhere(koff*M0/partition/(koff+kon),RHO0)
  uh_ρ0_old = uh_ρ0
  uh_x = zero(X)
  uh_x_old = uh_x
  uh_v = zero(V)
  f1(x) = 0.01*x[1]
  f2(x) = α₀+0.5*exp(-x[1]^2/wrac^2)
  uh_rac = interpolate_everywhere(0.0 ,Rac) #0.95
  uh_rho = interpolate_everywhere(0.62 ,Rho) #0.62
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  Minitial=sum(get_free_dof_values(uh_ρ0))+sum(get_free_dof_values(uh_ρ))

  print("Start1 sum(ρ_b) "*string(sum(get_free_dof_values(uh_ρ)))*", sum(ρ_u) "*string(sum(get_free_dof_values(uh_ρ0)))*", and sum(ρ_b+ρ_u)"*string(sum(get_free_dof_values(uh_ρ0))+sum(get_free_dof_values(uh_ρ)))*"\n")
  
  #variables to store temporal information that we would like to plot later
  tensiont = zeros(trunc(Int,T/Δt)+1,partition-2)
  vt = zeros(trunc(Int,T/Δt)+1,partition-1)
  xt = zeros(trunc(Int,T/Δt)+1,partition-1)
  αt = zeros(trunc(Int,T/Δt)+1,partition)
  βt = zeros(trunc(Int,T/Δt)+1,partition)
  Mt = zeros(trunc(Int,T/Δt))
  ract = zeros(trunc(Int,T/Δt)+1,partition+1)
  rhot = zeros(trunc(Int,T/Δt)+1,partition+1)
  #We start wriiting down the values of some parameters at time=0
  λ = 0.0
  ρt = zeros(trunc(Int,T/Δt)+1,partition+1)
  ρt[1,:] = get_free_dof_values(uh_ρ)[:]
  ρ0t = zeros(trunc(Int,T/Δt)+1,partition+1)
  ρ0t[1,:] = get_free_dof_values(uh_ρ0)[:]  
  
  #we write down the weak form of the membrane equation
  # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a 
  a(k,x,w) = ∫(k*(∇(x)⋅∇(w)))dΩ
  m(ρ_b,Δt,x,w) = ∫( (χ*ρ_b)*(x*w)/Δt )dΩ
  aₓ(ρ_b,x,w) = m(ρ_b,Δt,x,w) + a(k,x,w)
  bₓ(ρ_b,v,w) = ∫(χ*(ρ_b*v)*w)dΩ
  
  # Now for MAC bound
  # mass term for the temporal evolution ρ_b
  mρ(Δt,ρ_b,w) = ∫( (ρ_b*w)/Δt )dΩ
  aρ(ρ_b,v,w) = mρ(Δt,ρ_b,w) + ∫(0.001*(∇(ρ_b)⋅∇(w)))dΩ  + ∫(w*(v*(∇(ρ_b)⋅VectorValue(1.0))+ρ_b*(∇(v)⋅VectorValue(1.0))))dΩ + ∫(koff*(ρ_b*w))dΩ +∫(λᵇ*((ρ_b*ρ_b*ρ_b)*w))dΩ
  bρ(w,ρ_u) = ∫(kon*(ρ_u*w)+ λ/L*(kon/(kon+koff)))dΩ + mρ(Δt,uh_ρ_old,w)  

  #Now for MCA unbound
  # mass term for the temporal evolution ρ_b0
  aρ0(ρ_u,x,x_old,w) = mρ(Δt,ρ_u,w) + ∫(kon*(ρ_u*w))dΩ +∫(D*(∇(ρ_u)⋅∇(w)))dΩ  + ∫(w*(((x-x_old)/Δt)*(∇(ρ_u)⋅VectorValue(1.0))+ρ_u*(((∇(x)-∇(x_old))/Δt)⋅VectorValue(1.0))))dΩ  
  bρ0(w,ρ_b) = ∫(koff*(ρ_b*w)+ λ/L*(koff/(kon+koff)))dΩ + mρ(Δt,uh_ρ0_old,w) 

  #Now for v
  aᵥ(ρ_b,v,w) = ∫(η*(∇(v)⋅∇(w)))dΩ + ∫( (χ*ρ_b)*(v*w) )dΩ
  bᵥ(w,uh_rho) = m(uh_ρ,Δt,uh_x,w) - m(uh_ρ,Δt,uh_x_old,w) +  ∫( w*σₐ₀*∇(uh_rho)⋅VectorValue(1.0) )dΩ #no feedback is ∫( w*∇σₐ )dΩ
  Aᵥ(v,w) = aᵥ(uh_ρ,v,w)
  Bᵥ(w) = bᵥ(w,uh_rho)
  #We can now use ρ_b and x to solve v
  op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)

  #SOLVING X AT t=0
  aₓ_0(x,w) = aₓ(0,x,w)
  b_0(w) =  bₓ(0,0,w)
  op_x = AffineFEOperator(aₓ_0,b_0,X,WD0)
  uh_x = solve(op_x)
  uh_x_old =  uh_x

  #SOLVE ρ_b AT t=0 vien initial velocity zero
  Aρ(ρ_b,w) = aρ(ρ_b,uh_v,w) 
  Bρ(w) = bρ(w,uh_ρ0)
  op_ρ= AffineFEOperator(Aρ,Bρ,RHO,Q0)
  uh_ρ=solve(op_ρ)
  uh_ρ_old=uh_ρ

  #SOLVE ρ_u AT t=0
  Aρ0(ρ_u,w) = aρ0(ρ_u,uh_x,uh_x_old,w)
  Bρ0(w) = bρ0(w,uh_ρ_old)
  op_ρ0 = AffineFEOperator(Aρ0,Bρ0,RHO0,Q0)
  uh_ρ0 = solve(op_ρ0)
  uh_ρ0_old = uh_ρ0

  #SOLVING V AT t=0 using x, xold and rho
  uh_v=zero(V)
  
  #define the production rates as vectors so that we can have heterogeneous α and β due to mechanics
  α₀v = zeros(partition)
  α₀v[:] .= α₀
  β₀v = zeros(partition)
  β₀v[:] .= β₀

  #we define the tension for a spring
  lx=get_free_dof_values(uh_x)
  x1=circshift(lx,1)
  tension=k*(lx-x1)/Δx 
  splice!(tension,  1) 

  vt[1,:] = get_free_dof_values(uh_v)[:]
  xt[1,:] = get_free_dof_values(uh_x)[:]    
  tensiont[1,:] = tension[:]
  #ten will be used in the computation of the model, we instill the Boundary Conditions and remove negative values
  ten = zeros(partition, 1)
  ten[1:end-2] = tension[:]
  ten[end-1] = tension[end]
  ten[end] = tension[end]
  
  #storing information of how alpha an beta behave spatially over time
  for j in 1:1:(partition)
    αt[1,j]= α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_ρ)[j])/rho0-ρth/rho0))) 
  end
  βt[1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 


  #DEFINING the equations for Rac and Rho
  a_rac(rac,w,v) = (1/dᵃ)*mρ(Δt,rac,w) +∫(Drac*(∇(rac)⋅∇(w)))dΩ  + ∫(w*rac)dΩ + ∫(mTorc*w*rac*threshold(ten,ten0,10))dΩ # + ∫(1/dᵃ*w*((v)*(∇(rac)⋅VectorValue(1.0))+rac*(∇(v)⋅VectorValue(1.0))))dΩ 
  b_rac(w,rho,ρ_b) =  (1/dᵃ)*mρ(Δt,uh_rac_old,w) + ∫(w*(α₀v/(1+rho*rho) + α*(0.5*(1-tanh∘(ρ_b/rho0-ρth/rho0)))/(1+rho*rho)))dΩ  #*(tanh∘(0.9 /(ρ_b) - 0.9/0.07)+1)
 
  a_rho(rho,w,v) = (1/dᵇ)*mρ(Δt,rho,w) +∫(Drho*(∇(rho)⋅∇(w)))dΩ  + ∫(w*rho)dΩ  +∫(λʳᴬ*((rho*rho*rho)*w))dΩ  #+ ∫(1/dᵇ*w*((v)*(∇(rho)⋅VectorValue(1.0))+rho*(∇(v)⋅VectorValue(1.0))))dΩ
  b_rho(w,ten,rac) =  (1/dᵇ)*mρ(Δt,uh_rho_old,w) + ∫( w*(β₀v/(1+rac*rac) + β*threshold(ten,sig0,tenth)/(1+rac*rac)))dΩ  

 
  #SOLVE Rac AT t=0
  Arac(rac,w) = a_rac(rac,w,uh_v)
  Brac(w) = b_rac(w,uh_rho,uh_ρ)
  op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_rho(rho,w,uh_v)
  Brho(w) = b_rho(w,ten,uh_rac)
  op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
  uh_rho = solve(op_rho)
  uh_rho_old = uh_rho


  i = 0
  t=0
  dummyx0=0
 #give steady state as initial conditions for rac and rho 
  for ti in 1:250
    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho
  end
  #threshold value to start protrusion, 1.3 times the initial condition
  racth = 1.3*get_free_dof_values(uh_rac)[1] 

  
  ract[1,:] = get_free_dof_values(uh_rac)[:]
  rhot[1,:] = get_free_dof_values(uh_rho)[:]
  
  write_outputs && writevtk(Ω,pVTU*"VTU$i",cellfields=["x"=>uh_x,"v"=>uh_v,"rho"=>uh_ρ,"rho0"=>uh_ρ0]) 
  for ti in t₀:Δt:(T-Δt)
    #HERE WE DEFINE WHETHER THE CODE IS FRONT TO BACK OR BACK TO FRONT, DEPENDING IN WHERE WE ACTIVATE OPTO
    if t==topto
      α₀v[:] .= α₀ .+ α₀opto[:]
      β₀v[:] .= β₀ .+ β₀opto[:]
    end
    if t==3*topto #at time=3*topto a while the input dies down
      α₀v[:] .= α₀
      β₀v[:] .= β₀ 
    end
    i1 = sum(get_free_dof_values(uh_ρ))
    i2 = sum(get_free_dof_values(uh_ρ0))
    λ = conservation(i2,i1,Minitial)
    i3 = i1+i2
    i4 = trunc(t)
    i = i + 1
    t = t + Δt  

    #boundary conditions for v and x
    rac1=get_free_dof_values(uh_rac)[1]
    vCTErac = vCTE*threshold(rac1,rac0,racth)#velocity polimerization
    v₀ = vCTErac/(1+ten[1]*ten[1]/ten0)
    #we introduce a slight relaxation for the membrane, decreases 2% x at the Boundary condition only
    dummyx0=dummyx0*0.98+Δt*vCTErac/(1+ten[1]*ten[1]/ten0)
    x₀ = dummyx0
    X = TrialFESpace(WD0,[x₀,xₗ])

    @info "Time step $i/$nΔt, time $i4, sum(ρ_b+ρ_u) $i3"
    Mt[i]=i3

    ρt[i+1,:] = get_free_dof_values(uh_ρ)[:]
    ρ0t[i+1,:] = get_free_dof_values(uh_ρ0)[:]
    ract[i+1,:] = get_free_dof_values(uh_rac)[:]
    rhot[i+1,:] = get_free_dof_values(uh_rho)[:]

    # Updating v to solve ρ and x
    A(x,w) = m(uh_ρ,Δt,x,w) + a(k,x,w)
    B(w) = m(uh_ρ,Δt,uh_x,w) + bₓ(uh_ρ,uh_v,w) 
    op_x = AffineFEOperator(A,B,X,WD0)
    uh_x = solve(op_x)

    po=get_free_dof_values(uh_ρ)
    pa=get_free_dof_values(uh_rac)
    paa=get_free_dof_values(uh_rho)
    poo = get_free_dof_values(uh_ρ0)

    op_ρ= AffineFEOperator(Aρ,Bρ,RHO,Q0)
    uh_ρ=solve(op_ρ)
    uh_ρ_old=uh_ρ
  
    op_ρ0= AffineFEOperator(Aρ0,Bρ0,RHO0,Q0)
    uh_ρ0=solve(op_ρ0)
    uh_ρ0_old=uh_ρ0
 
    V = TrialFESpace(WD0,[v₀,vₗ])
    op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)
    uh_v=solve(op_v)

    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho

    #updating x_old per time derivarive
    uh_x_old =  uh_x
    write_outputs && writevtk(Ω,pVTU*"VTU$i",cellfields=["x"=>uh_x,"v"=>uh_v,"rho"=>uh_ρ,"rho0"=>uh_ρ0])

    write_outputs && plotting("ρ_b",po,pPNG*"ezrin_time/","$i")
    write_outputs && plotting("rac",pa,pPNG*"Rac_time/","$i")
    write_outputs && plotting("rac",paa,pPNG*"Rho_time/","$i")
    write_outputs && plotting("ρ_u",poo,pPNG*"ezrin_unbound_time/","$i")

    lx=get_free_dof_values(uh_x)
    x1=circshift(lx,1)
    tension=k*(lx-x1)/Δx
    splice!(tension,  1)
    ten[1:end-2] = tension[:]
    ten[end-1] = tension[end]
    ten[end] = tension[end] 
    vt[i+1,:] = get_free_dof_values(uh_v)[:]
    xt[i+1,:] = get_free_dof_values(uh_x)[:]
    tensiont[i+1,:] = tension[:]
    βt[i+1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 
    for j in 1:1:(partition)
      αt[i+1,j] = α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_ρ)[j])/rho0-ρth/rho0))) 
    end
  end 
  write_outputs && plots_run(nΔt,Δx,vt,xt,tensiont,ρt,ρ0t,ract,rhot,Mt,λ⁻²,pPNG,α, β, dᵃ, dᵇ,αt,βt,αopto,βopto,topto)
  return tensiont,ρt,vt,ract,rhot,αt,βt
end


###############CONSTANT MECHANICAL PARAMETERS ##############
partition=parse(Int, get(ENV, "DEBELLY_PARTITION", "100")) # size of FE system
const write_outputs = get(ENV, "DEBELLY_WRITE_OUTPUTS", "1") == "1"
χ = 100               # [pN s/um^3] friction coefficient
η=10000               # [pN s/ um] 2D viscosity of the cortex
const L = 20.           # [um] System size length
k=100.0            # [pN /um] elastic constant for the membrane
σₐ₀ =  100        # [pN /um^2] maximum active "pressure"

const koff = 1.4  # ezrin koff [1/s] toff=0.7s Fritzsche et al
const kon  = 5.0  # ezrin kon [1/s] ton =0.2s Fritzsche et al
D = 0.3   # diffusion of the membrane [um^2/s] 0.003  Fritzsche et al
M0 = 10.0 #total amount of ezrin
const wrac=L/5 #if we want to start polarized we use this for the Gaussian distribution
#in case one wants nonlinear deactication if the system would be unstable
λᵇ = 0   # for bound ezrin 
λʳᴬ= 0   #for Rho

#adimensional numbers that define the system
τₐ = η/σₐ₀ #for σₐ₀=100
τₑ = η/k #almost always?
Pe⁻¹ = D*η/(σₐ₀*L^2) #D*η/(σₐ₀*L^2) 0.00075 for σₐ₀=100
λ = sqrt( η*L /(χ*L*L*M0) ) #\bar λ in the suppl material. Adimensional lenthscale

rac0=0.01# ADIMENSINAlizes rac switch for protrusion, changes the slope of the threshold function
vCTE=-0.15 #Scale for the polimerizatio velocity
ten0=10.0 #denominator for protrusion dependence on tension
tenth=5.0 #threshold for tension to activate rho
sig0 = 5.00  #adimensionalizes the tension term for rho
ρth = 0.068 #MCA/ezrin value below which protrusion starts
rho0 = 0.01 #adimensionalizes the MCA concentration term for rac

#Rac Coefficients
dᵃ = 0.04  #deactivation
α₀ = 1.00  #basal activation
α =0.2*α₀ #activation through mechanics, zero for only local inhibition
Drac = 0.5 #diffusion rate
mTorc = 0.0#1
αopto = 1 # opto input for Rac

#Rho Coefficients
dᵇ = 0.04 #deactivation
β₀ = 1.00 #basal activation
β =0.2*β₀ #activation through mechanics, zero for only local inhibition
Drho = 0.5 #diffusion rate 
βopto = 0 # opto input for Rho

#time variables
topto=parse(Float64, get(ENV, "DEBELLY_OPTO_TIME", "50"))
Δt=parse(Float64, get(ENV, "DEBELLY_DT", "1"))
T=parse(Float64, get(ENV, "DEBELLY_T", "350"))
#if we wanna do multiple simulations in a row with different variables
setsv=1
setsx=2 

#output folder name
len=trunc(λ,digits=2)
const simulation= "Opto_activation/T=$T ten0=$ten0 db=$dᵇ da=$dᵃ sig0=$sig0 rth=$ρth len=$len/bopto=$βopto aopto=$αopto rho0=$rho0 mTorc=$mTorc D=$D a0=$α₀ b0=$β₀ te=$τₑ ta=$τₐ M0=$M0 deltat=$Δt vCTE=$vCTE tenth=$tenth lb=$λᵇ/"
#store VTUs in one folder
pVTU="./VTU/"*simulation
write_outputs && mkpath(pVTU)
#Store PNGs in another
pPNG="./PNG/"*simulation
write_outputs && mkpath(pPNG)

#IF we run multiple simulations with changing parameters, we use this arrays to store the results of each simulation
xχ = zeros(setsx)
xη = zeros(setsv)
tension = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
ρ_b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
a = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
tension2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
ρ2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
αt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)
βt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)

#FOR for the possible multiple simulations
for i in 2:1:setsx
  for j in 1:1:setsv
    # χ = 2*10^((i))    # [pN s/um^3] friction coefficient
    # χ = trunc(χ,digits=2)
    λ⁻² = 1/(λ*λ)
    xχ[i]=χ 
    #η= 0+0*j+10^(1*(j+1))              # [pN s/ um] 2D viscosity of the cortex
    lη = η
    xη[j]=lη
    tension[i,j,:,:],ρ_b[i,j,:,:],v[i,j,:,:],a[i,j,:,:],b[i,j,:,:],αt[i,j,:,:],βt[i,j,:,:]= run(χ,λ⁻²,lη,T,Δt,partition)
  end
end

println("DEBELLY_SUMMARY protocol=mechanochemical rac_front=$(a[2,1,end,1]) rac_back=$(a[2,1,end,end]) rho_front=$(b[2,1,end,1]) rho_back=$(b[2,1,end,end])")
