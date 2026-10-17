module DeBellyAxisymmetric

using LinearAlgebra
using SparseArrays

export ModelParams, Simulation, run_axisymmetric, write_csv

Base.@kwdef struct ModelParams
    radius::Float64 = 6.0
    theta_max::Float64 = acos(-0.9)
    viscosity::Float64 = 10_000.0
    friction::Float64 = 2.0
    basal_friction::Float64 = 0.0
    sigma0::Float64 = 100.0
    sigma_rho::Float64 = 200.0
    prestress::Float64 = 0.0
    elastic_lambda::Float64 = 17.32
    elastic_mu::Float64 = 17.32
    koff::Float64 = 1.4
    kon::Float64 = 5.0
    mca_total::Float64 = 1.0
    mca_diffusion_bound::Float64 = 5e-4
    mca_diffusion_unbound::Float64 = 0.3
    rac_diffusion::Float64 = 0.05
    rho_diffusion::Float64 = 0.05
    rac_timescale::Float64 = 0.04
    rho_timescale::Float64 = 0.04
    alpha0::Float64 = 1.0
    beta0::Float64 = 1.0
    alpha_mechanical::Float64 = 1.0
    beta_mechanical::Float64 = 1.0
    alpha_opto::Float64 = 2.0
    beta_opto::Float64 = 0.0
    opto_start_step::Int = 50
    opto_stop_step::Int = 150
    opto_width::Float64 = 6.0
    rac_boundary_scale::Float64 = 0.01
    rho_threshold_scale::Float64 = 0.001
    tension_boundary_scale::Float64 = 0.02
    tension_threshold::Float64 = 0.025
    tension_threshold_scale::Float64 = 0.01
    mca_threshold::Float64 = 0.035
    protrusion_speed::Float64 = -0.3
    boundary_relaxation::Float64 = 0.9
    dt::Float64 = 1.0
end

struct MeridianGrid
    nodes::Int
    s::Vector{Float64}
    theta::Vector{Float64}
    tangent_x::Vector{Float64}
    tangent_y::Vector{Float64}
end

struct Simulation
    time::Vector{Float64}
    arclength::Vector{Float64}
    rac::Matrix{Float64}
    rho::Matrix{Float64}
    mcab::Matrix{Float64}
    mcau::Matrix{Float64}
    displacement_x::Matrix{Float64}
    displacement_y::Matrix{Float64}
    tension::Matrix{Float64}
    velocity_x::Matrix{Float64}
    velocity_y::Matrix{Float64}
    mca_relative_error::Float64
end

@inline threshold(x, scale, value) = 0.5 * (tanh(x / scale - value / scale) + 1.0)

function make_grid(elements::Int, p::ModelParams)
    elements >= 8 || throw(ArgumentError("n must be at least 8"))
    theta = collect(range(0.0, p.theta_max; length=elements + 1))
    s = p.radius .* theta
    MeridianGrid(elements + 1, s, theta, sin.(theta), cos.(theta))
end

const GAUSS_X = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
const GAUSS_W = (5 / 9, 8 / 9, 5 / 9)

@inline function geometry_at(s, p)
    theta = s / p.radius
    st, ct = sincos(theta)
    (theta=theta, y=p.radius * st, tx=st, ty=ct, nx=-ct, ny=st)
end

@inline interp(f, a, b, n1, n2) = n1*f[a] + n2*f[b]
@inline deriv(f, a, b, d1, d2) = d1*f[a] + d2*f[b]

function scalar_geometry_matrices(g, p)
    n=g.nodes; rows=Int[]; cols=Int[]; mass=Float64[]; stiff=Float64[]
    for e in 1:n-1
        a,b=e,e+1; h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2; d1,d2=-1/h,1/h
            q=geometry_at(n1*g.s[a]+n2*g.s[b],p); jac=h*wg/2
            for (ii,ni,di) in ((a,n1,d1),(b,n2,d2)), (jj,nj,dj) in ((a,n1,d1),(b,n2,d2))
                push!(rows,ii);push!(cols,jj)
                push!(mass,ni*nj*q.y*jac);push!(stiff,di*dj*q.y*jac)
            end
        end
    end
    sparse(rows,cols,mass,n,n),sparse(rows,cols,stiff,n,n)
end

function weighted_load(f,g,p)
    out=zeros(g.nodes)
    for e in 1:g.nodes-1
        a,b=e,e+1;h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2
            q=geometry_at(n1*g.s[a]+n2*g.s[b],p);jac=h*wg/2
            value=f(a,b,n1,n2,q)
            out[a]+=n1*value*q.y*jac;out[b]+=n2*value*q.y*jac
        end
    end
    out
end

function unweighted_integral(g,field)
    total=0.0
    for e in 1:g.nodes-1
        a,b=e,e+1;h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2
            total+=interp(field,a,b,n1,n2)*h*wg/2
        end
    end
    total
end

function reaction_step(old,antagonist,modulator,basal,coupling,scale,value,timescale,
                       diffusion,M,K,g,p;invert=false)
    source=weighted_load(g,p) do a,b,n1,n2,q
        ant=interp(antagonist,a,b,n1,n2);mod=interp(modulator,a,b,n1,n2)
        gate=threshold(mod,scale,value)
        activation=interp(basal,a,b,n1,n2)+coupling*(invert ? 1-gate : gate)
        activation/(1+ant^2)
    end
    A=(1/(timescale*p.dt)+1)*M+diffusion*K
    A\(M*old/(timescale*p.dt)+source)
end

function active_alpha(g,p,opto)
    fill(p.alpha0,g.nodes) .+
        (opto ? p.alpha_opto .* exp.(-0.5 .* ((g.s .- pi*p.radius)./p.opto_width).^2) : zeros(g.nodes))
end

function solve_velocity(rho,mcab,g,p,front_velocity)
    nd=2g.nodes;rows=Int[];cols=Int[];vals=Float64[];rhs=zeros(nd)
    penalty=10_000.0/(2/40)^2
    for e in 1:g.nodes-1
        a,b=e,e+1;h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2;d1,d2=-1/h,1/h
            q=geometry_at(n1*g.s[a]+n2*g.s[b],p);jac=h*wg/2
            mb=interp(mcab,a,b,n1,n2);rh=interp(rho,a,b,n1,n2)
            sigma=max(0.0,p.sigma0+p.sigma_rho*rh);fr=p.basal_friction+p.friction*mb
            for (ii,ni,di) in ((a,n1,d1),(b,n2,d2)), c in 1:2
                tc=c==1 ? q.tx : q.ty;nc=c==1 ? q.nx : q.ny
                au=di*tc;bu=c==2 ? ni/q.y : 0.0;row=2ii-2+c
                rhs[row]+=-(au+bu)*sigma*q.y*jac
                for (jj,nj,dj) in ((a,n1,d1),(b,n2,d2)), d in 1:2
                    td=d==1 ? q.tx : q.ty;nn=d==1 ? q.nx : q.ny
                    av=dj*td;bv=d==2 ? nj/q.y : 0.0
                    visc=p.viscosity*(2au*av+2bu*bv+au*bv+av*bu)*q.y
                    push!(rows,row);push!(cols,2jj-2+d)
                    push!(vals,(visc+fr*ni*nj*q.y+penalty*ni*nj*nc*nn)*jac)
                end
            end
        end
    end
    A=sparse(rows,cols,vals,nd,nd);bc=[1,2,nd-1,nd]
    # GridapEmbedded orients nΓ into the negative level-set phase, so the
    # authors' τ = J*nΓ is the increasing-arclength tangent. Their boundary
    # value is -vₗτ.
    values=[0.0,0.0,-front_velocity*g.tangent_x[end],-front_velocity*g.tangent_y[end]]
    free=setdiff(1:nd,bc);sol=zeros(nd);sol[bc].=values
    sol[free]=A[free,free]\(rhs[free]-A[free,bc]*values)
    sol
end

function solve_displacement(xold,velocity,mcab,g,p,front_displacement)
    nd=2g.nodes;rows=Int[];cols=Int[];vals=Float64[];rhs=zeros(nd)
    penalty=10.0/(2/40)^2
    for e in 1:g.nodes-1
        a,b=e,e+1;h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2;d1,d2=-1/h,1/h
            q=geometry_at(n1*g.s[a]+n2*g.s[b],p);jac=h*wg/2
            mb=interp(mcab,a,b,n1,n2)
            xo1=interp(xold,2a-1,2b-1,n1,n2);xo2=interp(xold,2a,2b,n1,n2)
            dxo1=deriv(xold,2a-1,2b-1,d1,d2);dxo2=deriv(xold,2a,2b,d1,d2)
            xo_t=xo1*q.tx+xo2*q.ty;ao=dxo1*q.tx+dxo2*q.ty
            # The released helper is cotθ(x)=x[1]/x[2].  Our meridian is
            # x=(-R*cos(theta), R*sin(theta)), hence the leading minus sign.
            cot=-cos(q.theta)/sin(q.theta);cot2=cot^2;cot3=cot^3;csc2=1+cot2
            v1=interp(velocity,2a-1,2b-1,n1,n2);v2=interp(velocity,2a,2b,n1,n2)
            for (ii,ni,di) in ((a,n1,d1),(b,n2,d2)), c in 1:2
                tc=c==1 ? q.tx : q.ty;nc=c==1 ? q.nx : q.ny
                aw=di*tc;wt=ni*tc;row=2ii-2+c
                rhs[row]+=p.friction*mb*(xo1*(c==1)+xo2*(c==2))/p.dt*p.radius*q.y*ni*jac
                rhs[row]+=p.friction*mb*(v1*(c==1)+v2*(c==2))*p.radius*q.y*ni*jac
                for (jj,nj,dj) in ((a,n1,d1),(b,n2,d2)), d in 1:2
                    td=d==1 ? q.tx : q.ty;nn=d==1 ? q.nx : q.ny
                    ax=dj*td;xt=nj*td;dotxw=ni*nj*(c==d);xocomp=d==1 ? xo1 : xo2
                    mass=p.friction*mb/p.dt*p.radius*q.y*dotxw
                    aM=2p.elastic_mu*(dotxw/2+p.radius^2*ax*aw+cot2*dotxw)
                    aM+=(p.elastic_mu/p.radius)*(p.radius*(nj*xocomp+p.radius^2*ao*ax)*aw+cot3*xo_t*dotxw)
                    aL=p.elastic_lambda*(p.radius^2*ax*aw+cot*(aw*xt+ax*wt)+cot2*dotxw)
                    aL+=(0.5p.elastic_lambda/p.radius)*(csc2*nj*xocomp+p.radius^2*ax*ao)*(p.radius*aw+cot*wt)
                    push!(rows,row);push!(cols,2jj-2+d)
                    push!(vals,(mass+(aM+aL)*sin(q.theta)+penalty*ni*nj*nc*nn)*jac)
                end
            end
        end
    end
    A=sparse(rows,cols,vals,nd,nd);bc=[1,2,nd-1,nd]
    values=[0.0,0.0,-front_displacement*g.tangent_x[end],-front_displacement*g.tangent_y[end]]
    free=setdiff(1:nd,bc);sol=zeros(nd);sol[bc].=values
    sol[free]=A[free,free]\(rhs[free]-A[free,bc]*values)
    sol
end

function project_tension(x,g,p,M)
    load=weighted_load(g,p) do a,b,n1,n2,q
        h=g.s[b]-g.s[a];d1,d2=-1/h,1/h
        dx=deriv(x,2a-1,2b-1,d1,d2);dy=deriv(x,2a,2b,d1,d2)
        p.prestress+(p.elastic_lambda+p.elastic_mu)*(dx*q.tx+dy*q.ty)
    end
    M\load
end

function bound_mca_step(old,unbound,velocity,lambda,M,K,g,p)
    n=g.nodes;rows=Int[];cols=Int[];vals=Float64[]
    for e in 1:n-1
        a,b=e,e+1;h=g.s[b]-g.s[a]
        for (z,wg) in zip(GAUSS_X,GAUSS_W)
            n1,n2=(1-z)/2,(1+z)/2;d1,d2=-1/h,1/h
            q=geometry_at(n1*g.s[a]+n2*g.s[b],p);jac=h*wg/2
            vx=interp(velocity,2a-1,2b-1,n1,n2);vy=interp(velocity,2a,2b,n1,n2)
            dvx=deriv(velocity,2a-1,2b-1,d1,d2);dvy=deriv(velocity,2a,2b,d1,d2)
            vt=vx*q.tx+vy*q.ty;av=dvx*q.tx+dvy*q.ty
            for (ii,ni) in ((a,n1),(b,n2)),(jj,nj,dj) in ((a,n1,d1),(b,n2,d2))
                push!(rows,ii);push!(cols,jj);push!(vals,ni*(vt*dj+av*nj)*q.y*jac)
            end
        end
    end
    C=sparse(rows,cols,vals,n,n)
    source=weighted_load(g,p) do a,b,n1,n2,q
        lambda/(pi*p.radius)*(p.kon/(p.kon+p.koff))
    end
    A=M/p.dt+p.mca_diffusion_bound*K+p.koff*M+C
    A\(M*old/p.dt+p.kon*M*unbound+source)
end

function unbound_mca_step(old,bound,lambda,M,K,g,p)
    source=weighted_load(g,p) do a,b,n1,n2,q
        lambda/(pi*p.radius)*(p.koff/(p.kon+p.koff))
    end
    A=M/p.dt+p.mca_diffusion_unbound*K+p.kon*M
    A\(M*old/p.dt+p.koff*M*bound+source)
end

function equilibrate_chemistry(g,p,M,K,mcab,tension)
    rac=zeros(g.nodes);rho=ones(g.nodes)
    alpha=fill(p.alpha0,g.nodes);beta=fill(p.beta0,g.nodes)
    rac=reaction_step(rac,rho,mcab,alpha,p.alpha_mechanical,p.rho_threshold_scale,
                      p.mca_threshold,200p.rac_timescale,p.rac_diffusion,M,K,g,p;invert=true)
    rho=reaction_step(rho,rac,tension,beta,p.beta_mechanical,p.tension_threshold_scale,
                      p.tension_threshold,200p.rho_timescale,p.rho_diffusion,M,K,g,p)
    for _ in 1:100
        rho=reaction_step(rho,rac,tension,beta,p.beta_mechanical,p.tension_threshold_scale,
                          p.tension_threshold,200p.rho_timescale,p.rho_diffusion,M,K,g,p)
        rac=reaction_step(rac,rho,mcab,alpha,p.alpha_mechanical,p.rho_threshold_scale,
                          p.mca_threshold,200p.rac_timescale,p.rac_diffusion,M,K,g,p;invert=true)
    end
    rac,rho
end

function run_axisymmetric(;n::Int=80,duration::Float64=350.0,saveat::Float64=1.0,
                           params::ModelParams=ModelParams(),kwargs...)
    saveat==params.dt||throw(ArgumentError("the published scheme saves every dt"))
    g=make_grid(n,params);M,K=scalar_geometry_matrices(g,params)
    density=params.mca_total/(pi*params.radius)
    mcab=fill(params.kon*density/(params.koff+params.kon),g.nodes)
    mcau=fill(params.koff*density/(params.koff+params.kon),g.nodes)
    displacement=zeros(2g.nodes);tension=zeros(g.nodes)
    rac,rho=equilibrate_chemistry(g,params,M,K,mcab,tension)
    minitial=unweighted_integral(g,mcab)+unweighted_integral(g,mcau)
    steps=round(Int,duration/params.dt)+1
    times=params.dt.*collect(1:steps)
    racout=zeros(steps,g.nodes);rhoout=similar(racout);bout=similar(racout);uout=similar(racout)
    dxout=similar(racout);dyout=similar(racout);tenout=similar(racout)
    vxout=similar(racout);vyout=similar(racout)
    front_v=0.0;front_x=0.0;opto=false;max_mass_error=0.0
    for step in 0:steps-1
        step==params.opto_start_step&&(opto=true)
        step==params.opto_stop_step&&(opto=false)
        velocity=solve_velocity(rho,mcab,g,params,front_v)
        displacement=solve_displacement(displacement,velocity,mcab,g,params,front_x)
        tension=project_tension(displacement,g,params,M)
        j=step+1
        racout[j,:].=rac;rhoout[j,:].=rho;bout[j,:].=mcab;uout[j,:].=mcau
        dxout[j,:].=displacement[1:2:end];dyout[j,:].=displacement[2:2:end]
        tenout[j,:].=tension;vxout[j,:].=velocity[1:2:end];vyout[j,:].=velocity[2:2:end]
        front_tension=tension[end]*threshold(tension[end],0.0001,0.0)
        front_v=params.protrusion_speed*threshold(rac[end],params.rac_boundary_scale,1.3rac[1])/
                (1+front_tension^2/params.tension_boundary_scale)
        front_x=params.boundary_relaxation*front_x+params.dt*front_v
        alpha=active_alpha(g,params,opto);beta=fill(params.beta0,g.nodes)
        rac=reaction_step(rac,rho,mcab,alpha,params.alpha_mechanical,params.rho_threshold_scale,
                          params.mca_threshold,params.rac_timescale,params.rac_diffusion,M,K,g,params;invert=true)
        rho=reaction_step(rho,rac,tension,beta,params.beta_mechanical,params.tension_threshold_scale,
                          params.tension_threshold,params.rho_timescale,params.rho_diffusion,M,K,g,params)
        lambda=minitial-unweighted_integral(g,mcab)-unweighted_integral(g,mcau)
        mcab=bound_mca_step(mcab,mcau,velocity,lambda,M,K,g,params)
        mcau=unbound_mca_step(mcau,mcab,lambda,M,K,g,params)
        mass=unweighted_integral(g,mcab)+unweighted_integral(g,mcau)
        max_mass_error=max(max_mass_error,abs(mass-minitial)/minitial)
    end
    Simulation(times,g.s,racout,rhoout,bout,uout,dxout,dyout,tenout,vxout,vyout,max_mass_error)
end

function write_csv(path::AbstractString,sim::Simulation)
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,"time,arclength,rac,rho,mcab,mcau,displacement_x,displacement_y,tension,velocity_x,velocity_y")
        for j in eachindex(sim.time),i in eachindex(sim.arclength)
            println(io,join((sim.time[j],sim.arclength[i],sim.rac[j,i],sim.rho[j,i],
                sim.mcab[j,i],sim.mcau[j,i],sim.displacement_x[j,i],sim.displacement_y[j,i],
                sim.tension[j,i],sim.velocity_x[j,i],sim.velocity_y[j,i]),','))
        end
    end
    path
end

end
