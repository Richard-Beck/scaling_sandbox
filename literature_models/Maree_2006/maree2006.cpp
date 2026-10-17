// Close reconstruction of Maree et al., Bull. Math. Biol. 68 (2006) 1169-1211.
// See README.md and spec.md for documented and inferred implementation choices.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

constexpr double PI = 3.14159265358979323846;
constexpr int ORI = 6;
enum Field { C, R, RH, CI, RI, RHI, ARP, NF };
enum class AreaMode { PublishedLiteral, CPMLatticeArea, PhysicalAreaEquivalent };

const char* area_mode_name(AreaMode mode){
    if(mode==AreaMode::PublishedLiteral)return "published_literal";
    if(mode==AreaMode::CPMLatticeArea)return "cpm_lattice_area";
    return "physical_area_equivalent";
}

struct Config {
    int nx = 500, ny = 500, initial_sites = 30000, steps = 1200;
    uint64_t seed = 1;
    double dx = 100.0, mcs_dt = 0.1, pde_dt = 2e-4, reaction_dt = 0.002;
    double active_diffusion_dt = 0.02, arp_diffusion_dt = 0.01, inactive_diffusion_dt = 0.01;
    double cg_tolerance = 1e-9;
    int cg_max_iterations = 200;
    double pulse_seconds = 10.0, pulse_fraction = 0.15;
    double xi = 0.06, snapshot_seconds = 10.0;
    double target_area_um2 = -1, lambda_cpm = 1e-9, contact_j = 7.5e-4;
    double temperature = 0.008, yield_energy = 0.046;
    double initial_p_total = 0;
    bool pulse = true, frozen = false, snapshots = true, implicit = true;
    bool binary_snapshots = false;
    AreaMode area_mode = AreaMode::CPMLatticeArea;
    std::string output = "run";
};

struct Model {
    struct CopyDiagnostics {
        long candidates=0,connectivity_rejects=0,ext_trials=0,ret_trials=0,ext_accepted=0,ret_accepted=0;
        double ext_expected=0,ret_expected=0,ext_min_dh=std::numeric_limits<double>::infinity(),ret_min_dh=std::numeric_limits<double>::infinity();
        double ext_area=0,ret_area=0,ext_contact=0,ret_contact=0,ext_push=0,ret_push=0,ext_rho=0,ret_rho=0;
    } copydiag;
    Config cfg;
    int n = 0, cell_sites = 0;
    double time = 0.0, target_area = 3e8, target_sites = 30000;
    std::vector<uint8_t> cell;
    std::array<std::vector<int>, ORI> nb;
    std::vector<std::array<int, 36>> third;
    std::array<std::vector<double>, NF> q, qnew;
    std::array<std::vector<double>, ORI> F, B, P, Fnew, Bnew, Pnew;
    std::vector<int> active, order;
    std::vector<double> cg_r, cg_p, cg_Ap;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniform{0.0, 1.0};
    std::ofstream metrics, snapshots;
    double initial_cx = 0, initial_cy = 0;
    double pulse_lo = 0, pulse_hi = 1;
    std::array<double, NF> max_scalar_geom_error{};
    std::array<double, ORI> max_f_geom_error{}, max_bp_geom_error{};
    std::uint64_t cg_solves = 0, cg_iterations = 0;
    int cg_max_used = 0;
    double cg_worst_relative_residual = 0;

    explicit Model(Config c) : cfg(std::move(c)), n(cfg.nx * cfg.ny), rng(cfg.seed) {
        if (cfg.nx < 16 || cfg.ny < 16 || cfg.initial_sites < 7 || cfg.initial_sites >= n / 2)
            throw std::runtime_error("invalid lattice or initial cell size");
        if (!cfg.implicit && (cfg.pde_dt <= 0 || cfg.pde_dt > cfg.dx * cfg.dx / (4.0e7) + 1e-15))
            throw std::runtime_error("pde_dt violates the positivity CFL bound dx^2/(4 D_max)");
        if (cfg.reaction_dt <= 0 || cfg.active_diffusion_dt <= 0 || cfg.arp_diffusion_dt <= 0 || cfg.inactive_diffusion_dt <= 0)
            throw std::runtime_error("all multirate timesteps must be positive");
        cell.assign(n, 0);
        for (auto &v : q) v.assign(n, 0.0);
        for (auto &v : qnew) v.assign(n, 0.0);
        for (auto *group : {&F, &B, &P, &Fnew, &Bnew, &Pnew})
            for (auto &v : *group) v.assign(n, 0.0);
        cg_r.assign(n,0.0);cg_p.assign(n,0.0);cg_Ap.assign(n,0.0);
        build_neighbors();
        initialize_circle();
        initialize_fields();
        order.resize(n);
        std::iota(order.begin(), order.end(), 0);
    }

    int index(int x, int y) const {
        x %= cfg.nx; if (x < 0) x += cfg.nx;
        y %= cfg.ny; if (y < 0) y += cfg.ny;
        return y * cfg.nx + x;
    }
    std::pair<int,int> xy(int i) const { return {i % cfg.nx, i / cfg.nx}; }

    void build_neighbors() {
        const int dq[ORI] = {1, 0, -1, -1, 0, 1};
        const int dr[ORI] = {0, 1, 1, 0, -1, -1};
        for (int m = 0; m < ORI; ++m) nb[m].resize(n);
        third.resize(n);
        std::vector<std::pair<int,int>> offsets;
        for (int dr0 = -3; dr0 <= 3; ++dr0)
            for (int dq0 = -3; dq0 <= 3; ++dq0)
                if ((dq0 || dr0) && std::max({std::abs(dq0), std::abs(dr0), std::abs(dq0 + dr0)}) <= 3)
                    offsets.push_back({dq0, dr0});
        if (offsets.size() != 36) throw std::runtime_error("internal third-neighbor construction error");
        for (int i = 0; i < n; ++i) {
            auto [x,y] = xy(i);
            for (int m = 0; m < ORI; ++m) nb[m][i] = index(x + dq[m], y + dr[m]);
            for (int k = 0; k < 36; ++k) third[i][k] = index(x + offsets[k].first, y + offsets[k].second);
        }
    }

    std::pair<double,double> physical(int i) const {
        auto [x,y] = xy(i);
        return {cfg.dx * (x + 0.5 * y), cfg.dx * (std::sqrt(3.0) * 0.5 * y)};
    }

    void initialize_circle() {
        struct Candidate { double d2; int i; };
        std::vector<Candidate> candidates;
        candidates.reserve(n);
        const int cx = cfg.nx / 2, cy = cfg.ny / 2;
        for (int i = 0; i < n; ++i) {
            auto [x,y] = xy(i);
            double dq = x - cx, dr = y - cy;
            double px = dq + 0.5 * dr, py = std::sqrt(3.0) * 0.5 * dr;
            candidates.push_back({px*px + py*py, i});
        }
        std::nth_element(candidates.begin(), candidates.begin() + cfg.initial_sites, candidates.end(),
                         [](const Candidate &a, const Candidate &b){ return a.d2 < b.d2; });
        for (int k = 0; k < cfg.initial_sites; ++k) cell[candidates[k].i] = 1;
        cell_sites = cfg.initial_sites;
        target_area = cfg.target_area_um2>0 ? cfg.target_area_um2*1e6
                    : (cfg.initial_sites == 30000 && std::abs(cfg.dx - 100.0) < 1e-12)
                    ? 3e8 : cfg.initial_sites * cfg.dx * cfg.dx;
        target_sites=target_area/(cfg.dx*cfg.dx);
        rebuild_active();
    }

    void initialize_fields() {
        for (int i : active) {
            q[C][i] = 1.0; q[R][i] = 3.0; q[RH][i] = 1.25;
            q[CI][i] = 1.4; q[RI][i] = 4.5; q[RHI][i] = 1.85; q[ARP][i] = 2.0;
            for (int m = 0; m < ORI; ++m) {
                F[m][i] = 0.278 / 6.0;
                B[m][i] = 1.7e-5 / 6.0;
                P[m][i] = 0.0;
            }
        }
        if(cfg.initial_p_total>0)for(int i:active)if(boundary(i)){
            int exposed=0;for(int m=0;m<ORI;++m)exposed+=!cell[nb[m][i]];
            for(int m=0;m<ORI;++m)if(!cell[nb[m][i]])P[m][i]=cfg.initial_p_total/exposed;
        }
        auto c = centroid(); initial_cx = c.first; initial_cy = c.second;
    }

    void rebuild_active() {
        active.clear(); active.reserve(cell_sites + 1024);
        for (int i = 0; i < n; ++i) if (cell[i]) active.push_back(i);
        if ((int)active.size() != cell_sites) throw std::runtime_error("cell site count mismatch");
    }

    bool boundary(int i) const {
        if (!cell[i]) return false;
        for (int m = 0; m < ORI; ++m) if (!cell[nb[m][i]]) return true;
        return false;
    }

    double laplacian(const std::vector<double> &v, int i) const {
        double s = 0;
        for (int m = 0; m < ORI; ++m) {
            int j = nb[m][i];
            if (cell[j]) s += v[j] - v[i]; // reflected ghost value gives no flux
        }
        return (2.0 / (3.0 * cfg.dx * cfg.dx)) * s;
    }

    double pulse_scale(int i) const {
        if (!cfg.pulse || time >= cfg.pulse_seconds) return 1.0;
        double s = (pulse_hi > pulse_lo) ? (physical(i).first - pulse_lo) / (pulse_hi - pulse_lo) - 0.5 : 0.0;
        return 1.0 + cfg.pulse_fraction * s;
    }

    void field_substep(double dt, bool explicit_diffusion) {
        constexpr double D[NF] = {1e5,1e5,1e5,1e7,1e7,1e7,1e6};
        #pragma omp parallel for schedule(static)
        for (std::size_t ai=0; ai<active.size(); ++ai) {
            int i=active[ai];
            double c=q[C][i], r0=q[R][i], h=q[RH][i];
            double hs=h/1.25, c3=c*c*c;
            double fc = 3.4*pulse_scale(i)/(1.0+hs*hs*hs)*q[CI][i]/2.4-c;
            double fr = (0.5+4.5*c)*q[RI][i]/7.5-r0;
            double fh = (3.3+0.3*r0)/(1.0+c3)*q[RHI][i]/3.1-h;
            double ft=0; for(int m=0;m<ORI;++m) ft+=F[m][i];
            double eta=60.0*q[ARP][i]/(2.0+q[ARP][i]+255.0*ft);
            const double reac[NF] = {fc,fr,fh,-fc,-fr,-fh,0.16*c+0.16*r0-eta*ft-0.1*q[ARP][i]};
            for (int f=0; f<NF; ++f) qnew[f][i]=q[f][i]+dt*(reac[f]+(explicit_diffusion?D[f]*laplacian(q[f],i):0.0));
            double kap=2.8-2.1*r0/(3.0+r0);
            for(int m=0;m<ORI;++m){
                int up=nb[(m+3)%ORI][i], down=nb[m][i];
                double adv=-500.0/cfg.dx*(B[m][i]-(cell[up]?B[m][up]:0.0));
                double branch=0.5*1.06e-4*eta*(F[(m+5)%ORI][i]+F[(m+1)%ORI][i]);
                Fnew[m][i]=F[m][i]+dt*(500.0*B[m][i]-0.03*F[m][i]);
                Bnew[m][i]=B[m][i]+dt*(adv+branch-kap*B[m][i]);
                double hit=cell[down]?0.0:500.0*B[m][i];
                Pnew[m][i]=P[m][i]+dt*(hit-0.14*kap*P[m][i]);
                if (Pnew[m][i] < 0 && Pnew[m][i] > -1e-15) Pnew[m][i]=0;
            }
        }
        for (int f=0;f<NF;++f) q[f].swap(qnew[f]);
        for (int m=0;m<ORI;++m){F[m].swap(Fnew[m]);B[m].swap(Bnew[m]);P[m].swap(Pnew[m]);}
    }

    void apply_implicit_matrix(const std::vector<double>& in,std::vector<double>& out,double alpha) const {
        #pragma omp parallel for schedule(static)
        for(std::size_t ai=0;ai<active.size();++ai){
            int i=active[ai],degree=0;double sum=0;
            for(int m=0;m<ORI;++m){int j=nb[m][i];if(cell[j]){++degree;sum+=in[j];}}
            out[i]=(1.0+alpha*degree)*in[i]-alpha*sum;
        }
    }

    void implicit_diffuse(int field,double dt,double diffusivity){
        if(dt<=0)return;
        double before=0;
        #pragma omp parallel for reduction(+:before) schedule(static)
        for(std::size_t ai=0;ai<active.size();++ai){int i=active[ai];qnew[field][i]=q[field][i];before+=q[field][i];}
        double alpha=2.0*diffusivity*dt/(3.0*cfg.dx*cfg.dx);
        apply_implicit_matrix(q[field],cg_Ap,alpha);
        double rs=0,bs=0;
        #pragma omp parallel for reduction(+:rs,bs) schedule(static)
        for(std::size_t ai=0;ai<active.size();++ai){int i=active[ai];cg_r[i]=qnew[field][i]-cg_Ap[i];cg_p[i]=cg_r[i];rs+=cg_r[i]*cg_r[i];bs+=qnew[field][i]*qnew[field][i];}
        double relative=std::sqrt(rs/std::max(bs,1e-300));int used=0;
        for(;used<cfg.cg_max_iterations&&relative>cfg.cg_tolerance;++used){
            apply_implicit_matrix(cg_p,cg_Ap,alpha);double denom=0;
            #pragma omp parallel for reduction(+:denom) schedule(static)
            for(std::size_t ai=0;ai<active.size();++ai){int i=active[ai];denom+=cg_p[i]*cg_Ap[i];}
            if(!(denom>0))throw std::runtime_error("implicit diffusion CG lost positive definiteness");
            double step=rs/denom,rsnew=0;
            #pragma omp parallel for reduction(+:rsnew) schedule(static)
            for(std::size_t ai=0;ai<active.size();++ai){int i=active[ai];q[field][i]+=step*cg_p[i];cg_r[i]-=step*cg_Ap[i];rsnew+=cg_r[i]*cg_r[i];}
            relative=std::sqrt(rsnew/std::max(bs,1e-300));
            if(relative<=cfg.cg_tolerance){rs=rsnew;++used;break;}
            double beta=rsnew/rs;rs=rsnew;
            #pragma omp parallel for schedule(static)
            for(std::size_t ai=0;ai<active.size();++ai){int i=active[ai];cg_p[i]=cg_r[i]+beta*cg_p[i];}
        }
        if(relative>cfg.cg_tolerance)throw std::runtime_error("implicit diffusion CG failed to converge");
        double after=0;
        #pragma omp parallel for reduction(+:after) schedule(static)
        for(std::size_t ai=0;ai<active.size();++ai)after+=q[field][active[ai]];
        double scale=before/after;
        #pragma omp parallel for schedule(static)
        for(std::size_t ai=0;ai<active.size();++ai)q[field][active[ai]]*=scale;
        ++cg_solves;cg_iterations+=used;cg_max_used=std::max(cg_max_used,used);cg_worst_relative_residual=std::max(cg_worst_relative_residual,relative);
    }

    void evolve_fields() {
        if (cfg.pulse && time < cfg.pulse_seconds) {
            pulse_lo = std::numeric_limits<double>::infinity(); pulse_hi = -pulse_lo;
            for (int i : active) { double x=physical(i).first; pulse_lo=std::min(pulse_lo,x); pulse_hi=std::max(pulse_hi,x); }
        }
        if(!cfg.implicit){
            int sub=(int)std::ceil(cfg.mcs_dt/cfg.pde_dt);double dt=cfg.mcs_dt/sub;
            for(int k=0;k<sub;++k)field_substep(dt,true);
            return;
        }
        constexpr double D[NF]={1e5,1e5,1e5,1e7,1e7,1e7,1e6};
        const double interval[NF]={cfg.active_diffusion_dt,cfg.active_diffusion_dt,cfg.active_diffusion_dt,
                                   cfg.inactive_diffusion_dt,cfg.inactive_diffusion_dt,cfg.inactive_diffusion_dt,cfg.arp_diffusion_dt};
        std::array<double,NF> accumulated{};
        int sub=(int)std::ceil(cfg.mcs_dt/cfg.reaction_dt);double dt=cfg.mcs_dt/sub;
        for(int k=0;k<sub;++k){
            field_substep(dt,false);
            for(int f=0;f<NF;++f){
                accumulated[f]+=dt;
                if(accumulated[f]>=interval[f]-1e-12||k==sub-1){implicit_diffuse(f,accumulated[f],D[f]);accumulated[f]=0;}
            }
        }
    }

    std::array<double,NF> scalar_mass() const {
        std::array<double,NF> z{};
        for(int i:active) for(int f=0;f<NF;++f) z[f]+=q[f][i];
        return z;
    }
    std::array<double,ORI> f_mass() const {
        std::array<double,ORI> z{}; for(int i:active)for(int m=0;m<ORI;++m)z[m]+=F[m][i]*cfg.dx*cfg.dx; return z;
    }
    std::array<double,ORI> bp_mass() const {
        std::array<double,ORI> z{}; for(int i:active)for(int m=0;m<ORI;++m)z[m]+=B[m][i]*cfg.dx*cfg.dx+P[m][i]*cfg.dx; return z;
    }

    bool removable(int i) const {
        std::array<bool,ORI> occ{}; int count=0;
        for(int m=0;m<ORI;++m){occ[m]=cell[nb[m][i]];count+=occ[m];}
        if(count==0)return false;
        int runs=0; for(int m=0;m<ORI;++m) if(occ[m]&&!occ[(m+5)%ORI])++runs;
        return runs<=1;
    }

    double local_contact_delta(int target, bool new_state) const {
        bool old=cell[target]; int before=0,after=0;
        for(int j:third[target]){before+=(old!=bool(cell[j]));after+=(new_state!=bool(cell[j]));}
        return cfg.contact_j*(after-before);
    }

    double push_toward(int site,int dir) const {
        // Dot products of the six unit lattice vectors: 1, 1/2, -1/2, -1, ...
        double sum=0;
        for(int m=0;m<ORI;++m){int d=std::abs(m-dir);d=std::min(d,ORI-d);if(d<2)sum+=P[m][site];}
        return sum;
    }

    int fallback_receiver(int removed,int preferred) const {
        if(cell[preferred])return preferred;
        for(int m=0;m<ORI;++m)if(cell[nb[m][removed]])return nb[m][removed];
        return -1;
    }

    void extend(int target,int source) {
        for(int f=0;f<NF;++f)q[f][target]=q[f][source];
        for(int m=0;m<ORI;++m){F[m][target]=F[m][source];B[m][target]=B[m][source];P[m][target]=0;}
        cell[target]=1;++cell_sites;
        // A membrane segment swallowed by the extension becomes free B again.
        for(int m=0;m<ORI;++m){
            int j=nb[m][target];
            if(cell[j]&&!boundary(j))for(int k=0;k<ORI;++k){B[k][j]+=P[k][j]/cfg.dx;P[k][j]=0;}
        }
    }

    void retract(int target) {
        cell[target]=0;--cell_sites;
        for(int m=0;m<ORI;++m){
            int recv=fallback_receiver(target,nb[(m+3)%ORI][target]);
            if(recv<0)continue;
            F[m][recv]+=F[m][target];
            // If the receiver now points into the removed site, it is a new pushing end.
            if(nb[m][recv]==target)P[m][recv]+=P[m][target]+B[m][target]*cfg.dx;
            else B[m][recv]+=B[m][target]+P[m][target]/cfg.dx;
        }
        for(int f=0;f<NF;++f)q[f][target]=0;
        for(int m=0;m<ORI;++m){F[m][target]=B[m][target]=P[m][target]=0;}
    }

    void renormalize_geometry(const std::array<double,NF>& qm,
                              const std::array<double,ORI>& fm,
                              const std::array<double,ORI>& bpm) {
        rebuild_active();
        auto qa=scalar_mass(); auto fa=f_mass(); auto ba=bp_mass();
        for(int f=0;f<NF;++f){
            if(qa[f]>0){double scale=qm[f]/qa[f];for(int i:active)q[f][i]*=scale;}
            max_scalar_geom_error[f]=std::max(max_scalar_geom_error[f],std::abs(qa[f]-qm[f])/std::max(std::abs(qm[f]),1e-300));
        }
        for(int m=0;m<ORI;++m){
            if(fa[m]>0){double scale=fm[m]/fa[m];for(int i:active)F[m][i]*=scale;}
            if(ba[m]>0){double scale=bpm[m]/ba[m];for(int i:active){B[m][i]*=scale;P[m][i]*=scale;}}
            max_f_geom_error[m]=std::max(max_f_geom_error[m],std::abs(fa[m]-fm[m])/std::max(std::abs(fm[m]),1e-300));
            max_bp_geom_error[m]=std::max(max_bp_geom_error[m],std::abs(ba[m]-bpm[m])/std::max(std::abs(bpm[m]),1e-300));
        }
    }

    double area_energy_delta(int sites,bool extension) const {
        double next=sites+(extension?1.0:-1.0);
        if(cfg.area_mode==AreaMode::CPMLatticeArea)
            return cfg.lambda_cpm*((next-target_sites)*(next-target_sites)-(sites-target_sites)*(sites-target_sites));
        double pixel_area=cfg.dx*cfg.dx,area=sites*pixel_area,anext=next*pixel_area;
        double coefficient=cfg.area_mode==AreaMode::PublishedLiteral ? cfg.lambda_cpm
                          : cfg.lambda_cpm/(pixel_area*pixel_area);
        return coefficient*((anext-target_area)*(anext-target_area)-(area-target_area)*(area-target_area));
    }

    void monte_carlo_step() {
        copydiag=CopyDiagnostics{};
        if(cfg.frozen)return;
        auto qm=scalar_mass(); auto fm=f_mass(); auto bpm=bp_mass();
        std::shuffle(order.begin(),order.end(),rng);
        for(int target:order){
            int dir=(int)(uniform(rng)*ORI); if(dir==ORI)dir=ORI-1;
            int source=nb[dir][target]; bool old=cell[target], neu=cell[source];
            if(old==neu)continue;
            ++copydiag.candidates;
            if(old&&(!removable(target)||cell_sites<=7)){++copydiag.connectivity_rejects;continue;}
            double contact=local_contact_delta(target,neu);
            double area_energy=area_energy_delta(cell_sites,neu);
            int site=old?target:source, move_dir=old?dir:(dir+3)%ORI;
            double pushing=push_toward(site,move_dir);
            double contraction=cfg.xi*std::max(0.0,q[RH][site]-1.25);
            double rho_energy=neu?contraction:-contraction, push_energy=neu?-pushing:pushing;
            double dH=contact+area_energy+push_energy+rho_energy;
            double probability=dH < -cfg.yield_energy ? 1.0 : std::min(1.0,std::exp(-(dH+cfg.yield_energy)/cfg.temperature));
            if(neu){++copydiag.ext_trials;copydiag.ext_expected+=probability;copydiag.ext_min_dh=std::min(copydiag.ext_min_dh,dH);
                copydiag.ext_area+=area_energy;copydiag.ext_contact+=contact;copydiag.ext_push+=pushing;copydiag.ext_rho+=contraction;}
            else{++copydiag.ret_trials;copydiag.ret_expected+=probability;copydiag.ret_min_dh=std::min(copydiag.ret_min_dh,dH);
                copydiag.ret_area+=area_energy;copydiag.ret_contact+=contact;copydiag.ret_push+=pushing;copydiag.ret_rho+=contraction;}
            bool accept=probability>=1.0||uniform(rng)<probability;
            if(accept){if(neu){++copydiag.ext_accepted;extend(target,source);}else{++copydiag.ret_accepted;retract(target);}}
        }
        renormalize_geometry(qm,fm,bpm);
    }

    std::pair<double,double> centroid() const {
        double x=0,y=0;for(int i:active){auto p=physical(i);x+=p.first;y+=p.second;}return{x/active.size(),y/active.size()};
    }
    double perimeter_nm() const {double z=0;for(int i:active)for(int m=0;m<ORI;++m)if(!cell[nb[m][i]])z+=cfg.dx/std::sqrt(3.0);return z;}
    double weighted_x(const std::vector<double>&v)const{double a=0,b=0;for(int i:active){a+=physical(i).first*v[i];b+=v[i];}return b?a/b:0;}

    void open_outputs() {
        fs::create_directories(cfg.output);
        metrics.open(fs::path(cfg.output)/"metrics.csv");
        if(cfg.snapshots&&!cfg.binary_snapshots)snapshots.open(fs::path(cfg.output)/"snapshots.csv");
        metrics<<std::setprecision(12)<<"step,time_s,sites,area_um2,perimeter_um,cx_um,cy_um,displacement_um,C_mass,R_mass,rho_mass,Ci_mass,Ri_mass,rhoi_mass,A_mass,C_minus_rho_x_um,copy_candidates,connectivity_rejects,extension_trials,retraction_trials,extensions,retractions,expected_extensions,expected_retractions,min_extension_dH,min_retraction_dH,mean_extension_area_dH,mean_retraction_area_dH,mean_extension_contact_dH,mean_retraction_contact_dH,mean_extension_push,mean_retraction_push,mean_extension_rho,mean_retraction_rho\n";
        if(cfg.snapshots)snapshots<<"step,time_s,q,r,x_um,y_um,C,R,rho,Ci,Ri,rhoi,A,F_total,B_total,P_total\n";
    }
    void write_output(int step,bool snap) {
        auto c=centroid();auto mass=scalar_mass();double disp=std::hypot(c.first-initial_cx,c.second-initial_cy)/1000;
        metrics<<step<<','<<time<<','<<cell_sites<<','<<cell_sites*cfg.dx*cfg.dx/1e6<<','<<perimeter_nm()/1000<<','<<c.first/1000<<','<<c.second/1000<<','<<disp;
        for(double x:mass) metrics<<','<<x;
        auto avg=[](double x,long count){return count?x/count:0.0;};
        auto finite=[](double x){return std::isfinite(x)?x:std::numeric_limits<double>::quiet_NaN();};
        metrics<<','<<(weighted_x(q[C])-weighted_x(q[RH]))/1000
               <<','<<copydiag.candidates<<','<<copydiag.connectivity_rejects<<','<<copydiag.ext_trials<<','<<copydiag.ret_trials
               <<','<<copydiag.ext_accepted<<','<<copydiag.ret_accepted<<','<<copydiag.ext_expected<<','<<copydiag.ret_expected
               <<','<<finite(copydiag.ext_min_dh)<<','<<finite(copydiag.ret_min_dh)
               <<','<<avg(copydiag.ext_area,copydiag.ext_trials)<<','<<avg(copydiag.ret_area,copydiag.ret_trials)
               <<','<<avg(copydiag.ext_contact,copydiag.ext_trials)<<','<<avg(copydiag.ret_contact,copydiag.ret_trials)
               <<','<<avg(copydiag.ext_push,copydiag.ext_trials)<<','<<avg(copydiag.ret_push,copydiag.ret_trials)
               <<','<<avg(copydiag.ext_rho,copydiag.ext_trials)<<','<<avg(copydiag.ret_rho,copydiag.ret_trials)<<'\n';
        metrics.flush();
        if(snap&&cfg.snapshots&&cfg.binary_snapshots){
            std::ofstream bin(fs::path(cfg.output)/("frame_"+std::to_string(step)+".bin"),std::ios::binary);
            int32_t header[3]={step,cell_sites,25};bin.write(reinterpret_cast<char*>(header),sizeof(header));
            for(int i:active){int32_t id=i;bin.write(reinterpret_cast<char*>(&id),sizeof(id));}
            for(int f=0;f<NF;++f)for(int i:active){double v=q[f][i];bin.write(reinterpret_cast<char*>(&v),sizeof(v));}
            for(auto group:{&F,&B,&P})for(int m=0;m<ORI;++m)for(int i:active){double v=(*group)[m][i];bin.write(reinterpret_cast<char*>(&v),sizeof(v));}
            if(!bin)throw std::runtime_error("snapshot write failed");
        }
        if(snap&&cfg.snapshots&&!cfg.binary_snapshots){for(int i:active){auto [x,y]=xy(i);auto p=physical(i);double ft=0,bt=0,pt=0;for(int m=0;m<ORI;++m){ft+=F[m][i];bt+=B[m][i];pt+=P[m][i];}
            snapshots<<step<<','<<time<<','<<x<<','<<y<<','<<p.first/1000<<','<<p.second/1000;
            for(int f=0;f<NF;++f) snapshots<<','<<q[f][i];
            snapshots<<','<<ft<<','<<bt<<','<<pt<<'\n';
        }
        snapshots.flush();}
    }

    void write_manifest() const {
        std::ofstream o(fs::path(cfg.output)/"manifest.txt");
        o<<std::setprecision(12)
         <<"implementation=Maree_2006_custom_cpp\nseed="<<cfg.seed<<"\nnx="<<cfg.nx<<"\nny="<<cfg.ny
         <<"\ninitial_sites="<<cfg.initial_sites<<"\nsteps="<<cfg.steps<<"\ndx_nm="<<cfg.dx
         <<"\nmcs_dt_s="<<cfg.mcs_dt<<"\nintegrator="<<(cfg.implicit?"implicit_multirate":"explicit_reference")
         <<"\npde_dt_s="<<cfg.pde_dt<<"\nreaction_dt_s="<<cfg.reaction_dt
         <<"\nactive_diffusion_dt_s="<<cfg.active_diffusion_dt<<"\narp_diffusion_dt_s="<<cfg.arp_diffusion_dt
         <<"\ninactive_diffusion_dt_s="<<cfg.inactive_diffusion_dt<<"\ncg_tolerance="<<cfg.cg_tolerance
         <<"\ncg_max_iterations="<<cfg.cg_max_iterations<<"\ncg_solves="<<cg_solves
         <<"\ncg_mean_iterations="<<(cg_solves?double(cg_iterations)/cg_solves:0)<<"\ncg_max_used="<<cg_max_used
         <<"\ncg_worst_relative_residual="<<cg_worst_relative_residual<<"\npulse="<<cfg.pulse
         <<"\npulse_seconds="<<cfg.pulse_seconds<<"\npulse_fraction="<<cfg.pulse_fraction
         <<"\nxi="<<cfg.xi<<"\ntarget_area_um2="<<target_area/1e6<<"\ntarget_sites="<<target_sites
         <<"\narea_mode="<<area_mode_name(cfg.area_mode)<<"\nlambda_cpm="<<cfg.lambda_cpm
         <<"\nlambda_physical_equivalent="<<cfg.lambda_cpm/std::pow(cfg.dx*cfg.dx,2)
         <<"\ncontact_j="<<cfg.contact_j<<"\ntemperature="<<cfg.temperature<<"\nyield_energy="<<cfg.yield_energy
         <<"\ninitial_p_total="<<cfg.initial_p_total
         <<"\nfrozen="<<cfg.frozen<<'\n';
        for(int f=0;f<NF;++f)o<<"max_pre_correction_scalar_"<<f<<'='<<max_scalar_geom_error[f]<<'\n';
        for(int m=0;m<ORI;++m)o<<"max_pre_correction_F_"<<m<<'='<<max_f_geom_error[m]<<"\nmax_pre_correction_BP_"<<m<<'='<<max_bp_geom_error[m]<<'\n';
    }

    void run() {
        open_outputs(); write_output(0,true);
        int snap_every=std::max(1,(int)std::llround(cfg.snapshot_seconds/cfg.mcs_dt));
        for(int step=1;step<=cfg.steps;++step){evolve_fields();monte_carlo_step();time=step*cfg.mcs_dt;write_output(step,step%snap_every==0||step==cfg.steps);}
        write_manifest();
    }
};

Config parse(int argc,char**argv){
    Config c;
    auto value=[&](int&i){if(i+1>=argc)throw std::runtime_error(std::string("missing value for ")+argv[i]);return std::string(argv[++i]);};
    for(int i=1;i<argc;++i){std::string a=argv[i];
        if(a=="--steps")c.steps=std::stoi(value(i));else if(a=="--seed")c.seed=std::stoull(value(i));
        else if(a=="--output")c.output=value(i);else if(a=="--nx")c.nx=std::stoi(value(i));else if(a=="--ny")c.ny=std::stoi(value(i));
        else if(a=="--initial-sites")c.initial_sites=std::stoi(value(i));else if(a=="--pde-dt")c.pde_dt=std::stod(value(i));
        else if(a=="--reaction-dt")c.reaction_dt=std::stod(value(i));else if(a=="--active-diffusion-dt")c.active_diffusion_dt=std::stod(value(i));
        else if(a=="--arp-diffusion-dt")c.arp_diffusion_dt=std::stod(value(i));else if(a=="--inactive-diffusion-dt")c.inactive_diffusion_dt=std::stod(value(i));
        else if(a=="--cg-tolerance")c.cg_tolerance=std::stod(value(i));else if(a=="--cg-max-iterations")c.cg_max_iterations=std::stoi(value(i));
        else if(a=="--integrator"){std::string mode=value(i);if(mode=="implicit")c.implicit=true;else if(mode=="explicit-reference")c.implicit=false;else throw std::runtime_error("integrator must be implicit or explicit-reference");}
        else if(a=="--pulse-fraction")c.pulse_fraction=std::stod(value(i));else if(a=="--pulse-seconds")c.pulse_seconds=std::stod(value(i));
        else if(a=="--xi")c.xi=std::stod(value(i));else if(a=="--snapshot-seconds")c.snapshot_seconds=std::stod(value(i));
        else if(a=="--target-area-um2")c.target_area_um2=std::stod(value(i));
        else if(a=="--lambda-cpm"||a=="--area-lambda")c.lambda_cpm=std::stod(value(i));
        else if(a=="--area-mode"){std::string mode=value(i);
            if(mode=="published_literal")c.area_mode=AreaMode::PublishedLiteral;
            else if(mode=="cpm_lattice_area")c.area_mode=AreaMode::CPMLatticeArea;
            else if(mode=="physical_area_equivalent")c.area_mode=AreaMode::PhysicalAreaEquivalent;
            else throw std::runtime_error("area mode must be published_literal, cpm_lattice_area, or physical_area_equivalent");}
        else if(a=="--contact-j")c.contact_j=std::stod(value(i));else if(a=="--temperature")c.temperature=std::stod(value(i));
        else if(a=="--yield-energy")c.yield_energy=std::stod(value(i));
        else if(a=="--initial-p-total")c.initial_p_total=std::stod(value(i));
        else if(a=="--no-pulse")c.pulse=false;else if(a=="--frozen")c.frozen=true;else if(a=="--no-snapshots")c.snapshots=false;
        else if(a=="--binary-snapshots")c.binary_snapshots=true;
        else if(a=="--self-test"){}else if(a=="--help"){std::cout<<"See README.md for usage.\n";std::exit(0);}
        else throw std::runtime_error("unknown argument: "+a);
    }return c;
}

void self_test(){
    double e=0;
    for(bool implicit:{false,true}){
        Config c;c.nx=48;c.ny=48;c.initial_sites=400;c.steps=0;c.pulse=false;c.frozen=true;c.snapshots=false;c.implicit=implicit;
        Model m(c);auto before=m.scalar_mass();m.evolve_fields();auto after=m.scalar_mass();
        for(int p:{C,R,RH})e=std::max(e,std::abs((after[p]+after[p+3])-(before[p]+before[p+3]))/std::max(1.0,before[p]+before[p+3]));
    }
    if(e>2e-12)throw std::runtime_error("GTPase pool conservation self-test failed");
    Config ac;ac.nx=48;ac.ny=48;ac.initial_sites=400;ac.pulse=false;ac.frozen=true;ac.snapshots=false;
    ac.area_mode=AreaMode::CPMLatticeArea;Model lattice(ac);
    ac.area_mode=AreaMode::PhysicalAreaEquivalent;Model converted(ac);
    for(int offset:{-100,-1,0,1,100})for(bool extension:{false,true}){
        double a=lattice.area_energy_delta(400+offset,extension),b=converted.area_energy_delta(400+offset,extension);
        if(std::abs(a-b)>1e-14*std::max({1.0,std::abs(a),std::abs(b)}))throw std::runtime_error("area-mode equivalence self-test failed");
    }
    ac.area_mode=AreaMode::PublishedLiteral;Model literal(ac);
    if(std::abs(literal.area_energy_delta(400,true)-0.1)>1e-12)throw std::runtime_error("published-literal area self-test failed");
    double act=std::abs(500*1.7e-5-0.03*0.278);
    if(act>2e-4)throw std::runtime_error("actin arithmetic self-test failed");
    std::cout<<"self-test passed; max GTPase relative error="<<std::setprecision(4)<<e<<"\n";
}

int main(int argc,char**argv){
    try{
        if (!std::getenv("OMP_NUM_THREADS")) omp_set_num_threads(1);
        for(int i=1;i<argc;++i)if(std::string(argv[i])=="--self-test"){self_test();return 0;}
        Config c=parse(argc,argv);Model m(c);m.run();
    }catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}return 0;
}
