// CUDA field-resident variant of the Maree et al. 2006 reconstruction.
// State remains on device; each MCS transfers only the shuffled proposal order
// and compact diagnostics. See GPU.md for numerical and benchmarking details.
#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs=std::filesystem;
constexpr int NX=500,NY=500,N=NX*NY,NF=7,NO=6;
enum Field{C,R,RH,CI,RI,RHI,ARP};

#define CUDA_OK(x) do{cudaError_t e=(x);if(e!=cudaSuccess)throw std::runtime_error(std::string(#x)+": "+cudaGetErrorString(e));}while(0)

struct Config{int steps=1200;uint64_t seed=1;std::string output="run_gpu";double lambda_cpm=1e-9;bool pulse=true,snapshot=false,frozen=false;};
struct CopyDiag{long long candidates,connectivity_rejects,ext_trials,ret_trials,ext_accepted,ret_accepted;double ext_expected,ret_expected;};
struct Totals{double q[NF],f[NO],bp[NO],x,y,perimeter;int sites;};

__device__ inline int ni(const int* nb,int m,int i){return nb[m*N+i];}
__device__ inline double cube(double x){return x*x*x;}

__global__ void field_step(const unsigned char* cell,const int* active,const int* sites,const int* nb,const double* q,double* qn,
                           const double* F,const double* B,const double* P,double* Fn,double* Bn,double* Pn,
                           double dt,double time,double pulse_lo,double pulse_hi,bool pulse){
 int ai=blockIdx.x*blockDim.x+threadIdx.x;if(ai>=*sites)return;int i=active[ai];
 double c=q[C*N+i],r=q[R*N+i],h=q[RH*N+i];int x=i%NX,y=i/NX;
 double px=100.0*(x+0.5*y),s=(pulse&&time<10.0)?1.0+0.15*((px-pulse_lo)/(pulse_hi-pulse_lo)-0.5):1.0;
 double fc=3.4*s/(1.0+cube(h/1.25))*q[CI*N+i]/2.4-c;
 double fr=(0.5+4.5*c)*q[RI*N+i]/7.5-r;
 double fh=(3.3+0.3*r)/(1.0+cube(c))*q[RHI*N+i]/3.1-h;
 double ft=0;for(int m=0;m<NO;++m)ft+=F[m*N+i];
 double eta=60.0*q[ARP*N+i]/(2.0+q[ARP*N+i]+255.0*ft);
 double reac[NF]={fc,fr,fh,-fc,-fr,-fh,0.16*c+0.16*r-eta*ft-0.1*q[ARP*N+i]};
 for(int f=0;f<NF;++f)qn[f*N+i]=q[f*N+i]+dt*reac[f];
 double kap=2.8-2.1*r/(3.0+r);
 for(int m=0;m<NO;++m){int up=ni(nb,(m+3)%NO,i),down=ni(nb,m,i);
  double adv=-5.0*(B[m*N+i]-(cell[up]?B[m*N+up]:0.0));
  double branch=0.5*1.06e-4*eta*(F[((m+5)%NO)*N+i]+F[((m+1)%NO)*N+i]);
  Fn[m*N+i]=F[m*N+i]+dt*(500.0*B[m*N+i]-0.03*F[m*N+i]);
  Bn[m*N+i]=B[m*N+i]+dt*(adv+branch-kap*B[m*N+i]);
  double hit=cell[down]?0.0:500.0*B[m*N+i];
  Pn[m*N+i]=fmax(0.0,P[m*N+i]+dt*(hit-0.14*kap*P[m*N+i]));
 }
}

__global__ void copy_field(const double* q,double* b,int f,const unsigned char* cell){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<N&&cell[i])b[i]=q[f*N+i];}
__global__ void matvec(const double* in,double* out,const unsigned char* cell,const int* nb,double alpha){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=N||!cell[i])return;int degree=0;double sum=0;
 for(int m=0;m<NO;++m){int j=ni(nb,m,i);if(cell[j]){++degree;sum+=in[j];}}out[i]=(1+alpha*degree)*in[i]-alpha*sum;
}
__global__ void init_cg(const double* b,const double* Ap,double* r,double* p,const unsigned char* cell,double* sums){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=N||!cell[i])return;double z=b[i]-Ap[i];r[i]=p[i]=z;atomicAdd(&sums[0],z*z);atomicAdd(&sums[1],b[i]*b[i]);
}
__global__ void dot_kernel(const double* a,const double* b,const unsigned char* cell,double* out){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<N&&cell[i])atomicAdd(out,a[i]*b[i]);}
__global__ void cg_update(double* x,double* r,const double* p,const double* Ap,const unsigned char* cell,double step,double* out){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<N&&cell[i]){x[i]+=step*p[i];r[i]-=step*Ap[i];atomicAdd(out,r[i]*r[i]);}}
__global__ void cg_direction(double* p,const double* r,const unsigned char* cell,double beta){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<N&&cell[i])p[i]=r[i]+beta*p[i];}

__global__ void cg_cooperative(double* x,double* b,double* r,double* p,double* Ap,const unsigned char* cell,const int* active,const int* sites,const int* nb,
                               double alpha,double tolerance,int max_iterations,double* work,int* keep_going){
 namespace cg=cooperative_groups;cg::grid_group grid=cg::this_grid();int tid=blockIdx.x*blockDim.x+threadIdx.x,stride=gridDim.x*blockDim.x;
 __shared__ double sh[512];
 double local_rs=0,local_bs=0;
 for(int ai=tid;ai<*sites;ai+=stride){int i=active[ai];b[i]=x[i];int degree=0;double sum=0;for(int m=0;m<NO;++m){int j=ni(nb,m,i);if(cell[j]){++degree;sum+=x[j];}}Ap[i]=(1+alpha*degree)*x[i]-alpha*sum;double z=b[i]-Ap[i];r[i]=p[i]=z;local_rs+=z*z;local_bs+=b[i]*b[i];}
 sh[threadIdx.x]=local_rs;sh[256+threadIdx.x]=local_bs;__syncthreads();for(int d=128;d;d>>=1){if(threadIdx.x<d){sh[threadIdx.x]+=sh[threadIdx.x+d];sh[256+threadIdx.x]+=sh[256+threadIdx.x+d];}__syncthreads();}if(threadIdx.x==0){atomicAdd(&work[0],sh[0]);atomicAdd(&work[1],sh[256]);}grid.sync();
 for(int it=0;it<max_iterations;++it){
  if(tid==0){*keep_going=(sqrt(work[0]/fmax(work[1],1e-300))>tolerance);work[2]=0;}grid.sync();if(!*keep_going)break;
  double denom=0;for(int ai=tid;ai<*sites;ai+=stride){int i=active[ai];int degree=0;double sum=0;for(int m=0;m<NO;++m){int j=ni(nb,m,i);if(cell[j]){++degree;sum+=p[j];}}Ap[i]=(1+alpha*degree)*p[i]-alpha*sum;denom+=p[i]*Ap[i];}sh[threadIdx.x]=denom;__syncthreads();for(int d=128;d;d>>=1){if(threadIdx.x<d)sh[threadIdx.x]+=sh[threadIdx.x+d];__syncthreads();}if(threadIdx.x==0)atomicAdd(&work[2],sh[0]);grid.sync();
  if(tid==0){work[3]=work[0]/work[2];work[2]=0;}grid.sync();double step=work[3],rn=0;
  for(int ai=tid;ai<*sites;ai+=stride){int i=active[ai];x[i]+=step*p[i];r[i]-=step*Ap[i];rn+=r[i]*r[i];}sh[threadIdx.x]=rn;__syncthreads();for(int d=128;d;d>>=1){if(threadIdx.x<d)sh[threadIdx.x]+=sh[threadIdx.x+d];__syncthreads();}if(threadIdx.x==0)atomicAdd(&work[2],sh[0]);grid.sync();
  if(tid==0){work[3]=work[2]/work[0];work[0]=work[2];}grid.sync();double beta=work[3];for(int ai=tid;ai<*sites;ai+=stride){int i=active[ai];p[i]=r[i]+beta*p[i];}grid.sync();
 }
}

__device__ inline uint64_t xr(uint64_t& s){s^=s>>12;s^=s<<25;s^=s>>27;return s*2685821657736338717ULL;}
__device__ inline double urand(uint64_t& s){return (xr(s)>>11)*0x1.0p-53;}
__device__ bool boundary(const unsigned char* cell,const int* nb,int i){if(!cell[i])return false;for(int m=0;m<NO;++m)if(!cell[ni(nb,m,i)])return true;return false;}
__device__ bool removable(const unsigned char* cell,const int* nb,int i){bool o[NO];int count=0,runs=0;for(int m=0;m<NO;++m){o[m]=cell[ni(nb,m,i)];count+=o[m];}if(!count)return false;for(int m=0;m<NO;++m)if(o[m]&&!o[(m+5)%NO])++runs;return runs<=1;}
__device__ double push_toward(const double* P,int site,int dir){double z=0;for(int m=0;m<NO;++m){int d=abs(m-dir);d=min(d,NO-d);if(d<2)z+=P[m*N+site];}return z;}
__device__ int fallback(const unsigned char* cell,const int* nb,int removed,int preferred){if(cell[preferred])return preferred;for(int m=0;m<NO;++m){int j=ni(nb,m,removed);if(cell[j])return j;}return-1;}

__global__ void cpm_serial(unsigned char* cell,int* active,int* position,const int* nb,const int* third,const int* order,
                           double* q,double* F,double* B,double* P,int* sites,double lambda,uint64_t seed,CopyDiag* d){
 if(blockIdx.x||threadIdx.x)return;CopyDiag z{};uint64_t state=seed;constexpr double J=7.5e-4,T=.008,Hb=.046,xi=.06;
 for(int oi=0;oi<N;++oi){int target=order[oi],dir=int(urand(state)*NO);if(dir>=NO)dir=NO-1;int source=ni(nb,dir,target);bool old=cell[target],neu=cell[source];if(old==neu)continue;++z.candidates;
  if(old&&(!removable(cell,nb,target)||*sites<=7)){++z.connectivity_rejects;continue;}
  int before=0,after=0;for(int k=0;k<36;++k){int j=third[k*N+target];before+=(old!=bool(cell[j]));after+=(neu!=bool(cell[j]));}
  double contact=J*(after-before),delta=*sites-30000.0,next=delta+(neu?1.0:-1.0);double area=lambda*(next*next-delta*delta);
  int site=old?target:source,move=old?dir:(dir+3)%NO;double pushing=push_toward(P,site,move),contraction=xi*fmax(0.0,q[RH*N+site]-1.25);
  double dh=contact+area+(neu?-pushing:pushing)+(neu?contraction:-contraction);double probability=dh<-Hb?1.0:fmin(1.0,exp(-(dh+Hb)/T));
  if(neu){++z.ext_trials;z.ext_expected+=probability;}else{++z.ret_trials;z.ret_expected+=probability;}
  if(probability<1&&urand(state)>=probability)continue;
  if(neu){++z.ext_accepted;for(int f=0;f<NF;++f)q[f*N+target]=q[f*N+source];for(int m=0;m<NO;++m){F[m*N+target]=F[m*N+source];B[m*N+target]=B[m*N+source];P[m*N+target]=0;}cell[target]=1;position[target]=*sites;active[*sites]=target;++*sites;
   for(int m=0;m<NO;++m){int j=ni(nb,m,target);if(cell[j]&&!boundary(cell,nb,j))for(int k=0;k<NO;++k){B[k*N+j]+=P[k*N+j]/100.;P[k*N+j]=0;}}
  }else{++z.ret_accepted;cell[target]=0;int slot=position[target],last=active[*sites-1];active[slot]=last;position[last]=slot;position[target]=-1;--*sites;for(int m=0;m<NO;++m){int recv=fallback(cell,nb,target,ni(nb,(m+3)%NO,target));if(recv<0)continue;F[m*N+recv]+=F[m*N+target];if(ni(nb,m,recv)==target)P[m*N+recv]+=P[m*N+target]+100*B[m*N+target];else B[m*N+recv]+=B[m*N+target]+P[m*N+target]/100.;}
   for(int f=0;f<NF;++f)q[f*N+target]=0;for(int m=0;m<NO;++m)F[m*N+target]=B[m*N+target]=P[m*N+target]=0;
  }
 }*d=z;
}

__global__ void totals_kernel(const unsigned char* cell,const int* active,const int* sites,const int* nb,const double* q,const double* F,const double* B,const double* P,Totals* t){
 int ai=blockIdx.x*blockDim.x+threadIdx.x;if(ai>=*sites)return;int i=active[ai],x=i%NX,y=i/NX;atomicAdd(&t->sites,1);atomicAdd(&t->x,100.0*(x+.5*y));atomicAdd(&t->y,86.60254037844386*y);
 for(int f=0;f<NF;++f)atomicAdd(&t->q[f],q[f*N+i]);for(int m=0;m<NO;++m){atomicAdd(&t->f[m],F[m*N+i]*1e4);atomicAdd(&t->bp[m],B[m*N+i]*1e4+P[m*N+i]*100);if(!cell[ni(nb,m,i)])atomicAdd(&t->perimeter,57.73502691896258);}
}
__global__ void scale_fields(const int* active,const int* sites,double* q,double* F,double* B,double* P,const Totals* before,const Totals* after){int ai=blockIdx.x*blockDim.x+threadIdx.x;if(ai>=*sites)return;int i=active[ai];for(int f=0;f<NF;++f)q[f*N+i]*=before->q[f]/after->q[f];for(int m=0;m<NO;++m){F[m*N+i]*=before->f[m]/after->f[m];double s=before->bp[m]/after->bp[m];B[m*N+i]*=s;P[m*N+i]*=s;}}

struct GPUModel{
 Config cfg;unsigned char *cell;int *active,*position,*nb,*third,*order,*sites,*cg_continue;double *q,*qn,*F,*Fn,*B,*Bn,*P,*Pn,*b,*r,*p,*Ap,*sums;CopyDiag* diag;Totals *before,*after;std::vector<int> horder;std::mt19937_64 rng;double pulse_lo=0,pulse_hi=1;int blocks=(N+255)/256,active_blocks=(50000+255)/256;
 GPUModel(Config c):cfg(c),horder(N),rng(c.seed){
  struct Candidate{double d2;int i;};std::vector<unsigned char> hc(N);std::vector<int> hnb(NO*N),hthird(36*N);std::vector<Candidate> cand;cand.reserve(N);int cx=NX/2,cy=NY/2;
  int dq[NO]={1,0,-1,-1,0,1},dr[NO]={0,1,1,0,-1,-1};auto idx=[](int x,int y){x=(x%NX+NX)%NX;y=(y%NY+NY)%NY;return y*NX+x;};
  std::vector<std::pair<int,int>> off;for(int rr=-3;rr<=3;++rr)for(int qq=-3;qq<=3;++qq)if((qq||rr)&&std::max({abs(qq),abs(rr),abs(qq+rr)})<=3)off.push_back({qq,rr});
  for(int i=0;i<N;++i){int x=i%NX,y=i/NX;double qq=x-cx,rr=y-cy,px=qq+.5*rr,py=std::sqrt(3.0)*.5*rr;cand.push_back({px*px+py*py,i});for(int m=0;m<NO;++m)hnb[m*N+i]=idx(x+dq[m],y+dr[m]);for(int k=0;k<36;++k)hthird[k*N+i]=idx(x+off[k].first,y+off[k].second);}
  std::nth_element(cand.begin(),cand.begin()+30000,cand.end(),[](const Candidate&a,const Candidate&b){return a.d2<b.d2;});for(int k=0;k<30000;++k)hc[cand[k].i]=1;std::iota(horder.begin(),horder.end(),0);
  std::vector<int> hactive;hactive.reserve(50000);std::vector<int> hposition(N,-1);for(int i=0;i<N;++i)if(hc[i]){hposition[i]=hactive.size();hactive.push_back(i);}
  std::vector<double> hq(NF*N),hF(NO*N),hB(NO*N),hP(NO*N);pulse_lo=1e300;pulse_hi=-1e300;
  for(int i=0;i<N;++i)if(hc[i]){int x=i%NX,y=i/NX;double px=100*(x+.5*y);pulse_lo=std::min(pulse_lo,px);pulse_hi=std::max(pulse_hi,px);double vals[NF]={1,3,1.25,1.4,4.5,1.85,2};for(int f=0;f<NF;++f)hq[f*N+i]=vals[f];for(int m=0;m<NO;++m){hF[m*N+i]=.278/6;hB[m*N+i]=1.7e-5/6;}}
  alloc(cell,N);alloc(active,50000);alloc(position,N);alloc(nb,NO*N);alloc(third,36*N);alloc(order,N);alloc(sites,1);alloc(cg_continue,1);alloc(q,NF*N);alloc(qn,NF*N);alloc(F,NO*N);alloc(Fn,NO*N);alloc(B,NO*N);alloc(Bn,NO*N);alloc(P,NO*N);alloc(Pn,NO*N);alloc(b,N);alloc(r,N);alloc(p,N);alloc(Ap,N);alloc(sums,4);alloc(diag,1);alloc(before,1);alloc(after,1);
  put(cell,hc.data(),N);put(active,hactive.data(),hactive.size());put(position,hposition.data(),N);put(nb,hnb.data(),NO*N);put(third,hthird.data(),36*N);put(q,hq.data(),NF*N);put(F,hF.data(),NO*N);put(B,hB.data(),NO*N);put(P,hP.data(),NO*N);int ns=30000;put(sites,&ns,1);CUDA_OK(cudaMemset(qn,0,NF*N*sizeof(double)));CUDA_OK(cudaMemset(Fn,0,NO*N*sizeof(double)));CUDA_OK(cudaMemset(Bn,0,NO*N*sizeof(double)));CUDA_OK(cudaMemset(Pn,0,NO*N*sizeof(double)));
 }
 template<class T>void alloc(T*&x,size_t n){CUDA_OK(cudaMalloc(&x,n*sizeof(T)));}template<class T>void put(T*x,const T*h,size_t n){CUDA_OK(cudaMemcpy(x,h,n*sizeof(T),cudaMemcpyHostToDevice));}
 ~GPUModel(){void* ptrs[]={(void*)cell,(void*)active,(void*)position,(void*)nb,(void*)third,(void*)order,(void*)sites,(void*)cg_continue,(void*)q,(void*)qn,(void*)F,(void*)Fn,(void*)B,(void*)Bn,(void*)P,(void*)Pn,(void*)b,(void*)r,(void*)p,(void*)Ap,(void*)sums,(void*)diag,(void*)before,(void*)after};for(void* x:ptrs)cudaFree(x);}
 double scalar(int j){double v;CUDA_OK(cudaMemcpy(&v,sums+j,sizeof(double),cudaMemcpyDeviceToHost));return v;}
 void diffuse(int f,double dt,double D){double alpha=2*D*dt/(3e4),tol=1e-9;int maxit=200;double* x=q+f*N;CUDA_OK(cudaMemset(sums,0,4*sizeof(double)));void* args[14]={&x,&b,&r,&p,&Ap,&cell,&active,&sites,&nb,&alpha,&tol,&maxit,&sums,&cg_continue};CUDA_OK(cudaLaunchCooperativeKernel((void*)cg_cooperative,32,256,args));}
 Totals totals(Totals* d){CUDA_OK(cudaMemset(d,0,sizeof(Totals)));totals_kernel<<<active_blocks,256>>>(cell,active,sites,nb,q,F,B,P,d);Totals h{};CUDA_OK(cudaMemcpy(&h,d,sizeof(h),cudaMemcpyDeviceToHost));return h;}
 void fields(double time){constexpr double D[NF]={1e5,1e5,1e5,1e7,1e7,1e7,1e6};double accum[NF]={};for(int k=0;k<50;++k){field_step<<<active_blocks,256>>>(cell,active,sites,nb,q,qn,F,B,P,Fn,Bn,Pn,.002,time+k*.002,pulse_lo,pulse_hi,cfg.pulse);std::swap(q,qn);std::swap(F,Fn);std::swap(B,Bn);std::swap(P,Pn);for(int f=0;f<NF;++f){double interval=f<3?.02:.01;accum[f]+=.002;if(accum[f]>=interval-1e-12||k==49){diffuse(f,accum[f],D[f]);accum[f]=0;}}}}
 CopyDiag cpm(int step){Totals hb=totals(before);std::shuffle(horder.begin(),horder.end(),rng);put(order,horder.data(),N);cpm_serial<<<1,1>>>(cell,active,position,nb,third,order,q,F,B,P,sites,cfg.lambda_cpm,cfg.seed+0x9e3779b97f4a7c15ULL*step,diag);Totals ha=totals(after);scale_fields<<<active_blocks,256>>>(active,sites,q,F,B,P,before,after);CopyDiag hd{};CUDA_OK(cudaMemcpy(&hd,diag,sizeof(hd),cudaMemcpyDeviceToHost));return hd;}
 void snapshot(){std::vector<unsigned char> hc(N);std::vector<double> hq(NF*N),hF(NO*N),hB(NO*N),hP(NO*N);CUDA_OK(cudaMemcpy(hc.data(),cell,N,cudaMemcpyDeviceToHost));CUDA_OK(cudaMemcpy(hq.data(),q,NF*N*sizeof(double),cudaMemcpyDeviceToHost));CUDA_OK(cudaMemcpy(hF.data(),F,NO*N*sizeof(double),cudaMemcpyDeviceToHost));CUDA_OK(cudaMemcpy(hB.data(),B,NO*N*sizeof(double),cudaMemcpyDeviceToHost));CUDA_OK(cudaMemcpy(hP.data(),P,NO*N*sizeof(double),cudaMemcpyDeviceToHost));std::ofstream o(fs::path(cfg.output)/"snapshots.csv");o<<"step,time_s,q,r,x_um,y_um,C,R,rho,Ci,Ri,rhoi,A,F_total,B_total,P_total\n"<<std::setprecision(12);for(int i=0;i<N;++i)if(hc[i]){int x=i%NX,y=i/NX;o<<cfg.steps<<','<<cfg.steps*.1<<','<<x<<','<<y<<','<<.1*(x+.5*y)<<','<<.08660254037844386*y;for(int f=0;f<NF;++f)o<<','<<hq[f*N+i];double ft=0,bt=0,pt=0;for(int m=0;m<NO;++m){ft+=hF[m*N+i];bt+=hB[m*N+i];pt+=hP[m*N+i];}o<<','<<ft<<','<<bt<<','<<pt<<'\n';}}
 void run(){fs::create_directories(cfg.output);std::ofstream out(fs::path(cfg.output)/"metrics.csv");out<<"step,time_s,sites,area_um2,perimeter_um,cx_um,cy_um,extensions,retractions,expected_extensions,expected_retractions\n";auto write=[&](int step,CopyDiag d){Totals t=totals(after);out<<step<<','<<step*.1<<','<<t.sites<<','<<t.sites*.01<<','<<t.perimeter/1000<<','<<t.x/t.sites/1000<<','<<t.y/t.sites/1000<<','<<d.ext_accepted<<','<<d.ret_accepted<<','<<d.ext_expected<<','<<d.ret_expected<<'\n';};write(0,{});for(int s=1;s<=cfg.steps;++s){fields((s-1)*.1);write(s,cfg.frozen?CopyDiag{}:cpm(s));}CUDA_OK(cudaDeviceSynchronize());if(cfg.snapshot)snapshot();}
};

int main(int argc,char**argv){try{Config c;for(int i=1;i<argc;++i){std::string a=argv[i];auto val=[&](){if(++i>=argc)throw std::runtime_error("missing CLI value");return std::string(argv[i]);};if(a=="--steps")c.steps=std::stoi(val());else if(a=="--seed")c.seed=std::stoull(val());else if(a=="--output")c.output=val();else if(a=="--lambda-cpm")c.lambda_cpm=std::stod(val());else if(a=="--no-pulse")c.pulse=false;else if(a=="--snapshot")c.snapshot=true;else if(a=="--frozen")c.frozen=true;else throw std::runtime_error("unknown argument "+a);}GPUModel m(c);m.run();}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 1;}return 0;}
