#include <Rcpp.h>
using namespace Rcpp;

struct ReactionParameters {
  double n, a1, a2, alpha, beta, ctot, ractot, rhotot;
  double dcdc, drac, drho, f, ip, ip1, ir, k21, kpi3k, kpi5k;
  double kpten, p3b, rb, rhob, delta_p1;
};

ReactionParameters unpack(const NumericVector& p) {
  return {p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9],
          p[10], p[11], p[12], p[13], p[14], p[15], p[16], p[17], p[18],
          p[19], p[20], p[21], p[22]};
}

double rates(const NumericMatrix& u, NumericMatrix& du, const NumericVector& x,
             const double time, const ReactionParameters& p) {
  const int nr = u.nrow();
  double sum_j5 = 0.0;
  for (int i = 0; i < nr; ++i) {
    const double ic = time <= 10.0 ? 2.6 + 0.05 * x[i] : 2.95;
    const double cdc_i = u(i, 0), cdc_a = u(i, 1);
    const double rac_i = u(i, 2), rac_a = u(i, 3);
    const double rho_i = u(i, 4), rho_a = u(i, 5);
    const double p1 = u(i, 6), p2 = u(i, 7), p3 = u(i, 8);
    const double lipid_feedback = (1.0 - p.f) + p.f * p3 / p.p3b;
    const double j0 = ic / (1.0 + std::pow(rho_a / p.a1, p.n)) *
      (cdc_i / p.ctot) * lipid_feedback - p.dcdc * cdc_a;
    const double j1 = (p.ir + p.alpha * cdc_a) * (rac_i / p.ractot) *
      lipid_feedback - p.drac * rac_a;
    const double j2 = (p.ip + p.beta * rac_a) /
      (1.0 + std::pow(cdc_a / p.a2, p.n)) * (rho_i / p.rhotot) -
      p.drho * rho_a;
    const double j3 = -p.k21 * p2 + p.kpi5k / 2.0 *
      (1.0 + rac_a / p.rb) * p1;
    const double j4 = p.kpi3k / 2.0 * (1.0 + rac_a / p.rb) * p2 -
      p.kpten / 2.0 * (1.0 + rho_a / p.rhob) * p3;
    const double j5 = p.ip1 - p.delta_p1 * p1;
    du(i, 0) = -j0; du(i, 1) = j0;
    du(i, 2) = -j1; du(i, 3) = j1;
    du(i, 4) = -j2; du(i, 5) = j2;
    du(i, 6) = -j3 + j5;
    du(i, 7) = j3 - j4;
    du(i, 8) = j4;
    sum_j5 += j5;
  }
  return sum_j5 / nr;
}

// [[Rcpp::export]]
List reaction_step_cpp(const NumericMatrix& state, const NumericVector& x,
                       const double time, const double dt,
                       const NumericVector& parameters) {
  const int nr = state.nrow(), nc = state.ncol();
  ReactionParameters p = unpack(parameters);
  NumericMatrix k1(nr, nc), k2(nr, nc), k3(nr, nc), k4(nr, nc), tmp(nr, nc);
  const double m1 = rates(state, k1, x, time, p);
  for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i)
    tmp(i,j) = state(i,j) + 0.5 * dt * k1(i,j);
  const double m2 = rates(tmp, k2, x, time + 0.5 * dt, p);
  for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i)
    tmp(i,j) = state(i,j) + 0.5 * dt * k2(i,j);
  const double m3 = rates(tmp, k3, x, time + 0.5 * dt, p);
  for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i)
    tmp(i,j) = state(i,j) + dt * k3(i,j);
  const double m4 = rates(tmp, k4, x, time + dt, p);
  NumericMatrix result(nr, nc);
  for (int j = 0; j < nc; ++j) for (int i = 0; i < nr; ++i)
    result(i,j) = state(i,j) + dt / 6.0 *
      (k1(i,j) + 2.0*k2(i,j) + 2.0*k3(i,j) + k4(i,j));
  const double pi_delta = -dt / 6.0 * (m1 + 2.0*m2 + 2.0*m3 + m4);
  return List::create(_["state"] = result, _["pi_delta"] = pi_delta);
}
