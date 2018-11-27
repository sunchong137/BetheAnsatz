/* This program is released under a "I hope no one has to deal with this
 * ever again"-license. You may do with it as you wish, including
 * distributing it in source- or binary form, as long as you do not blame
 * me if anything goes wrong.
 * 
 * If you find the program useful, I would appreciate a citation to
 * 
 *    Knizia, Chan - Density Matrix Embedding: A Simple Alternative
 *       to Dynamical Mean-Field Theory
 *    http://dx.doi.org/10.1103/PhysRevLett.109.186404
 * 
 * in published work using the program. That is the paper for which it
 * was written. In any case you should cite
 * 
 *    Shiba - Magnetic Susceptibility at Zero Temperature for the
 *       One-Dimensional Hubbard Model
 *    http://dx.doi.org/10.1103/PhysRevB.6.930
 * 
 * which describes the algorithm which was employed in this program.
 * -- Gerald Knizia, 2014-06-03
 */

// calculates 1D Hubbard model E(U,f) where f is the filling; based on
// [1] Shiba - Magnetic susceptibility at zero temperture for the
//     one-dimensional Hubbard model [Phys Rev B 6 930 (1972)]
//
// License: Use as you see fit, but don't blame me if something goes
//          wrong!
//
//                                   -- Gerald Knizia, 2011
//
// Compile:
//    c++ -O3 -DNDEBUG bethe_ansatz.cpp -o bethe_ansatz
// Run:
//    Meant to be run by bethe_ansatz.py.
//    This program just takes two arguments (U and Q) and calculates
//    the energy and the particle number resulting from this. Other properties
//    can be obtained as derivatives from those.
//
// See README.txt

#include <iostream>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <cstdlib> // atof/atoi
#include <vector>
#include <boost/format.hpp>

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

using std::cout;
using std::endl;
using boost::format;
using boost::str;

typedef double
   scalar;
typedef unsigned int
   uint;

#define USE_ALTERNATING_SUM_ACCEL
#if not defined(USE_ALTERNATING_SUM_ACCEL)
   // Very slow. This branch just for testing.

   // [0] eq. 2.7.
   static scalar R(scalar x)
   {
      scalar
         r = 0,
         sign = 1.0;  // start at + for n == 1.
      for (scalar n = 1; n < 2000000 ; ++ n) {
         scalar d = 2*n/(x*x + 4*n*n);
         r += sign * d;
         sign = -sign;
         if ((((long)n) & 0xfff) == 1)
            cout << (format("   n = %i   r = %.8f   d = %f") % n % r % d) << endl;
         if (std::abs(d/r) < 1e-16)
            break;
      }
      return r/M_PI;
   }
#else
   // Convergence acceleration of alternating series
   // Henri Cohen, Fernando Rodriguez Villegas, and Don Zagier
   // Experiment. Math. Volume 9, Issue 1 (2000), 3-12.
   typedef scalar (*sum_fn_t)(uint k, double x);

   inline scalar r_k(uint n, scalar x){
      n += 1; // summation starts at n=0 here with positive sign, that's why we shift from (2.7) by one.
      return 2*n/(x*x + 4*n*n);
   }

   static scalar alternating_sum(uint n, sum_fn_t an, scalar param)
   {
      scalar
         b = -1, c, d, s = 0;
      d = std::pow((3. + std::sqrt(8.)), n);
      d = (d+1/d)/2;
      c = -d;
      for (scalar k = 0; k < n; ++ k) {
         c = b - c;
         s = s + c * an(k, param);
         b *= (k + n)*(k - n)/((k+.5)*(k+1));
      }
      return s/d;
   }

   static scalar R(scalar x)
   {
      scalar e2 = alternating_sum(22, r_k, x);
   #ifdef _DEBUG
      scalar e1 = alternating_sum(20, r_k, x),
      if (abs(e1-e2) > 1e-15) throw std::runtime_error("N too small in R(x)!");
   #endif // _DEBUG
      return e2/M_PI;
   }
#endif // USE_ALTERNATING_SUM_ACCEL


uint
   // we use an equidistant grid. That is fine in this case, because the
   // function to integrate falls of quickly at the boundaries and is also
   // otherwise very well-behaved. Higher order integration rules are
   // unnecessary, and the 200pt grid is also plenty.
   nGrid = 200;

// grid positions x(k)
inline scalar x(uint k, scalar Q) {
   return Q*(2*(signed)k + 1 - (signed)nGrid)/static_cast<scalar>(nGrid);
}

// integration weight dx(k)
inline scalar dx(uint k, scalar Q) {
   return (2.*Q)/nGrid;
}

// inline scalar dx(uint k, scalar Q) {
//    scalar f0 = (1./3.)*(2.*Q)/nGrid;
//    assert(nGrid % 2 == 1 );
//    if ( k == 0 || k == nGrid -1 ) return f0;
//    if ( k % 2 == 1 ) return 4.*f0;
//    else return 2.*f0;
// };
// ^- hm... not very convincing.


static scalar eval_rho_rhs_for_Q(scalar &EperSite, scalar &NbyNa, scalar *rho_out, scalar *rho_in, scalar Q, scalar u, std::vector<scalar> &tmp)
{
   // evaluates one iteration of [0] eq. 2.8 for a given rho, on the
   // interval [-Q,Q] with nGrid grid points.
   scalar
      invu = 1./u,
      max_dev = 0; // maximal change in rho.
   EperSite = 0; // E/Na
   NbyNa = 0;
   tmp.resize(nGrid*2);
   scalar
      *cosxi = &tmp[0],
      *sinxi = &tmp[nGrid];
   for (uint i = 0; i < nGrid; ++ i) {
      scalar xi = x(i,Q);
      cosxi[i] = std::cos(xi);
      sinxi[i] = std::sin(xi);
   }

   for (uint i = 0; i < nGrid; ++ i) {
      scalar fint = 0;
      for (uint k = 0; k < nGrid; ++ k) {
         fint += dx(k,Q) * rho_in[k] * R(4.*invu*(sinxi[i] - sinxi[k]));
      }
      fint *= 8.*M_PI*invu;
      max_dev = std::max(std::abs(rho_out[i] - rho_in[i]), max_dev);
      rho_out[i] = (1.0/(2.*M_PI)) * (1 + cosxi[i] * fint);
      EperSite += -2 * dx(i,Q) * cosxi[i] * rho_out[i];
      NbyNa += dx(i,Q) * rho_out[i];
   }
   return max_dev;
}


struct result_entry_t {
   scalar u;
   scalar filling;
   scalar energy_per_site;
   scalar Q;
};

static void solve_for_rho(result_entry_t *r, scalar *rho, scalar Q, scalar u, int print = 1)
{
   double dabl = 0.;
   for (uint i = 0; i < nGrid; ++ i)
      dabl += dx(i,Q);
   if (print >= 2) {
      cout << format("sum[dx] = %.8f  2*Q=%.8f") % dabl % (2*Q) << endl;
      cout << format("x[0] = %.8f") % x(0,Q) << endl;
      cout << format("x[nGrid-1] = %.8f") % x(nGrid-1,Q) << endl;
   }
   r->Q = Q;
   r->u = u;

   bool
      converged = false;
   std::vector<scalar>
      dummy(nGrid, 0.0),
      tmp;
   scalar
      *pin = &dummy[0], *pout = rho, last_e = 1e99;
   for (uint it = 0; it < 1000; ++ it) {
      double max_dev_rho = eval_rho_rhs_for_Q(r->energy_per_site, r->filling, pout, pin, Q, u, tmp);
      double delta_e = last_e - r->energy_per_site;
      if (print >= 1) {
         cout << format("it = %-4i  E = %10.5f  N/Na = %10.5f  max_dev_rho = %7.2e  delta_e = %7.2e  rho = [")
            % it % r->energy_per_site % r->filling % max_dev_rho % delta_e;
         uint nPrint = 7;
         for ( uint i = 0; i < nPrint; ++ i )
            cout << format(" %10.5f") % pout[(i*nGrid)/(nPrint-1)];
         cout << "]" << endl;
      }
      std::swap(pin, pout);
      if (it > 2 && max_dev_rho < 1e-15 && std::abs(delta_e) < 1e-15) {
         converged = true;
         break;
      }
      last_e = r->energy_per_site;
   }
   if (!converged)
      throw std::runtime_error(str(format("sorry, bethe_ansatz.cpp failed to converge! Q=%.8f u=%8f.") % Q % u));
   if (rho != pin)
      std::memcpy(rho, pin, sizeof(rho[0])*nGrid);
}

int main(int argc, char **argv)
{
   // note: t is fixed at 1.0.

   if (argc != 3)
      throw std::runtime_error("usage: bethe <U> <Q>    with U, Q floating point numbers. Q \\in [0,1]. U < 0.05 should not be used.");
   scalar u = atof(argv[1]);
   scalar q = atof(argv[2]);

   // for small U we need more grid points---need to accurately
   // sample R(4*sin(x)/u) for x = -pi..pi.
   if (u < 1.) {
      nGrid *= 1./u;
      nGrid += (nGrid % 2);
      std::cout << format("WARNING: Increased number of grid points to %i due to small U. May get slow!") % nGrid << std::endl;
   }

   std::vector<scalar>
      rho(nGrid, 0.0);
   result_entry_t r;
   solve_for_rho(&r, &rho[0], q * M_PI, u, 0);
   std::cout << format("E= %.15f  n= %.15f") % r.energy_per_site % r.filling << std::endl;
   return 0;
}


// Note on the source formatting:
//    cgk: sometimes I try new coding styles on small projects, to see
//    if I like them. Turns out that I do not like this one at all.


//
// (note: see also Cole III - ``Excitation spectrum of the one-dimensional
//  Hubbard model'' for some more information.
//  On excited states: can't we just take higher eigenvectors of the
//  rho self-consistency equation? Since rho has a maximum at k=0, it seems
//  that we should just solve for higher distributions. But what about the
//  E = int cos(k) rho(k)  instead of E = int rho[k]?.. we need orthogonal
//  wave functions, not density distributions. However, Cole's approach
//  does not apply to the half filled case (Woynarovich, J Phys C 15 85 (1982)))


// kate: indent-width 3
