#ifndef TRANSFER_FUNCTION_H
#define TRANSFER_FUNCTION_H

#include "General.hpp"
#include "NFunction.hpp"
#include "HelmholtzUtil.hpp"

  /* helper for code location */
  #define LOC std::cerr << __FILE__ << ":" << __LINE__ << "  "
  /* macro for general debug print statements. */
  #define COUT_TXT(text) LOC; std::cerr << text << std::endl
  /* macro that dumps a variable name and its value */
  #define COUT_VAR(var) LOC; std::cerr << (#var) << " = " << var << std::endl

class Transfer_Function : public NFunction
{
 public:

  Transfer_Function( Quadrature* quad, const Vec3& r0)
      : NFunction(quad)
  {
    //cerr << "r0 = " << r0 << endl;

    double kappa = quad->kappa;
    int L = quad->L;
    int N_rows = quad->numRows();
    int N_phi = 2*L + 2;

    COUT_VAR(kappa);
    COUT_VAR(L);
    COUT_VAR(N_rows);
    COUT_VAR(N_phi);
    COUT_VAR(r0);

    // Get the low-pass modified transfer matrix N_rows x N_phi
    vector< vector<complex> > T = mod_transfer_fn(L, kappa, N_rows, N_phi, r0);

    // maxRowSize can be larger than Nphi -> shouldn't interpolate in place
    vector<complex> t( max(quad->maxRowSize(),N_phi) );

    // Anterpolate the rows for an optimized quadrature
    for( int n = 0; n < N_rows; ++n ) {
      QuadratureRow& row = quad->getRow(n);
      int N_phi_n = row.size();

      // Anterpolate to N_phi_k
      for( int m = 0; m < N_phi; ++m ) t[m] = T[n][m];
      fftinterp( &t[0], N_phi, N_phi_n );

      // Copy into the NFunction
      for( int m = 0; m < N_phi_n; ++m ) {
        int index = row.getPoint(m).index;
        C[index] = t[m];
      }
    }

    if (r0.x > 0)
      std::cerr << C << std::endl;
  }

};

#endif
