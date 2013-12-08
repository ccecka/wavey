#ifndef HELMHOLTZUTIL_H
#define HELMHOLTZUTIL_H

#include "Vec3.hpp"

// Tightness of the error bound [0,1]
//const double MLFMM_ALPHA(1.0);
//const double MLFMM_ALPHA(0.95);
const double MLFMM_ALPHA(0.9);
//const double MLFMM_ALPHA(0.85);
//const double MLFMM_ALPHA(0.8);
//const double MLFMM_ALPHA(0.7071);
//const double MLFMM_ALPHA(0.57735);


inline complex gegenbauer_series( int L, double kappa,
				  const Vec3& r, const Vec3& r0 )
{
  double kr0 = kappa * r0.mag();
  double kr  = kappa * r.mag();
  double rr0 = r0.dot(r) / (r0.mag() * r.mag());

  complex G = 0;
  int one = 1;
  for( int n = 0; n <= L; ++n ) {
    G += one * (2*n+1) * bessel_h(n,kr0) * bessel_j(n,kr) * legendre_P(n,rr0);
    one *= -1;
  }

  return CI * kappa * G;
}

// A function to evaluate the transfer function T_{L,r_0}(s)
inline complex transfer_function( int L, double kappa,
				  const Vec3& r0, const Vec3& s )
{
  assert( abs(s.mag()-1) < 1e-15 );
  double sdotr0 = r0.dot(s) / r0.mag();
  double knormr0 = kappa * r0.mag();

  complex t = 0;
  complex i_k = 1;
  for( int n = 0; n <= L; ++n ) {
    t += i_k * (2*n+1) * bessel_h(n,knormr0) * legendre_P(n,sdotr0);
    i_k *= CI;
  }

  return (CI*kappa)/(4*PI) * t;
}


// A class to quickly evaluate the transfer function T_{ell,r_0}(s)
// for many possible directions s
struct Transfer_Function_Eval
{
  int L;
  Vec3 r0hat;
  vector<complex> H;   // Coefficients of the series
  vector<double> P;

Transfer_Function_Eval( double kappa, int L_, Vec3 r0 )
: L(L_), r0hat(r0/r0.mag()), H(L+1), P(L+1)
  {
    double knormr0 = kappa * r0.mag();

    complex i_k = (CI*kappa)/(4*PI);
    for( int n = 0; n <= L; ++n ) {
      H[n] = i_k * (2*n+1) * bessel_h(n,knormr0);
      i_k *= CI;
    }
  }

  inline complex operator()(const double sdotr0)
  {
    gsl_sf_legendre_Pl_array(L, sdotr0, &P[0]);

    complex sum = H[0] * P[0];
    for( int n = 1; n <= L; ++n )
      sum += H[n] * P[n];

    return sum;
  }

  inline complex operator()(const Vec3& s)
  {
    assert( abs(s.mag()-1) < 1e-15 );
    return operator()( r0hat.dot(s) );
  }
};


inline int ell_ebf(double kappa, double norm_r, double eps)
{
  double knormr = kappa * norm_r;
  return ceil(knormr + 1.8*pow(-log10(eps), 2.0/3) * pow(knormr, 1.0/3));
}


// Get the Gegenbauer truncation, L, via a number of possible methods
inline static int get_Truncature( double kappa,
				  double norm_r, double norm_r0,
				  double eps, int method = 4 )
{
  int ell = 0;

  switch( method ) {
    case 1: // Old Log Method
      ell = (int) (kappa*norm_r - log10(eps)*log(PI + kappa*norm_r));
      break;

    case 2: // EBF Formula (Chew)
      ell = ell_ebf(kappa, norm_r, eps);
      break;

    case 3: {// Direct method (Collino) with worst case cos = -1
      // Compute the EBF_L as an upperbound
      double knorm_r  = kappa*norm_r;
      double knorm_r0 = kappa*norm_r0;

      // A slightly more accurate ebf_ell to use as a starting point
      int ebf_ell = get_Truncature( kappa, norm_r, norm_r0, 1e-2 * eps, 2 );
      //int ebf_ell = (int) (knorm_r + 1.8*pow(2-log10(eps), 2.0/3.0)
      //	 		              *pow(knorm_r, 1.0/3.0));

      ell = ebf_ell;

      double eM = knorm_r*knorm_r0 / abs(norm_r0-norm_r);
      double eP = knorm_r*knorm_r0 / abs(norm_r0+norm_r);

      complex hl = bessel_h(ell, knorm_r0), hlp1 = bessel_h(ell+1, knorm_r0);
      double  jl = bessel_j(ell, knorm_r ), jlp1 = bessel_j(ell+1, knorm_r );

      // Decrease ell until we hit the cutoff
      double error = max( eM*abs(hlp1*jl-hl*jlp1), eP*abs(hlp1*jl+hl*jlp1) );
      //cout << "Trunc Error: " << error << "    ell: " << ell << endl;
      while( error < eps && ell > 1 ) {
        --ell;
        hlp1 = hl;  hl = bessel_h( ell, knorm_r0 );
        jlp1 = jl;  jl = bessel_j( ell, knorm_r  );
        error = max( eM*abs(hlp1*jl-hl*jlp1), eP*abs(hlp1*jl+hl*jlp1) );
        //cerr << "Trunc Error: " << error << "    ell: " << ell << endl;
      }
      // Went one step over
      ++ell;

      // Check this against the true Gegenbauer error
      // If it agrees, we're good
      // If it doesn't agree, then |r| is too small and we should
      //          do an ultradirect check (using the true Geg error)

      break; }
    case 4: {// Ultradirect method - check again gegnbauer series itself
      int ell_ebf = get_Truncature(kappa, norm_r, norm_r0, eps, 2);
      ell = ell_ebf;

      double max_error = 1e100;
      while (max_error > eps) {
        ++ell;
        max_error = 0;

        int N = 20;
        for (int n = 0; n <= N; ++n) {
          Vec3 r  = Vec3(0,0,1) * norm_r;
          Vec3 r0 = Vec3(sin((PI*n)/N),0,cos((PI*n)/N)) * norm_r0;

          double R = (r+r0).mag();
          complex Ie = exp(CI * kappa* R) / R;

          complex GS = gegenbauer_series(ell, kappa, r, r0);

          //cout << Ie << "\t" << GS << "\t" << abs(Ie-GS) << endl;

          max_error = max(max_error, abs(Ie - GS));
        }

        //cout << ell << ": " << max_error << endl;
      }

      //cout << "Truncation - EBF: " << ell_ebf << "  UD: " << ell << endl;
      //ell = min(ell,ell_ebf);

      break; }
  }
  return ell;
}



/* Computes a low-pass |sin(phi)|
 *
 * nF: number of desired frequencies of abs(sin)
 * nR: number of desired real space points of abs(sin), nR >= 2*nF+1
 */
inline vector<complex> FabsSin( int nF, int nR )
{
  assert( nR >= 2*nF+1 );
  vector<complex> abssin(nR,0);

  // Compute all the Fourier coefficients
  abssin[0] = 2.0/PI;
  for( int k = 2; k <= nF; k += 2 )
    abssin[nR - k] = abssin[k] = 2.0/(PI*(1 - k*k));

  fft( &abssin[0], nR );
  return abssin;
}


// Get the nth coefficient of |sin|
inline double sincoef( int n )
{
  if( ISODD(n) )  return 0;
  else            return 2.0/(PI*(1 - n*n));
}


// Takes a [-L,L] and convolves with |sin| to get [-M,M]
inline vector<complex> slowConv( vector<complex> a, int L, int M )
{
  ifft( &a[0], 2*L+1 );

  vector<complex> F(2*M+2,0);

  for( int m = 0; m <= M; ++m ) {
    for( int k = 0; k <= L; ++k ) {
      F[m] += a[k] * sincoef(m-k);
    }
    for( int k = -L; k < 0; ++k ) {
      F[m] += a[2*L+1+k] * sincoef(m-k);
    }
  }

  for( int m = -M; m < 0; ++m ) {
    for( int k = 0; k <= L; ++k ) {
      F[2*M+2+m] += a[k] * sincoef(m-k);
    }
    for( int k = -L; k < 0; ++k ) {
      F[2*M+2+m] += a[2*L+1+k] * sincoef(m-k);
    }
  }

  fft( &F[0], 2*M+2 );
  return F;
}


// Computes the low-pass modified transfer function matrix N_rows x N_phi
// L:  Gegenbauer truncation
// kappa: Wavenumber
// N_phi: the number of quadrature points in phi
// N_rows: the number of rows of the quadrature N_theta = 2*N_rows - 2
//
// Returns: T[n][m] nth row, mth phi value of the low-pass modified trans fn
inline vector<vector<complex> > mod_transfer_fn( const int L,
						 const double kappa,
						 const int N_rows,
						 const int N_phi,
						 const Vec3& r0 )
{
  assert( !ISODD(N_phi) );

  const Vec3 r0hat = r0 / r0.mag();

  // Construct a Transfer_Function_Eval T_{L,r_0}
  Transfer_Function_Eval Tfn( kappa, L, r0 );

  const int N_theta = 2*N_rows - 2;           // # points used for ModTransFn
  const int N_t = 2*L + 1;                    // # points needed for TransFn

  // Get the low-pass |sin(theta)|
  const int N_sinF = (N_theta/2-1) + L;       // # freqs needed in |sin|
  //int N_conv = 2*(L+N_sinF) + 1;
  const int N_conv = N_theta + 2*L - 1;       // size of convolution with |sin|
  vector<complex> abssin = FabsSin(N_sinF, N_conv);

  const double scale = 0.5 / (N_t*N_conv);    // 1/2 and interp/anterp coeffs

  // The N_rows x N_phi transfer function
  vector< vector<complex> > T(N_rows, vector<complex>(N_phi,0));

  // Temp storage for t_phi(theta)
  vector<complex> t(N_conv);

  // FFTs for t(theta)    (Optimized)
  fftw_complex* a = reinterpret_cast<fftw_complex*>(&t[0]);
  fftw_plan tFFT1 = fftw_plan_dft_1d(N_t,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan tIFFT1 = fftw_plan_dft_1d(N_conv,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_plan tFFT2 = fftw_plan_dft_1d(N_conv,a,a,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan tIFFT2 = fftw_plan_dft_1d(N_theta,a,a,FFTW_BACKWARD,FFTW_ESTIMATE);


  // For each phi [0,PI)
  for( int m = 0; m < N_phi/2; ++m ) {
    const double phi = m * (2*PI)/N_phi;
    const double cosP = cos( phi );
    const double sinP = sin( phi );
    const double xydotr0 = cosP*r0hat.x + sinP*r0hat.y;

    // Sample the transfer function at 2L+1 points in theta [0,2PI)
    for( int n = 0; n < N_t; ++n ) {
      const double theta = n * (2*PI)/N_t;
      const double sinT = sin( theta );
      const double cosT = cos( theta );

      //Vec3 s( sintheta*cosphi, sintheta*sinphi, costheta );
      //t[n] = scale * transfer_function(L, kappa, r0, s);
      const double sdotr0 = sinT*xydotr0 + cosT*r0hat.z;
      t[n] = scale * Tfn( sdotr0 );    // Optimized
    }

    // Interpolate the transfer function to N_conv
    fftw_execute( tFFT1 );
    f_cut( &t[0], N_t, N_conv );
    fftw_execute( tIFFT1 );

    // Multiply by low-pass |sin| to get the modified transfer function
    for( int k = 0; k < N_conv; ++k )
      t[k] *= abssin[k];

    // Anterpolate to N_theta
    fftw_execute( tFFT2 );
    f_cut( &t[0], N_conv, N_theta );
    t[N_theta/2] = 0;                    // Set the oddball freq to zero
    fftw_execute( tIFFT2 );

    // Explicit convolution to check the answer
    //vector<complex> exact = slowConv( t, L, N_theta/2-1 );

    // Unwrap into the N_rows x N_phi matrix
    T[0][m] = T[0][N_phi/2+m] = t[0];       // North Poles
    for( int n = 1; n < N_rows; ++n ) {
      T[n][m]         = t[n];
      T[n][N_phi/2+m] = t[N_theta-n];
    }
  }

  fftw_destroy_plan( tFFT1 );
  fftw_destroy_plan( tIFFT1 );
  fftw_destroy_plan( tFFT2 );
  fftw_destroy_plan( tIFFT2 );

  return T;
}

/*
inline vector< vector<complex> > translation_fn( double kappa,
						 int N_rows, int N_phi,
						 const Vec3& r )
{
  // The N_rows x N_phi translation function
  vector< vector<complex> > E(N_rows, vector<complex>(N_phi,0));

  // For each theta [0,PI]
  for( int n = 0; n < N_rows; ++n ) {
    double theta = n * PI/(N_rows-1);
    double sintheta = sin( theta );
    double costheta = cos( theta );

    // For each phi [0,2PI)
    for( int m = 0; m < N_phi; ++m ) {
      double phi = m * (2*PI)/N_phi;
      double cosphi = cos( phi );
      double sinphi = sin( phi );

      Vec3 s( sintheta*cosphi, sintheta*sinphi, costheta );
      E[n][m] = exp( CI * kappa * r.dot(s) );
    }
  }

  return E;
}

// Applies the symmetry
// f( theta_n, phi_m ) = f( theta_{N_theta-n}, phi_{N_phi/2+m} )
// to get the whole spherical data
inline vector< vector<complex> > unwrap( vector< vector<complex> >& f )
{
  int N_rows = f.size();
  int N_theta = 2*N_rows - 2;
  int N_phi = f[0].size();    // Needs to equal all f[n].size()

  vector< vector<complex> > f_full( N_theta, vector<complex>(N_phi,0) );

  // For each phi
  for( int m = 0; m < N_phi; ++m ) {
    vector<complex> temp(N_theta, 0);

    // Unwrap
    f_full[0][m] = f[0][m];                        // North Pole
    for( int n = 1; n < N_rows; ++n ) {
      f_full[n][m] = f[n][m];
      if( m < N_phi/2 )
	f_full[N_theta-n][m] = f[n][m+N_phi/2];
      else
	f_full[N_theta-n][m] = f[n][m-N_phi/2];
    }
  }

  return f_full;
}

// Takes a N_rows x N_phi sampled function and returns the full 2D IFFT
inline vector< vector<complex> > ifft2_sphere( vector< vector<complex> >& F )
{
  int N_rows = F.size();
  int N_theta = 2*N_rows - 2;
  int N_phi = F[0].size();    // Needs to equal all f[n].size()

  // Compute the 2D FFT
  vector< vector<complex> > f(N_theta, vector<complex>(N_phi,0));

  // For each phi
  for( int m = 0; m < N_phi; ++m ) {
    vector<complex> temp(N_theta, 0);

    // Unwrap
    temp[0] = F[0][m];                        // North Pole
    for( int n = 1; n < N_rows; ++n ) {
      temp[n] = F[n][m];
      if( m < N_phi/2 )
	temp[N_theta-n] = F[n][m+N_phi/2];
      else
	temp[N_theta-n] = F[n][m-N_phi/2];
    }

    ifft( &temp[0], N_theta );

    // Unwrap into the N_theta x N_phi matrix
    for( int n = 0; n < N_theta; ++n ) {
      f[n][m] = temp[n];
    }
  }

  // Take the FFT along each theta to get 2D FFT
  // For each theta
  for( int n = 0; n < N_theta; ++n ) {
    ifft( &f[n][0], N_phi );
  }

  return f;
}
*/

#endif
