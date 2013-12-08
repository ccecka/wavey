#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Quadrature.hpp"
#include "NFunction.hpp"
#include "General.hpp"

#include <stdio.h>
#include <string.h>
#include <memory.h>

#include <gsl/gsl_cblas.h>

// A reterpolator between two quadratures using FFT interpolation
// Optimized for use with single-point poles
class FFTInterp_Poles
{
  Quadrature* q1;           // Pointer to quadrature 1 (from quad)
  Quadrature* q2;           // Pointer to quadrature 2 (to quad)

  // q1 Data
  int N_rows1, N_theta1, N_phi1;

  // q2 Data
  int N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  int N_rows, N_phi;

  // Intermediate Storage  N_rows x N_phi
  vector<complex> T;
  // Intermediate Storage  N_theta
  vector<complex> t;

  // First stage phi (row) FFTs
  vector<fftw_plan> rFFT1;
  fftw_plan rIFFT1;

  // Second stage theta (column) FFTs
  fftw_plan tFFT;
  fftw_plan tIFFT;

  // Third stage phi (row) FFTs
  fftw_plan rFFT2;
  vector<fftw_plan> rIFFT2;

 public:

  // Constructor
  FFTInterp_Poles( Quadrature* q1_, Quadrature* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_theta1(2*N_rows1-2), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_theta2(2*N_rows2-2), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ), N_phi( max(N_phi1,N_phi2) ),
        T( N_phi * N_rows ), t( max(N_theta1,N_theta2) )
  {
    fftw_complex* T0 = reinterpret_cast<fftw_complex*>(&T[0]);
    fftw_complex* t0 = reinterpret_cast<fftw_complex*>(&t[0]);

    // Define an phi1 FFT for each row
    rFFT1 = vector<fftw_plan>( N_rows1 );
    for( int n = 1; n < N_rows1-1; ++n ) {
      int N_phi_n = q1->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rFFT1[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_FORWARD,FFTW_MEASURE);
    }

    // Define a phi block IFFT for each row
    rIFFT1 = fftw_plan_many_dft( 1, &N_phi, N_rows1-2,
                                 T0 + N_phi, &N_phi, 1, N_phi,
                                 T0 + N_phi, &N_phi, 1, N_phi,
                                 FFTW_BACKWARD, FFTW_MEASURE );

    // Define a theta FFT and IFFT (could also do these as blocks...)
    tFFT = fftw_plan_dft_1d(N_theta1, t0, t0, FFTW_FORWARD, FFTW_MEASURE);
    tIFFT = fftw_plan_dft_1d(N_theta2, t0, t0, FFTW_BACKWARD, FFTW_MEASURE);

    // Define a phi2 block FFT for each row (except poles)
    rFFT2 = fftw_plan_many_dft( 1, &N_phi, N_rows2-2,
                                T0 + N_phi, &N_phi, 1, N_phi,
                                T0 + N_phi, &N_phi, 1, N_phi,
                                FFTW_FORWARD, FFTW_MEASURE );

    // Define an phi2 IFFT for each row (except poles)
    rIFFT2 = vector<fftw_plan>( N_rows2 );
    for( int n = 1; n < N_rows2-1; ++n ) {
      int N_phi_n = q2->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rIFFT2[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_BACKWARD,FFTW_MEASURE);
    }
  }

  // Destructor
  ~FFTInterp_Poles() {
    for( int n = 1; n < N_rows1-1; ++n )
      fftw_destroy_plan( rFFT1[n] );
    fftw_destroy_plan( rIFFT1 );
    fftw_destroy_plan( tFFT );
    fftw_destroy_plan( tIFFT );
    fftw_destroy_plan( rFFT2 );
    for( int n = 1; n < N_rows2-1; ++n )
      fftw_destroy_plan( rIFFT2[n] );
  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( A.quad == q1 && B.quad == q2 );

    double scale = 1.0/(N_theta1*N_phi);

    // Copy from A and interpolate each row
    for( int n = 1; n < N_rows1-1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      // Copy into the kth row of T (scale for the upcoming interpolations)
      for( int m = 0; m < N_phi_n; ++m ) {
        Tn[m] = (scale/N_phi_n) * A.C[ rown.getPoint(m).index ];
      }

      fftw_execute( rFFT1[n] );                    // FFT_phi
      f_cut( Tn, N_phi_n, N_phi );                 // Smooth phi
    }

    fftw_execute( rIFFT1 );                        // IFFT_phi all rows

    // Interpolate from N_theta1 to N_theta2
    complex poleN = scale * A.C[ q1->getRow(0).getPoint(0).index ];
    complex poleS = scale * A.C[ q1->getRow(N_rows1-1).getPoint(0).index ];
    complex sumN = 0, sumS = 0;
    for( int m = 0; m < N_phi/2; ++m ) {

      // Unwrap into t
      t[0] = poleN;
      for( int n = 1; n < N_rows1-1; ++n ) {
        t[n]          = T[m         + n*N_phi];
        t[N_theta1-n] = T[m+N_phi/2 + n*N_phi];
      }
      t[N_rows1-1] = poleS;

      // Interpolate from N_theta1 to N_theta2
      fftw_execute( tFFT );                        // FFT_theta
      f_cut( &t[0], N_theta1, N_theta2 );          // Smooth theta
      //t[N_theta2/2] = 0;                           // Zero oddball?
      fftw_execute( tIFFT );                       // IFFT_theta

      // Unwrap back into T
      sumN += t[0];
      for( int n = 1; n < N_rows2-1; ++n ) {
        T[m         + n*N_phi] = t[n];
        T[m+N_phi/2 + n*N_phi] = t[N_theta2-n];
      }
      sumS += t[N_rows2-1];
    }

    fftw_execute( rFFT2 );                         // FFT_phi all rows

    // Interpolate each row and copy into B
    B.C[ q2->getRow(0).getPoint(0).index ] = 2.0 * sumN;
    for( int n = 1; n < N_rows2-1; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      f_cut( Tn, N_phi, N_phi_n );                 // Smooth phi
      //if( N_phi_n != 1 ) Tn[N_phi_n/2] = 0;      // Zero oddball?
      fftw_execute( rIFFT2[n] );                   // IFFT phi

      // Copy into B
      for( int m = 0; m < N_phi_n; ++m ) {
        B.C[ rown.getPoint(m).index ] = Tn[m];
      }
    }
    B.C[ q2->getRow(N_rows2-1).getPoint(0).index ] = 2.0 * sumS;
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_Poles(const FFTInterp_Poles& S) { (void) S; }
  void operator=(const FFTInterp_Poles& S) { (void) S; }
};


// A reterpolator between two quadratures using FFT interpolation on each path
// Generalized to handle any quadrature
class FFTInterp_NoPoles
{
  Quadrature* q1;
  Quadrature* q2;

  // q1 Data
  int N_rows1, N_theta1, N_phi1;

  // q2 Data
  int N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  int N_rows, N_phi;

  // Intermediate Storage  N_rows x N_phi
  vector<complex> T;
  // Intermediate Storage  N_theta
  vector<complex> t;

  // First stage phi (row) FFTs
  vector<fftw_plan> rFFT1;
  fftw_plan rIFFT1;

  // Second stage theta (column) FFTs
  fftw_plan tFFT;
  fftw_plan tIFFT;

  // Third stage phi (row) FFTs
  fftw_plan rFFT2;
  vector<fftw_plan> rIFFT2;

 public:

  // Constructor
  FFTInterp_NoPoles( Quadrature* q1_, Quadrature* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_theta1(2*N_rows1-2), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_theta2(2*N_rows2-2), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ), N_phi( max(N_phi1,N_phi2) ),
        T( N_phi * N_rows ), t( max(N_theta1,N_theta2) )
  {
    fftw_complex* T0 = reinterpret_cast<fftw_complex*>(&T[0]);
    fftw_complex* t0 = reinterpret_cast<fftw_complex*>(&t[0]);

    // Define an phi1 FFT for each row
    rFFT1 = vector<fftw_plan>( N_rows1 );
    for( int n = 0; n < N_rows1; ++n ) {
      int N_phi_n = q1->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rFFT1[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_FORWARD,FFTW_MEASURE);
    }

    // Define a phi block IFFT for each row
    rIFFT1 = fftw_plan_many_dft( 1, &N_phi, N_rows1,
                                 T0, &N_phi, 1, N_phi,
                                 T0, &N_phi, 1, N_phi,
                                 FFTW_BACKWARD, FFTW_MEASURE );

    // Define a theta FFT and IFFT (could also do these as blocks...)
    tFFT = fftw_plan_dft_1d(N_theta1, t0, t0, FFTW_FORWARD, FFTW_MEASURE);
    tIFFT = fftw_plan_dft_1d(N_theta2, t0, t0, FFTW_BACKWARD, FFTW_MEASURE);

    // Define a phi2 block FFT for each row
    rFFT2 = fftw_plan_many_dft( 1, &N_phi, N_rows2,
                                T0, &N_phi, 1, N_phi,
                                T0, &N_phi, 1, N_phi,
                                FFTW_FORWARD, FFTW_MEASURE );

    // Define an phi2 IFFT for each row
    rIFFT2 = vector<fftw_plan>( N_rows2 );
    for( int n = 0; n < N_rows2; ++n ) {
      int N_phi_n = q2->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rIFFT2[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_BACKWARD,FFTW_MEASURE);
    }
  }

  // Destructor
  ~FFTInterp_NoPoles() {
    for( int n = 0; n < N_rows1; ++n )
      fftw_destroy_plan( rFFT1[n] );
    fftw_destroy_plan( rIFFT1 );
    fftw_destroy_plan( tFFT );
    fftw_destroy_plan( tIFFT );
    fftw_destroy_plan( rFFT2 );
    for( int n = 0; n < N_rows2; ++n )
      fftw_destroy_plan( rIFFT2[n] );
  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( A.quad == q1 && B.quad == q2 );

    // Copy from A and interpolate each row
    for( int n = 0; n < N_rows1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      // Copy into the kth row of T (scale for the upcoming interpolations)
      for( int m = 0; m < N_phi_n; ++m ) {
        Tn[m] = A.C[ rown.getPoint(m).index ] / double(N_phi_n*N_theta1*N_phi);
      }

      fftw_execute( rFFT1[n] );                    // FFT_phi
      f_cut( Tn, N_phi_n, N_phi );                 // Smooth phi
    }

    fftw_execute( rIFFT1 );                        // IFFT_phi all rows

    // Interpolate from N_theta1 to N_theta2
    for( int m = 0; m < N_phi/2; ++m ) {

      // Unwrap into t
      t[0] = T[m         + 0*N_phi];
      for( int n = 1; n < N_rows1; ++n ) {
        t[n]          = T[m         + n*N_phi];
        t[N_theta1-n] = T[m+N_phi/2 + n*N_phi];
      }

      // Interpolate from N_theta1 to N_theta2
      fftw_execute( tFFT );                        // FFT_theta
      f_cut( &t[0], N_theta1, N_theta2 );          // Smooth theta
      //t[N_theta2/2] = 0;                           // Zero oddball?
      fftw_execute( tIFFT );                       // IFFT_theta

      // Unwrap back into T
      T[m + 0*N_phi] = T[m+N_phi/2 + 0*N_phi] = t[0];
      for( int n = 1; n < N_rows2; ++n ) {
        T[m         + n*N_phi] = t[n];
        T[m+N_phi/2 + n*N_phi] = t[N_theta2-n];
      }
    }

    fftw_execute( rFFT2 );                         // FFT_phi all rows

    // Interpolate each row and copy into B
    for( int n = 0; n < N_rows2; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      f_cut( Tn, N_phi, N_phi_n );                 // Smooth phi
      //if( N_phi_n != 1 ) Tn[N_phi_n/2] = 0;      // Zero oddball?
      fftw_execute( rIFFT2[n] );                   // IFFT_phi

      // Copy into B
      for( int m = 0; m < N_phi_n; ++m ) {
        B.C[ rown.getPoint(m).index ] = Tn[m];
      }
    }
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_NoPoles(const FFTInterp_NoPoles& S) { (void) S; }
  void operator=(const FFTInterp_NoPoles& S) { (void) S; }
};


// A reterpolator between two quadratures using FFT interpolation
// Optimized to reduce the number of FFTs
// Generalized to handle any quadrature
class FFTInterp_FFT2
{
  Quadrature* q1;
  Quadrature* q2;

  // q1 Data
  int N_rows1, N_theta1, N_phi1;

  // q2 Data
  int N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  int N_rows, N_phi, N_phi0;

  // Intermediate Storage  N_rows x N_phi
  vector<complex> T;
  // Intermediate Storage  N_theta
  vector<complex> t;

  // First stage phi (row) FFTs
  vector<fftw_plan> rFFT1;

  // Second stage theta (column) FFTs
  fftw_plan tFFT;
  fftw_plan tIFFT;

  // Third stage phi (row) FFTs
  vector<fftw_plan> rIFFT2;

 public:

  // Constructor
  FFTInterp_FFT2( Quadrature* q1_, Quadrature* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_theta1(2*N_rows1-2), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_theta2(2*N_rows2-2), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ),
        N_phi( max(N_phi1,N_phi2) ), N_phi0( min(N_phi1,N_phi2) ),
        T( N_phi * N_rows ), t( max(N_theta1,N_theta2) )
  {
    fftw_complex* T0 = reinterpret_cast<fftw_complex*>(&T[0]);
    fftw_complex* t0 = reinterpret_cast<fftw_complex*>(&t[0]);

    // Define an phi1 FFT for each row
    rFFT1 = vector<fftw_plan>( N_rows1 );
    for( int n = 0; n < N_rows1; ++n ) {
      int N_phi_n = q1->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rFFT1[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_FORWARD,FFTW_MEASURE);
    }

    // Define a theta FFT and IFFT (could also do these as blocks...)
    tFFT = fftw_plan_dft_1d(N_theta1, t0, t0, FFTW_FORWARD, FFTW_MEASURE);
    tIFFT = fftw_plan_dft_1d(N_theta2, t0, t0, FFTW_BACKWARD, FFTW_MEASURE);

    // Define an phi2 IFFT for each row
    rIFFT2 = vector<fftw_plan>( N_rows2 );
    for( int n = 0; n < N_rows2; ++n ) {
      int N_phi_n = q2->getRow(n).size();
      fftw_complex* Tn = T0 + n*N_phi;
      rIFFT2[n] = fftw_plan_dft_1d(N_phi_n,Tn,Tn,FFTW_BACKWARD,FFTW_MEASURE);
    }
  }

  // Destructor
  ~FFTInterp_FFT2() {
    for( int n = 0; n < N_rows1; ++n )
      fftw_destroy_plan( rFFT1[n] );
    fftw_destroy_plan( tFFT );
    fftw_destroy_plan( tIFFT );
    for( int n = 0; n < N_rows2; ++n )
      fftw_destroy_plan( rIFFT2[n] );
  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( A.quad == q1 && B.quad == q2 );

    // Copy from A and interpolate each row
    for( int n = 0; n < N_rows1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      // Copy into the kth row of T (scale for the upcoming interpolations)
      for( int m = 0; m < N_phi_n; ++m ) {
        Tn[m] = A.C[ rown.getPoint(m).index ] / double(N_phi_n*N_theta1);
      }

      fftw_execute( rFFT1[n] );                    // FFT_phi
      f_cut( Tn, N_phi_n, N_phi0 );                // Smooth phi
    }

    // Interpolate from N_theta1 to N_theta2
    for( int m = 0; m < N_phi0; ++m ) {

      int one = ISODD(m) ? -1 : 1;

      // Unwrap into t
      t[0] = T[m + 0*N_phi];
      for( int n = 1; n < N_rows1; ++n ) {
        t[n]          =       T[m + n*N_phi];
        t[N_theta1-n] = one * T[m + n*N_phi];
      }

      //cerr << m << "\t" << t[N_theta1/2] << endl;
      //assert( !ISODD(m) || abs(t[N_theta1/2]) == 0 );
      //if( ISODD(m) ) t[N_theta1/2] = 0;

      // Interpolate from N_theta1 to N_theta2
      fftw_execute( tFFT );                        // FFT_theta
      f_cut( &t[0], N_theta1, N_theta2 );          // Smooth theta
      //t[N_theta2/2] = 0;                           // Zero oddball?
      fftw_execute( tIFFT );                       // IFFT_theta

      // Unwrap back into T
      for( int n = 0; n < N_rows2; ++n ) {
        T[m + n*N_phi] = t[n];
      }
    }

    // Interpolate each row and copy into B
    for( int n = 0; n < N_rows2; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T[0] + n*N_phi;

      f_cut( Tn, N_phi0, N_phi_n );                // Smooth phi
      //if( N_phi_n != 1 ) Tn[N_phi_n/2] = 0;      // Zero oddball?
      fftw_execute( rIFFT2[n] );                   // IFFT_phi

      // Copy into B
      for( int m = 0; m < N_phi_n; ++m ) {
        B.C[ rown.getPoint(m).index ] = Tn[m];
      }
    }
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_FFT2(const FFTInterp_FFT2& S) { (void) S; }
  void operator=(const FFTInterp_FFT2& S) { (void) S; }
};



// A reterpolator between two quadratures using FFT interpolation
// Optimized to reduce the number of FFTs
// Generalized to handle any quadrature
class FFTInterp_FFT2_Fast
{
  Quadrature* q1;
  Quadrature* q2;

  // q1 Data
  int N_rows1, N_theta1, N_phi1;

  // q2 Data
  int N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  int N_rows, N_theta, N_phi, N_phi0;

  // Intermediate Storage  N_rows x N_theta
  vector<complex> T1;
  vector<complex> T2;

  // First stage phi (row) FFTs
  vector<fftw_plan> rFFT;

  // Second stage theta (column) FFTs
  fftw_plan cFFT;
  fftw_plan cIFFT;

  // Third stage phi (row) FFTs
  vector<fftw_plan> rIFFT;

 public:

  // Constructor
  FFTInterp_FFT2_Fast( Quadrature* q1_, Quadrature* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_theta1(2*N_rows1-2), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_theta2(2*N_rows2-2), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ),
        N_theta( max(N_theta1,N_theta2) ),
        N_phi( max(N_phi1,N_phi2) ), N_phi0( min(N_phi1,N_phi2) ),
        T1( N_phi * N_theta ), T2( N_phi * N_theta )
  {
    fftw_complex* t1 = reinterpret_cast<fftw_complex*>(&T1[0]);
    fftw_complex* t2 = reinterpret_cast<fftw_complex*>(&T2[0]);

    // Define an phi1 FFT for each row
    rFFT = vector<fftw_plan>( N_rows1 );
    for( int n = 0; n < N_rows1; ++n ) {
      int N_phi_n = q1->getRow(n).size();
      rFFT[n] = fftw_plan_dft_1d(N_phi_n, t1+n*N_phi, t2+n*N_phi,
                                 FFTW_FORWARD, FFTW_PATIENT|FFTW_DESTROY_INPUT);
    }

    // Define a theta FFT and IFFT as blocks
    cFFT = fftw_plan_many_dft( 1, &N_theta1, N_phi0,
                               t2, &N_theta, N_phi, 1,
                               t1, &N_theta, N_phi, 1,
                               FFTW_FORWARD, FFTW_PATIENT|FFTW_DESTROY_INPUT);

    cIFFT = fftw_plan_many_dft( 1, &N_theta2, N_phi0,
                                t2, &N_theta, N_phi, 1,
                                t1, &N_theta, N_phi, 1,
                                FFTW_BACKWARD, FFTW_PATIENT|FFTW_DESTROY_INPUT);

    // Define an phi2 IFFT for each row
    rIFFT = vector<fftw_plan>( N_rows2 );
    for( int n = 0; n < N_rows2; ++n ) {
      int N_phi_n = q2->getRow(n).size();
      rIFFT[n] = fftw_plan_dft_1d(N_phi_n, t1+n*N_phi, t2+n*N_phi,
                                  FFTW_BACKWARD, FFTW_PATIENT|FFTW_DESTROY_INPUT);
    }
  }

  // Destructor
  ~FFTInterp_FFT2_Fast() {
    for( int n = 0; n < N_rows1; ++n )
      fftw_destroy_plan( rFFT[n] );
    fftw_destroy_plan( cFFT );
    fftw_destroy_plan( cIFFT );
    for( int n = 0; n < N_rows2; ++n )
      fftw_destroy_plan( rIFFT[n] );
  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( A.quad == q1 && B.quad == q2 );

    // Copy from A and interpolate each row
    for( int n = 0; n < N_rows1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T1[0] + n*N_phi;

      // Copy into the kth row of T (scale for the upcoming interpolations)
      for( int m = 0; m < N_phi_n; ++m ) {
        Tn[m] = A.C[ rown.getPoint(m).index ] / double(N_phi_n*N_theta1);
      }

      fftw_execute( rFFT[n] );                    // FFT_phi   t1->t2
      f_cut( &T2[0]+n*N_phi, N_phi_n, N_phi0 );   // Smooth phi
    }

    // Interpolate from N_theta1 to N_theta2

    // Apply Spherical Fourier Symmetry
    for( int n = 0; n < N_rows1; ++n ) {
      complex* Tn = &T2[0] + n*N_phi;
      complex* TN = &T2[0] + (N_theta1-n)*N_phi;
      // N_phi0 is always even - unroll the loop
      for( int m = 0; m < N_phi0; m += 2 ) {
        TN[m  ] =        Tn[m  ];
        TN[m+1] = -1.0 * Tn[m+1];
      }
    }

    // Block FFT   t2 -> t1
    fftw_execute( cFFT );

    // Zero T2
    //memset( &T2[0], 0, N_phi*N_theta*sizeof(complex) );

    // Interp each column
    if( N_theta1 > N_theta2 ) {
      int M = N_theta2/2;                       // Note integer division
      int k = M + ISODD(N_theta2);
      int index = N_theta1 - M;
      memcpy( &T2[0], &T1[0], k*N_phi*sizeof(complex) );
      /*
        for( int k1 = 0; k1 < k; ++k1 ) {
        complex* Tk = &T2[0] + k1*N_phi;
        complex* Ti = &T1[0] + k1*N_phi;
        for( int m = 0; m < N_phi0; ++m ) {
        Tk[m] = Ti[m];
        }
        }
      */
      for( ; k < N_theta2; ++k, ++index ) {
        //a[k*S] = a[index*S];
        //a[index*S] = 0;
        complex* Tk = &T2[0] + k*N_phi;
        complex* Ti = &T1[0] + index*N_phi;
        for( int m = 0; m < N_phi0; ++m ) {
          Tk[m] = Ti[m];
        }
      }
    } else if( N_theta1 < N_theta2 ) {
      int end = N_theta1/2 - !ISODD(N_theta1);       // Note integer division
      int index = N_theta1 - 1;
      int k = N_theta2 - 1;
      memcpy( &T2[0], &T1[0], (end+1)*N_phi*sizeof(complex) );
      memset( &T2[0] + (end+1)*N_phi, 0, (k-end-1)*N_phi*sizeof(complex) );
      /*
        for( int k1 = 0; k1 <= end; ++k1 ) {
        complex* Tk = &T2[0] + k1*N_phi;
        complex* Ti = &T1[0] + k1*N_phi;
        for( int m = 0; m < N_phi0; ++m ) {
        Tk[m] = Ti[m];
        }
        }
      */
      for( ; index > end; --k, --index ) {
        //a[k*S] = a[index*S];
        //a[index*S] = 0;
        complex* Tk = &T2[0] + k*N_phi;
        complex* Ti = &T1[0] + index*N_phi;
        for( int m = 0; m < N_phi0; ++m ) {
          Tk[m] = Ti[m];
        }
      }
    } // else they're equal and do nothing

    // Block IFFT  t2 -> t1
    fftw_execute( cIFFT );

    // Interpolate each row and copy into B
    for( int n = 0; n < N_rows2; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      int N_phi_n = rown.size();


      f_cut( &T1[0] + n*N_phi, N_phi0, N_phi_n );  // Smooth phi
      //if( N_phi_n != 1 ) Tn[N_phi_n/2] = 0;      // Zero oddball?
      fftw_execute( rIFFT[n] );                   // IFFT_phi   t2 -> t1

      // Copy into B
      complex* Tn = &T2[0] + n*N_phi;
      for( int m = 0; m < N_phi_n; ++m ) {
        B.C[ rown.getPoint(m).index ] = Tn[m];
      }
    }
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_FFT2_Fast(const FFTInterp_FFT2_Fast& S) { (void) S; }
  void operator=(const FFTInterp_FFT2_Fast& S) { (void) S; }
};


// A reterpolator between two quadratures using FFT interpolation
// Optimized to reduce the number of FFTs
// Generalized to handle any quadrature
class FFTInterp_FFT2_Fast2
{
  Quadrature* q1;
  Quadrature* q2;

  // q1 Data
  int N_rows1, N_theta1, N_phi1;

  // q2 Data
  int N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  int N_rows, N_theta, N_phi, N_phi0;

  // Intermediate Storage  N_rows x N_theta
  vector<complex> T1;
  vector<complex> T2;

  // First stage phi (row) FFTs
  vector<fftw_plan> rFFT;

  // Second stage theta (column) FFTs
  fftw_plan cFFT;
  fftw_plan cIFFT;

  // Third stage phi (row) FFTs
  vector<fftw_plan> rIFFT;

 public:

  // Constructor
  FFTInterp_FFT2_Fast2( Quadrature* q1_, Quadrature* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_theta1(2*N_rows1-2), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_theta2(2*N_rows2-2), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ),
        N_theta( max(N_theta1,N_theta2) ),
        N_phi( max(N_phi1,N_phi2) ), N_phi0( min(N_phi1,N_phi2) ),
        T1( N_phi * N_theta ),
        T2( N_phi * N_theta )
  {
    fftw_complex* t1 = reinterpret_cast<fftw_complex*>(&T1[0]);
    fftw_complex* t2 = reinterpret_cast<fftw_complex*>(&T2[0]);

    // Define an phi1 FFT for each row
    rFFT = vector<fftw_plan>( N_rows1 );
    for( int n = 0; n < N_rows1; ++n ) {
      int N_phi_n = q1->getRow(n).size();
      rFFT[n] = fftw_plan_dft_1d( N_phi_n, t1+n*N_phi, t2+n*N_phi,
                                  FFTW_FORWARD,
                                  FFTW_PATIENT | FFTW_DESTROY_INPUT );
    }

    // Define a theta FFT and IFFT as blocks (with transpose)
    cFFT = fftw_plan_many_dft( 1, &N_theta1, N_phi0,
                               t2, NULL, N_phi, 1,
                               t1, NULL, 1, N_theta,
                               FFTW_FORWARD,
                               FFTW_PATIENT | FFTW_DESTROY_INPUT );

    cIFFT = fftw_plan_many_dft( 1, &N_theta2, N_phi0,
                                t1, NULL, 1, N_theta,
                                t2, NULL, N_phi, 1,
                                FFTW_BACKWARD,
                                FFTW_PATIENT | FFTW_DESTROY_INPUT );

    // Define an phi2 IFFT for each row
    rIFFT = vector<fftw_plan>( N_rows2 );
    for( int n = 0; n < N_rows2; ++n ) {
      int N_phi_n = q2->getRow(n).size();
      rIFFT[n] = fftw_plan_dft_1d( N_phi_n, t2+n*N_phi, t1+n*N_phi,
                                   FFTW_BACKWARD,
                                   FFTW_PATIENT | FFTW_DESTROY_INPUT );
    }
  }

  // Destructor
  ~FFTInterp_FFT2_Fast2() {
    for( int n = 0; n < N_rows1; ++n )
      fftw_destroy_plan( rFFT[n] );
    fftw_destroy_plan( cFFT );
    fftw_destroy_plan( cIFFT );
    for( int n = 0; n < N_rows2; ++n )
      fftw_destroy_plan( rIFFT[n] );
  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( A.quad == q1 && B.quad == q2 );

    // Copy from A and interpolate each row
    for( int n = 0; n < N_rows1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      int N_phi_n = rown.size();
      complex* Tn = &T1[0] + n*N_phi;

      // Copy into the kth row of T (scale for the upcoming interpolations)
      // Assume row-major storage
      int start_index = rown.getPoint(0).index;
      complex* AC = &A.C[0] + start_index;
      for( int m = 0; m < N_phi_n; ++m ) {
        //assert( m+start_index == rown.getPoint(m).index );
        Tn[m] = AC[m] / double(N_phi_n*N_theta1);
      }

      fftw_execute( rFFT[n] );                // FFT_phi   T1->T2
      f_cut( &T2[0]+n*N_phi, N_phi_n, N_phi0 );   // Smooth phi
    }

    // Interpolate from N_theta1 to N_theta2

    // Apply Spherical Fourier Symmetry
    for( int n = 0; n < N_rows1; ++n ) {
      complex* Tn = &T2[0] + n*N_phi;
      complex* TN = &T2[0] + (N_theta1-n)*N_phi;
      // N_phi0 is always even - unroll the loop
      for( int m = 0; m < N_phi0; m += 2 ) {
        TN[m  ] =        Tn[m  ];
        TN[m+1] = -1.0 * Tn[m+1];
      }
    }

    // Block FFT   T2 -> T1 (with transpose)
    fftw_execute( cFFT );

    // N0 = N_theta1
    // N = N_theta2
    if( N_theta2 > N_theta1 ) {
      // Cut all thetas   (inlined and optimized)
      int end = N_theta1/2 - !ISODD(N_theta1);
      for( int m = 0; m < N_phi0; ++m ) {
        //f_extend( &T1[0] + m*N_theta, N_theta1, N_theta2 );
        complex* T = &T1[0] + m*N_theta;
        int index = N_theta1 - 1;
        int q = N_theta2 - 1;
        for( ; index > end; --q, --index ) {
          T[q] = T[index];
          T[index] = 0;
        }
        // Zero out anything that was missed
        for( ; q >= N_theta1; --q ) {
          T[q] = 0;
        }
      }
    } else if( N_theta2 < N_theta1 ) {
      // Cut all thetas   (inlined and optimized)
      int M = N_theta2/2;
      for( int m = 0; m < N_phi0; ++m ) {
        //f_trunc( &T1[0] + m*N_theta, N_theta1, N_theta2 );
        complex* T = &T1[0] + m*N_theta;
        int k = M + ISODD(N_theta2);
        int index = N_theta1 - M;
        for( ; k < N_theta2; ++k, ++index ) {
          T[k] = T[index];
          T[index] = 0;
        }
        // Zero out anything that was missed
        for( ; k < N_theta1 - M; ++k ) {
          T[k] = 0;
        }
      }
    } // else they're equal and do nothing



    // Block IFFT  T1 -> T2 (with transpose)
    fftw_execute( cIFFT );

    // Interpolate each row and copy into B
    for( int n = 0; n < N_rows2; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      int N_phi_n = rown.size();


      f_cut( &T2[0] + n*N_phi, N_phi0, N_phi_n );  // Smooth phi
      //if( N_phi_n != 1 ) Tn[N_phi_n/2] = 0;      // Zero oddball?
      fftw_execute( rIFFT[n] );                   // IFFT_phi   T2 -> T1

      // Copy into B
      complex* Tn = &T1[0] + n*N_phi;
      int start_index = rown.getPoint(0).index;
      complex* BC = &B.C[0] + start_index;
      for( int m = 0; m < N_phi_n; ++m ) {
        //assert( m+start_index == rown.getPoint(m).index );
        BC[m] = Tn[m];
      }
    }
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_FFT2_Fast2(const FFTInterp_FFT2_Fast2& S) { (void) S; }
  void operator=(const FFTInterp_FFT2_Fast2& S) { (void) S; }
};




class S2Interp
{
  Quadrature_S2* q1;
  Quadrature_S2* q2;

  // q1 Data
  int N_rows1, N_phi1;

  // q2 Data
  int N_rows2, N_phi2;

  // Intermediate Data
  int N_rows, N_phi;

  // Temp Storage
  vector<complex> T1;
  vector<complex> T2;

  // FFTs
  fftw_plan rIFFT1;
  fftw_plan rFFT2;

  // Polynomial interpolation matrix
  vector<complex> P;

 public:

  S2Interp( Quadrature_S2* q1_, Quadrature_S2* q2_ )
      : q1(q1_), q2(q2_),
        N_rows1(q1->numRows()), N_phi1(q1->maxRowSize()),
        N_rows2(q2->numRows()), N_phi2(q2->maxRowSize()),
        N_rows( max(N_rows1,N_rows2) ), N_phi( max(N_phi1,N_phi2) ),
        T1( N_phi * N_rows ), T2( N_phi * N_rows ),
        P( N_rows1 * N_rows2 )
  {
    fftw_complex* Ta = reinterpret_cast<fftw_complex*>(&T1[0]);
    fftw_complex* Tb = reinterpret_cast<fftw_complex*>(&T2[0]);

    // Define a phi1 block IFFT for each row
    rIFFT1 = fftw_plan_many_dft( 1, &N_phi, N_rows1,
                                 Ta, &N_phi, 1, N_phi,
                                 Tb, &N_phi, 1, N_phi,
                                 FFTW_BACKWARD, FFTW_MEASURE );

    // Define a phi2 block FFT for each row
    rFFT2 = fftw_plan_many_dft( 1, &N_phi, N_rows2,
                                Ta, &N_phi, 1, N_phi,
                                Tb, &N_phi, 1, N_phi,
                                FFTW_FORWARD, FFTW_MEASURE );

    // Define the polynomial interpolation matrix Pij
    // We interpolate from x_i = cos(theta_i) to x_j = cos(theta_j)

    // Define x_i
    vector<double> x1(N_rows1);
    for( int n = 0; n < N_rows1; ++n )
      x1[n] = cos( q1->getRow(n).getTheta() );

    // Define x_j
    vector<double> x2(N_rows2);
    for( int n = 0; n < N_rows2; ++n )
      x2[n] = cos( q2->getRow(n).getTheta() );

    // Compute Pij (slow way)
    for( int i = 0; i < N_rows2; ++i ) {
      for( int j = 0; j < N_rows1; ++j ) {
        double Pij = 1;
        for( int k = 0; k < j; ++k ) {
          Pij *= (x2[i] - x1[k]) * (x1[j] - x1[k]);
        }
        for( int k = j+1; k < N_rows1; ++k ) {
          Pij *= (x2[i] - x1[k]) * (x1[j] - x1[k]);
        }
        P[i*N_rows1 + j] = Pij;
      }
    }

  }

  inline void apply( NFunction& A, NFunction& B )
  {
    assert( (Quadrature_S2*) A.quad == q1 && (Quadrature_S2*) B.quad == q2 );

    // Copy the data from A into a matrix
    for( int n = 0; n < N_rows1; ++n ) {
      QuadratureRow& rown = q1->getRow(n);
      complex* Tn = &T1[0] + n*N_phi;

      // Copy into the kth row of T
      for( int m = 0; m < N_phi1; ++m )
        Tn[m] = A.C[ rown.getPoint(m).index ] / N_phi1;
    }

    // Take the FFT along every row and T1 -> T2
    fftw_execute( rIFFT1 );

    // Fourier truncation and polynomial interpolation
    if( N_phi1 > N_phi2 ) {
      // Truncate the Fourier coefficients
      for( int n = 0; n < N_rows1; ++n ) {
        complex* Tn = &T2[0] + n*N_phi;
        f_cut( Tn, N_phi1, N_phi2 );
      }

      // Apply the polynomial interpolation
      /*
        double scale = 1.0;
        cblas_zgemm(CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        N_rows2, N_rows1, &scale,
        &P[0], N_rows2,
        &T2[0], 1, &scale,
        &T1[0], 1);
      */
      for( int i = 0; i < N_rows2; ++i ) {
        for( int j = 0; j < N_phi2; ++j ) {
          complex Tij = 0;
          for( int k = 0; k < N_rows1; ++k ) {
            Tij += P[i*N_rows1+k] * T2[k*N_phi+j];
          }
          T1[i*N_phi+j] = Tij;
        }
      }
    } else {
      // Apply the polynomial interpolation
      /*
        double scale = 1.0;
        cblas_zgemv(CblasRowMajor,
        CblasNoTrans,
        N_rows2, N_rows1, &scale,
        &P[0], N_rows2,
        &T2[0], 1, &scale,
        &T1[0], 1);
      */
      for( int i = 0; i < N_rows2; ++i ) {
        for( int j = 0; j < N_phi1; ++j ) {
          complex Tij = 0;
          for( int k = 0; k < N_rows1; ++k ) {
            Tij += P[i*N_rows1+k] * T2[k*N_phi+j];
          }
          T1[i*N_phi+j] = Tij;
        }
      }

      // Pad the Fourier coefficients
      for( int n = 0; n < N_rows1; ++n ) {
        complex* Tn = &T1[0] + n*N_phi;
        f_cut( Tn, N_phi1, N_phi2 );
      }
    }

    // Take the IFFT along every row and T1 -> T2
    fftw_execute( rFFT2 );

    // Copy the data into B
    for( int n = 0; n < N_rows2; ++n ) {
      QuadratureRow& rown = q2->getRow(n);
      complex* Tn = &T2[0] + n*N_phi;

      // Copy into B
      for( int m = 0; m < N_phi2; ++m ) {
        B.C[ rown.getPoint(m).index ] = Tn[m];
      }
    }


  }

 private:
  // Disable Copy and Assignment
  S2Interp(const S2Interp& S) { (void) S; }
  void operator=(const S2Interp& S) { (void) S; }
};



typedef FFTInterp_Poles Interpolator;
typedef FFTInterp_Poles Anterpolator;

//typedef FFTInterp_NoPoles Interpolator;
//typedef FFTInterp_NoPoles Anterpolator;

//typedef FFTInterp_FFT2 Interpolator;
//typedef FFTInterp_FFT2 Anterpolator;


#endif
