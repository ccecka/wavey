#ifndef GENERAL_H
#define GENERAL_H

// A general library of useful methods and support packages

// Define to remove all assert statements
//#define NDEBUG
#include <assert.h>


/////////
// STL //
/////////

#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

using namespace std;


// Quick StopWatch
#include <sys/time.h>
struct StopWatch {
  timeval startTime, stopTime, result;
  inline void start() { gettimeofday(&startTime,NULL); }
  inline double stop() { return elapsed(); }
  inline double elapsed() {
    gettimeofday(&stopTime,NULL);
    timersub(&stopTime,&startTime,&result);
    return result.tv_sec + result.tv_usec/1000000.0;   // 10^6 uSec per Sec
  }
};


/////////////////////
// Complex numbers //
/////////////////////

#include <complex>

#define complex complex<double>

// Kill that stupid complex -op- int error
inline complex operator*(const complex& c, const int& n) {return c*double(n);}
inline complex operator*(const int& n, const complex& c) {return c*double(n);}
inline complex operator/(const complex& c, const int& n) {return c/double(n);}
inline complex operator/(const int& n, const complex& c) {return double(n)/c;}


///////////////////////
// Common and Helper //
///////////////////////

// Program Constants
const double  PI(M_PI);
const complex CI(0,1);

#define ISODD(x) ((x) & 1)
#define ISNAN(x) ((x) != (x))
// Ceil/Round/Floor to nearest multiple of K
inline int ceil ( double x, unsigned int K ) { return K*((int)ceil(x/K));  }
inline int round( double x, unsigned int K ) { return K*((int)round(x/K)); }
inline int floor( double x, unsigned int K ) { return K*((int)floor(x/K)); }
// Random number in (0,1)
inline double getRandom() { return drand48(); }
// Random number in (A,B)
inline double getRandom( double A, double B ) { return A + (B-A)*getRandom(); }


//////////
/// IO ///
//////////

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

// Overloaded complex output
inline ostream& operator<<(ostream& os, const complex a)
{
  //ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
  //ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

  int oldp = os.precision(6);

  os << real(a);
  if( imag(a) != 0 ) {
    if( imag(a) < 0 )
      os << " - " << -imag(a) << "*i";
    else
      os << " + " << imag(a) << "*i";
  }

  //os.setf(olda,ios::adjustfield);
  //os.setf(oldf,ios::floatfield);
  os.precision(oldp);
  os << "";

  return os;
}

template <class T>
ostream& operator<<(ostream& os, const vector<T>& a)
{
  ios::fmtflags olda = os.setf(ios::right,ios::adjustfield);
  ios::fmtflags oldf = os.setf(ios::scientific,ios::floatfield);

  int oldp = os.precision(8);

  int N = a.size();
  for( int k = 0; k < N; ++k ) {
    os << k << "\t" << a[k] << "\n";
  }

  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp);
  os << "";

  return os;
}

template <class T>
inline string to_string(const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template <class T>
inline void saveArray(T* a, int N, const char* filename)
{
  fstream myFile(filename, ios::out | ios::binary);
  myFile.write((char*)a, sizeof(T)*N);
  myFile.close();
}

template <class T>
inline void saveVector( const vector<T> a, const char* filename )
{
  saveArray( &a[0], a.size(), filename );
}

template <class T>
inline T* readArray(int N, const char* filename)
{
  T* a = new T[N];
  fstream myFile(filename, ios::in | ios::binary);
  myFile.read((char*)a, sizeof(T)*N);
  myFile.close();
  return a;
}

template <class T>
inline vector<T> readVector(int N, const char* filename)
{
  vector<T> a(N);
  fstream myFile(filename, ios::in | ios::binary);
  myFile.read((char*)&a[0], sizeof(T)*N);
  myFile.close();
  return a;
}



//////////////////////
// FFTs and Helpers //
//////////////////////

#include <fftw3.h>

// Shifts the zero-frequency component to the center of spectrum.
template <class T>
inline T* fftshift( T* a, int N ) {
  int n1 = N/2;
  T temp;
  if( ISODD(N) ) {      // N is odd
    int index = n1;
    T last = a[0];
    while( index != 0 ) {
      temp = a[index];
      a[index] = last;
      last = temp;
      index = (index + n1) % N;
    }
    a[0] = last;
  } else {          // N is even, just swaps
    for( int k = 0; k < n1; ++k ) {
      temp = a[k];
      a[k] = a[n1+k];
      a[n1+k] = temp;
    }
  }
  return a;
}

// Inverse fftshift
template <class T>
inline T* ifftshift( T* a, int N ) {
  int n1 = N/2;
  T temp;
  if( ISODD(N) ) {      // N is odd
    int index = ++n1;
    T last = a[0];
    while( index != 0 ) {
      temp = a[index];
      a[index] = last;
      last = temp;
      index = (index + n1) % N;
    }
    a[0] = last;
  } else {          // N is even, just swaps
    for( int k = 0; k < n1; ++k ) {
      temp = a[k];
      a[k] = a[n1+k];
      a[n1+k] = temp;
    }
  }
  return a;
}

// Truncates the sequence    N0 > N
// Use stride S
template <class T>
inline T* f_trunc( T* a, const int N0, const int N, const int S = 1 ) {
  int M = N/2;                       // Note integer division
  int k = M + ISODD(N);
  int index = N0 - M;
  for( ; k < N; ++k, ++index ) {
    a[k*S] = a[index*S];
    a[index*S] = 0;
  }
  // Zero out anything that was missed
  for( ; k < N0 - M; ++k ) {
    a[k*S] = 0;
  }
  return a;
}


// Extend the sequence       N > N0
// Use stride S
template <class T>
inline T* f_extend( T* a, const int N0, const int N, const int S = 1 ) {
  int end = N0/2 - !ISODD(N0);       // Note integer division
  int index = N0 - 1;
  int q = N - 1;
  for( ; index > end; --q, --index ) {
    a[q*S] = a[index*S];
    a[index*S] = 0;
  }
  // Zero out anything that was missed
  for( ; q > N0 - 1; --q ) {
    a[q*S] = 0;
  }
  return a;
}

template <class T, class Ts>
inline T* scale( T* a, int N, Ts scale ) {
  for( int k = 0; k < N; ++k ) {
    a[k] *= scale;
  }
  return a;
}

template <class T>
inline T* f_smooth( T* a, const int N0, const int N ) {
  if( N0 > N ) {
    return scale( f_trunc(a, N0, N), N, double(N)/N0 );
  } else if( N0 < N ) {
    return f_extend( scale(a, N0, double(N)/N0), N0, N );
  } // else they're equal and do nothing
  return a;
}

template <class T>
inline T* f_cut( T* a, const int N0, const int N, const int S = 1 ) {
  if( N0 > N ) {
    return f_trunc( a, N0, N, S );
  } else if( N0 < N ) {
    return f_extend( a, N0, N, S );
  } // else they're equal and do nothing
  return a;
}


// Perform a (slow) forward fft
inline complex* fft( complex* in, complex* out, int N )
{
  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_1d( N, a, b, FFTW_FORWARD, FFTW_ESTIMATE );
  fftw_execute( p );
  fftw_destroy_plan( p );
  return out;
}

// Perform a (slow) in-place forward fft
inline complex* fft( complex* in, int N )
{
  return fft(in,in,N);
}

// Perform a (slow) backward fft
inline complex* ifft( complex* in, complex* out, int N )
{
  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_1d( N, a, b, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_execute( p );
  fftw_destroy_plan( p );
  scale( out, N, 1.0/N );
  return out;
}

// Perform a (slow) in-place backward fft
inline complex* ifft( complex* in, int N )
{
  return ifft(in,in,N);
}

inline complex* fftinterp( complex* in, int N0, int N )
{
  if( N0 == N ) return in;

  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_plan pFFT = fftw_plan_dft_1d(N0, a, a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan pIFFT = fftw_plan_dft_1d(N, a, a, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute( pFFT );
  if( N0 > N )       scale( f_trunc(in, N0, N), N, 1.0/N0 );
  else if( N0 < N )  f_extend( scale(in, N0, 1.0/N0), N0, N );
  fftw_execute( pIFFT );

  fftw_destroy_plan( pFFT );
  fftw_destroy_plan( pIFFT );
  return in;
}



//////////////////
// GSL Wrappers //
//////////////////

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


// Override GSL underflow error exit
struct _GSLError_ {
  _GSLError_() { gsl_set_error_handler(&_GSLError_::Handler); }
  static void Handler(const char* msg, const char* file, int line, int err) {
    if( err != GSL_EUNDRFLW ) {
      printf("GSLError %d in %s at %d : %s\n",err,file,line,msg);
      exit(1);
    }
  }
};
// Re-define GSL default error handler when loading the library
static _GSLError_ __GSLError__;


// Cylindrical Bessel function
// Returns 0 in the case of underflow
inline double bessel_J( int n, double x )
{
  gsl_sf_result result;
  int status = gsl_sf_bessel_Jn_e(n,x,&result);
  if( status == GSL_EUNDRFLW ) return 0;
  return result.val;
}

// Spherical Bessel function j
// Returns 0 in the case of underflow
inline double bessel_j( int n, double x )
{
  gsl_sf_result result;
  int status = gsl_sf_bessel_jl_e(n,x,&result);
  if( status == GSL_EUNDRFLW ) return 0;
  return result.val;
}

// Cylindral Bessel function y
// Returns 0 in the case of underflow
inline double bessel_Y( int n, double x )
{
  gsl_sf_result result;
  int status = gsl_sf_bessel_Yn_e(n,x,&result);
  if( status == GSL_EUNDRFLW ) return 0;
  return result.val;
}

// Spherical Bessel function y
// Returns 0 in the case of underflow
inline double bessel_y( int n, double x )
{
  gsl_sf_result result;
  int status = gsl_sf_bessel_yl_e(n,x,&result);
  if( status == GSL_EUNDRFLW ) return 0;
  return result.val;
}

// Spherical Hankel function h
inline complex bessel_h( int n, double x )
{
  return complex( bessel_j(n,x), bessel_y(n,x) );
}

// Legendre Polynomial
// Returns 0 in the case of underflow
inline double legendre_P( int n, double x )
{
  gsl_sf_result result;
  int status = gsl_sf_legendre_Pl_e(n,x,&result);
  if( status == GSL_EUNDRFLW ) return 0;
  return result.val;
}

inline void getGaussLegendreQuad(int N, vector<double>& x, vector<double>& w)
{
  gsl_integration_glfixed_table* GLTable =
    gsl_integration_glfixed_table_alloc(N);

  for( int n = 0; n < (N+1)/2; ++n ) {
    x[n + N/2] = GLTable->x[n];
    w[n + N/2] = GLTable->w[n];
  }
  for( int n = 0; n < N/2; ++n ) {
    x[n] = -x[N-1-n];
    w[n] =  w[N-1-n];
  }
}


#endif
