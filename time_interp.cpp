#include "General.hpp"

#include "MLFMM_Env.hpp"
#include "Direct.hpp"

double printResults( const vector<complex>& mlfmm,
                     const vector<complex>& exact )
{
  double TotRelError = 0;

  double TotErrorSq = 0;
  double TotNormSq = 0;

  double MaxRelErr = 0;

  int N = exact.size();
  for( int k = 0; k < N; ++k ) {
    //cerr << mlfmm[k] << "\t\t" << exact[k] << "\t\t" << mlfmm[k]-exact[k] << endl;

    // Individual Relative
    TotRelError += abs(mlfmm[k] - exact[k]) / abs(exact[k]);

    // Total Relative
    TotErrorSq += norm(mlfmm[k] - exact[k]);
    TotNormSq  += norm(exact[k]);

    // Max Absolute
    MaxRelErr = max(abs(mlfmm[k] - exact[k]) / abs(exact[k]), MaxRelErr);
  }

  // Relative Error
  double RelError = sqrt(TotErrorSq/TotNormSq);
  cout << "Tot Rel Error: " << RelError << endl;

  // Average Relative Error
  double AveRelError = TotRelError/N;
  cout << "Ave Rel Error: " << AveRelError << endl;

  // Max Relative Error
  cout << "Max Rel Error: " << MaxRelErr << endl;

  return AveRelError;
}

struct HelmKernel
{
  double kappa;
  inline complex operator()( double r ) { return exp( CI * kappa * r ) / r; }
};

template <typename Interp>
double test_interp(NFunction& f1, NFunction& f2) {
  Interp interp(f1.quad, f2.quad);

  StopWatch timer;
  timer.start();
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  interp.apply(f1, f2);
  double time = timer.stop() / 10;

  return time;
}





// Timing comparison of interpolation and anterpolation stage
int main( int argc, char* argv[] )
{
  (void) argc; (void) argv;

  double boxSize = 1;
  double eps = 1e-6;
  // Set MLFMM_ALPHA = 0.8

  HelmKernel K;

  fstream myFile("time_interp.dat", ios::out);
  myFile << "kappa\tell\tPolyInterp\tFFTInterp\tFFTInterp2\tFFTInterp3" << endl;

  for( double kappa = 0.95; kappa < 10000; kappa *= 1.20 ) {

    K.kappa = kappa;

    // Uniform Quadrature Fourier version
    Quadrature_Uniform* quad_L = new Quadrature_Uniform(K, boxSize, eps);
    Quadrature_Uniform* quad_H = new Quadrature_Uniform(K, 2*boxSize, eps);

    std::cout << "Quadratures Created. Kappa = " << kappa << "\n";

    int ell = quad_H->getTruncation();

    // Create a NFunction
    NFunction F_U_L = Translation_Function( quad_L, Vec3(boxSize,0,0) );


    NFunction F_U_H = NFunction(quad_H);
    double timeOrig = test_interp<FFTInterp_Poles>(F_U_L, F_U_H);
    std::cout << "Done FFTInterp_Poles\n";

    NFunction F_U_H_FFT2 = NFunction(quad_H);
    double time_fft2 = test_interp<FFTInterp_FFT2>(F_U_L, F_U_H_FFT2);
    std::cout << "Done FFTInterp_FFT2\n";

    NFunction F_U_H_Fast = NFunction(quad_H);
    double time_fast = test_interp<FFTInterp_FFT2_Fast>(F_U_L, F_U_H_Fast);
    std::cout << "Done FFTInterp_FFT2_Fast\n";

    NFunction F_U_H_Fast2 = NFunction(quad_H);
    double time_fast2 = test_interp<FFTInterp_FFT2_Fast2>(F_U_L, F_U_H_Fast2);
    std::cout << "Done FFTInterp_FFT2_Fast2\n";

    double error  = 0;
    double error2 = 0;
    double error3 = 0;
    for( int k = 0; k < quad_H->size(); ++k ) {
      error   = max(abs(F_U_H.C[k] - F_U_H_FFT2.C[k]), error);
      error2  = max(abs(F_U_H.C[k] - F_U_H_Fast.C[k]), error2);
      error3  = max(abs(F_U_H.C[k] - F_U_H_Fast2.C[k]), error3);
    }
    cerr << "Error:\t" << error << "\t" << error2 << "\t" << error3 << endl;

    delete quad_L;
    delete quad_H;


    /*
    // Spherical Harmonics version
    Quadrature_S2* quadS2_L = new Quadrature_S2(K, boxSize, eps);
    Quadrature_S2* quadS2_H = new Quadrature_S2(K, 2*boxSize, eps);

    // Create a NFunction
    NFunction F_S2_L = Translation_Function(quadS2_L, Vec3(0.5*boxSize,0,0));

    S2Interp s2interp(quadS2_L, quadS2_H);
    NFunction F_S2_H = NFunction(quadS2_H);

    timer.start();
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    s2interp.apply(F_S2_L, F_S2_H);
    double timeS2 = timer.stop() / 10;

    delete quadS2_L;
    delete quadS2_H;
    */

    cout << kappa << "\t" << ell << "\t"
        //<< timeS2 << "\t"
         << timeOrig << "\t"
         << time_fft2 << "\t"
         << time_fast << "\t"
        << time_fast2 << "\t"
         << endl;

    myFile << kappa << "\t" << ell << "\t"
        //<< timeS2 << "\t"
           << timeOrig << "\t"
           << time_fft2 << "\t"
           << time_fast << "\t"
           << time_fast2 << "\t"
           << endl;
  }

  return 0;
}
