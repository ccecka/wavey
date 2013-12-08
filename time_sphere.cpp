const double MLFMM_ALPHA(1.0);

#include "General.hpp"

#include "MLFMM_Env.hpp"
#include "Direct.hpp"

double printResults(const vector<complex>& mlfmm,
		    const vector<complex>& exact)
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
  inline complex operator()(double r) {
    return exp(CI * kappa * r) / r;
  }
};

// Big Run - Many tests keeping points/wavelength constant
int main(int argc, char** argv)
{
  (void) argc; (void) argv;

  double a0 = 2;

  int Nnum = 15;
  int Nlist[15] = {1000, 2000, 4000, 8000, 16000,
		   32000, 64000, 128000, 256000, 512000,
		   1024000, 2048000, 4096000, 8192000, 16384000};

  fstream myFile("time_sphere.dat", ios::out);
  myFile << "N\tkappa\tDIRtime\tFMMtime2\tFMMtime3\tFMMtime4\tFMMtime5\tFMMtime6\tFMMtime7\tFMMtime8\tFMMtime9" << endl;
  //myFile << "N\ta0kappa\tDIRtime\tFMMtime6\tFMMtime7\tFMMtime8" << endl;

  for (int Nindex = 0; Nindex < Nnum; ++Nindex) {
    int N = Nlist[Nindex];
    srand48(0);

    vector<Vec3> p(N);
    vector<complex> psi(N);

    for (int k = 0; k < N;) {
      Vec3 pk = Vec3(getRandom(-a0/2,a0/2),
                     getRandom(-a0/2,a0/2),
                     getRandom(-a0/2,a0/2));
      // Sphere selection
      if (pk.mag() <= 1) {
        pk.normalize();
        p[k] = pk;
        psi[k] = 1;
        ++k;
      }
    }

    HelmKernel K;
    K.kappa = (400/a0) * pow(N/8192000.0, 1.0/2.0);

    myFile << N << "\t" << a0*K.kappa << "\t";

    vector<complex> exact(N,0);
    if (N <= 200000) {
    //if( N <= 1 ) {
      cout << "Computing Direct " << N << endl;

      StopWatch timer;
      timer.start();

      Direct(K, p, psi, exact);

      double DIRtime = timer.stop();
      myFile << DIRtime << "\t";
    } else {
      myFile << "NaN\t";
    }


    for (int levels = 2; levels <= 9; ++levels) {
      if(N > 500000 && levels <= 2) { myFile << "NaN\t"; continue; }
      if(N > 1000000 && levels <= 3) { myFile << "NaN\t"; continue; }
      if(N > 2000000 && levels <= 4) { myFile << "NaN\t"; continue; }
      if(N > 4000000 && levels <= 5) { myFile << "NaN\t"; continue; }
      if(N > 8000000 && levels <= 6) { myFile << "NaN\t"; continue; }
      if(N > 16000000 && levels <= 7) { myFile << "NaN\t"; continue; }

      cout << "\nStarting " << N << "\t" << K.kappa << "\t" << levels << endl;
      const static int DIM = 3;
      MLFMM_Env<HelmKernel,DIM> env(K, p, levels, 1e-4);

      vector<complex> mlfmm(N,0);

      StopWatch timer;
      timer.start();
      env.execute(psi, mlfmm);
      env.execute(psi, mlfmm);
      env.execute(psi, mlfmm);
      env.execute(psi, mlfmm);
      double MLFMMtime = timer.stop() / 4;

      double aveErr = 0;
      if(exact[0] != complex(0) || exact[1] != complex(0))
	aveErr = printResults(mlfmm, exact);

      cout << "N\ta0kappa\tlevels\tFMMtime\tAveErr" << endl;
      cout << N << "\t" << a0*K.kappa << "\t" << levels << "\t" << MLFMMtime << "\t" << "\t" << aveErr << endl;

      myFile << MLFMMtime << "\t";
    }

    myFile << endl;
  }

  myFile.close();
  return 0;
}
