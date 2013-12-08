#ifndef MLFMM_CPP
#define MLFMM_CPP

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
    cerr << mlfmm[k] << "\t\t" << exact[k] << "\t\t" << mlfmm[k]-exact[k] << endl;

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

#if 0
// Standard Test - Random Points in a Box
int main(int, char**)
{
  StopWatch timer;

  int N = 400;
  vector<Vec3> p(N);
  vector<complex> psi(N);

  HelmKernel K;
  K.kappa = 100;

  srand48(0);
  for( int k = 0; k < N; ++k ) {
    p[k] = Vec3( getRandom(-1,1), getRandom(-1,1), getRandom(-1,1) );
    psi[k] = 1;
    //cout << "p[" << k << "] = " << p[k] << endl;
  }

  int levels = 3;
  double eps = 1e-4;

  timer.start();
  MLFMM_Env<HelmKernel,3> MLFMMEnv(K, p, levels, eps);
  double mlfmm_pre_time = timer.stop();
  cout << endl;

  vector<complex> omega(N);
  timer.start();
  MLFMMEnv.execute( psi, omega );
  double mlfmm_time = timer.stop();
  cout << endl;

  vector<complex> exact(N);
  timer.start();
  Direct( K, p, psi, exact );
  double direct_time = timer.stop();
  cout << endl;

  cout << "Direct Time: " << direct_time << endl;
  cout << "FMMPre Time: " << mlfmm_pre_time << endl;
  cout << "MLFMM  Time: " << mlfmm_time << endl;

  double err = printResults( omega, exact );
  cout << "Error: " << err << endl;

  return 0;
}
#endif


// Two point test
int main( int argc, char* argv[] )
{
  HelmKernel K;
  K.kappa = 100;

  int N = 2;
  vector<Vec3> p(N);
  vector<complex> psi(N);

  p[0] = Vec3(0,0,0);
  psi[0] = 1;
  p[1] = Vec3(1,1,1);
  psi[1] = 1;

  // Create Tree
  const static int DIM = 3;
  int levels = 2;
  double eps = 1e-4;
  MLFMM_Env<HelmKernel,DIM> MLFMMEnv(K, p, levels, eps);

  // HACK: Move points to the center of the boxes
  // Note this is not worst case...
  p[0] = Vec3(0.1,0.1,0.1);
  p[1] = Vec3(0.9,0.9,0.9);
  // Overwrite Env
  MLFMMEnv.sourcefield = p;

  // Compute MLFMM
  vector<complex> omega(N);
  MLFMMEnv.execute( psi, omega );

  // Compute exact
  vector<complex> exact(N);
  Direct( K, p, psi, exact );

  printResults( omega, exact );

  return 0;
}


#if 0

// N Level Convergence test
// Show the convergence of the n-level method in the Quadrature size
int main( int argc, char* argv[] )
{
  HelmKernel K;
  K.kappa = 100;
  const static int DIM = 3;
  // Set MLFMM_ALPHA = ?

  int N = 4;
  vector<Vec3> p(N);
  vector<complex> psi(N);

  // Dummies for the bounding box
  p[0] = Vec3(0,0,0);
  psi[0] = 0;
  p[1] = Vec3(4,4,4);
  psi[1] = 0;

  // With 2 (active) levels, the box size is 0.5
  // Set up two points in the r0 = x worst case
  p[2] = Vec3(1,0,1);
  psi[2] = 1;
  p[3] = Vec3(2+1e-12,0,0);
  psi[3] = 1;

  Vec3 r0 = Vec3(2,0,0);
  Vec3 r = Vec3(-1+1e-12,0,-1);

  // Get exact value
  vector<complex> exact(N);
  Direct(K, p, psi, exact);
  complex Ie = exact[3];

  vector<complex> epsIroot;

  for( int levels = 2; levels < 9; ++levels ) {

    string filename = "Error_VS_Quadrature_" + to_string(levels) + ".dat";
    ofstream myfile;
    myfile.open(filename.c_str());

    myfile << "kappa\teps\teps_T\teps_I\teps_G\tell\tQSize\teps_A" << endl;
    int index = 0;

    for( double eps = 1e-0; eps >= 1e-16; eps *= 0.316227766 ) {

      // Create Tree
      MLFMM_Env<HelmKernel,DIM> MLFMMEnv(K, p, levels, eps);

      // Compute MLFMM
      vector<complex> omega(N,0);
      MLFMMEnv.execute( psi, omega );
      complex I = omega[3];

      int ell = MLFMMEnv.getLevel(2).getQuad()->getTruncation();

      // Gegenbauer value
      complex G = gegenbauer_series(ell, K.kappa, r, r0);

      complex eps_T = I - Ie;
      complex eps_I = G - I;
      complex eps_G = Ie - G;

      int Qsize = MLFMMEnv.getLevel(2).getQuad()->size();

      complex eps_A = 0;
      if( levels == 2 ) {
	epsIroot.push_back( eps_I );
      } else {
	eps_A = eps_I - epsIroot[index];
      }

      cout << "kappa\teps\teps_T\teps_I\teps_G\tell\tQsize\teps_A" << endl;
      cout << K.kappa << "\t" << eps << "\t" << abs(eps_T) << "\t" << abs(eps_I) << "\t" << abs(eps_G) << "\t" << ell << "\t" << Qsize << "\t" << abs(eps_A) << endl << endl;

      myfile << K.kappa << "\t" << eps << "\t" << abs(eps_T) << "\t" << abs(eps_I) << "\t" << abs(eps_G) << "\t" << ell << "\t" << Qsize << "\t" << abs(eps_A) << endl;

      ++index;
    }

    myfile.close();
  }

  return 0;
}

#endif

#if 0

// Error as a function of kappa
// Demonstrate accuracy of direct L and quadrature generation
int main( int argc, char* argv[] )
{
  StopWatch timer;
  timer.start();

  ofstream myfile;
  myfile.open("Error_VS_Kappa.dat");

  double eps = 1e-4;

  double boxSize = 1;

  double norm_r0 = 2 * boxSize;
  double norm_r  = MLFMM_ALPHA * sqrt(3) * boxSize;

  Vec3 r0hat( 1, 0, 0 );
  Vec3 r0 = r0hat * norm_r0;

  vector<Vec3> RCUBE(26);
  RCUBE[0]  = Vec3( 1, 0, 0);  RCUBE[1]  = Vec3(-1, 0, 0);
  RCUBE[2]  = Vec3( 0, 1, 0);  RCUBE[3]  = Vec3( 0,-1, 0);
  RCUBE[4]  = Vec3( 0, 0, 1);  RCUBE[5]  = Vec3( 0, 0,-1);
  RCUBE[6]  = Vec3( 1, 1, 0);  RCUBE[7]  = Vec3( 1,-1, 0);
  RCUBE[8]  = Vec3(-1, 1, 0);  RCUBE[9]  = Vec3(-1,-1, 0);
  RCUBE[10] = Vec3( 1, 0, 1);  RCUBE[11] = Vec3( 1, 0,-1);
  RCUBE[12] = Vec3(-1, 0, 1);  RCUBE[13] = Vec3(-1, 0,-1);
  RCUBE[14] = Vec3( 0, 1, 1);  RCUBE[15] = Vec3( 0, 1,-1);
  RCUBE[16] = Vec3( 0,-1, 1);  RCUBE[17] = Vec3( 0,-1,-1);
  RCUBE[18] = Vec3( 1, 1, 1);  RCUBE[19] = Vec3( 1, 1,-1);
  RCUBE[20] = Vec3( 1,-1, 1);  RCUBE[21] = Vec3(-1, 1, 1);
  RCUBE[22] = Vec3( 1,-1,-1);  RCUBE[23] = Vec3(-1, 1,-1);
  RCUBE[24] = Vec3(-1,-1, 1);  RCUBE[25] = Vec3(-1,-1,-1);

  myfile << "kappa\tTotalE\tIntE\tGegenE\tEBFE\tell\tell_ebf\tQSize" << endl;
  cout << "kappa\tTotalE\tIntE\tGegenE\tEBFE\tell\tell_ebf\tQsize" << endl;

  for( double kappa = 1; kappa < 10000; kappa *= 1.20 )
  //double kappa = 40;
  {
    HelmKernel K;
    K.kappa = kappa;

    // Construct a Quadrature
    Quadrature* quad = new Quadrature(K, boxSize, eps);

    // Get the EBF L
    int EBF_L = get_Truncature(kappa, norm_r, norm_r0, eps, 2);

    // Construct a Transfer Function
    Transfer_Function T(quad, r0);

    // Construct Scratch Space
    NFunction A(quad);

    double totE_k = 0, intE_k = 0, gegE_k = 0, ebfE_k = 0;
    Vec3 rMax;

    // For all the r-directions
    for( int k = 0; k < (int) RCUBE.size(); ++k ) {

      Vec3 rhat = RCUBE[k] / RCUBE[k].mag();
      Vec3 r = rhat * norm_r;

      // A = T * G(r)
      Translation_Function::times( T, r, A );

      // Integrate
      complex I = A.integrate();

      // Get exact value
      complex Ie = K( (r+r0).mag() );

      // Gegenbauer value
      complex G = gegenbauer_series(quad->getTruncation(), kappa, r, r0);

      // Gegenbauer value with EBF
      complex G_ebf = gegenbauer_series(EBF_L, kappa, r, r0);

      //cout << kappa << "\t" << abs(Ie-I) << "\t" << abs(G-I) << "\t" << abs(Ie-G) << endl;

      // Total error
      totE_k = max( totE_k, abs( Ie - I ) );
      // Integration error
      if( abs( G - I ) > intE_k ) {
	intE_k = abs(G-I);
	rMax = r;
      }
      intE_k = max( intE_k, abs( G - I ) );
      // Gegenbauer Error
      gegE_k = max( gegE_k, abs( Ie - G ) );
      // EBF Gegenbauer Error
      ebfE_k = max( ebfE_k, abs( Ie - G_ebf ) );
    }

    myfile << kappa << "\t"
	   << totE_k << "\t" << intE_k << "\t"
	   << gegE_k << "\t" << ebfE_k << "\t"
	   << quad->getTruncation() << "\t" << EBF_L << "\t"
	   << quad->size() << endl;

    cout << kappa << "\t"
	 << totE_k << "\t" << intE_k << "\t"
	 << gegE_k << "\t" << ebfE_k << "\t"
	 << rMax << "\t" << quad->getTruncation() << "\t" << EBF_L << "\t"
	 << quad->size() << endl;

    delete quad;
  }

  myfile.close();

  cout << timer.stop() << endl;
  return 0;
}

#endif

#endif
