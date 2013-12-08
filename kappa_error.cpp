const double MLFMM_ALPHA(1.0);

#include "General.hpp"

#include "MLFMM_Env.hpp"
#include "Direct.hpp"


struct HelmKernel
{
  double kappa;
  inline complex operator()(double r) {
    return exp(CI * kappa * r) / r;
  }
};


// Error as a function of kappa
// Demonstrate accuracy of direct L and quadrature generation
int main( int argc, char* argv[] )
{
  (void) argc; (void) argv;

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
