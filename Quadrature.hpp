#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "General.hpp"
#include "Vec3.hpp"
#include "HelmholtzUtil.hpp"

struct QuadraturePoint
{
  int index;
  double theta, phi, w;
  Vec3 s;
  QuadraturePoint( double theta_ = 0, double phi_ = 0, double w_ = 0 )
      : theta(theta_), phi(phi_), w(w_),
        s(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)) {}
  ~QuadraturePoint() {}
};


// A row of quadrature points on the sphere
struct QuadratureRow
{
  int index;
  double theta;
  vector<QuadraturePoint> point;

  // A row of N_phi quadrature points (of total weight w0) at constant theta
  QuadratureRow( double theta_ = 0, int N_phi = 0, double w0 = 0 )
      : theta(theta_), point(N_phi)
  {
    // Construct the quadrature points
    for( int n = 0; n < N_phi; ++n ) {
      double phi = n * 2*PI/N_phi;
      double w = w0 * 2*PI/N_phi;
      point[n] = QuadraturePoint(theta, phi, w);
    }
  }
  ~QuadratureRow() {}

  inline int size() const { return point.size(); }
  inline int numPoints() const { return size(); }
  inline double getTheta() const { return theta; }
  inline QuadraturePoint& getPoint( int k ) { return point[k]; }

  // Output
  friend ostream& operator<<(ostream& os, const QuadratureRow& quad) {
    (void) quad;
    return os;
  }
};


class QuadratureHelm
{
 public:
  const double kappa;

  vector<QuadratureRow> row;         // List of quadrature rows (const theta)
  vector<QuadraturePoint*> point;    // Enumerated quadrature points

  int L;                             // Gegenbauer truncation
  int max_N_phi;

  // Index maps for rotations and reflections of the quadrature
  vector< vector<int> > indexMap;

  template <typename FUNC>
  QuadratureHelm( FUNC K ) : kappa(K.kappa) {}
  ~QuadratureHelm() {}

  inline int size() const { return point.size(); }
  inline int numPoints() const { return size(); }
  inline int numRows() const { return row.size(); }
  inline int maxRowSize() const { return max_N_phi; }
  inline int getTruncation() const { return L; }
  inline QuadratureRow& getRow( int k ) { return row[k]; }
  inline const QuadratureRow& getRow( int k ) const { return row[k]; }
  inline QuadraturePoint& getPoint( int k ) { return *(point[k]); }

  // Output
  friend ostream& operator<<(ostream& os, const QuadratureHelm& q)
  {
    os << "(" << q.getTruncation()
       << "," << q.numRows()
       << "," << q.size() << ")" << endl;
    for (int k = 0; k < q.numRows(); ++k) {
      os << "\t" << q.getRow(k).size() << endl;
    }
    return os;
  }

 protected:

  inline void makePointerAccess()
  {
    // Quadrature is done being constructed
    // Construct the enumerated point vector
    int rowNumber = 0, pointNumber = 0;
    max_N_phi = 0;
    for( int k = 0; k < numRows(); ++k ) {
      // Assign each a row number
      row[k].index = rowNumber++;
      // Find the largest row
      max_N_phi = max( max_N_phi, row[k].size() );

      for( int p = 0; p < row[k].numPoints(); ++p ) {
        // Assign each a point number
        row[k].point[p].index = pointNumber++;
        // Insert point
        point.push_back( &(row[k].point[p]) );
      }
    }
    //cout << pointNumber << endl;
  }
};


// A uniform quadrature in phi and theta
class Quadrature_Uniform : public QuadratureHelm
{
 public:

  template <typename FUNC>
  Quadrature_Uniform( FUNC K, double boxSize, double eps )
      : QuadratureHelm(K)
  {
    // Characteristic lengths (worse case)
    double norm_r0 = 2 * boxSize;
    double norm_r = MLFMM_ALPHA * sqrt(3) * boxSize;

    /* Determine the Gegenbauer truncation (L) */
    L = get_Truncature(kappa, norm_r, norm_r0, eps);
    int maxL = L + 20;

    cerr << "norm_r0 = " << norm_r0 << endl;
    cerr << "norm_r = " << norm_r << endl;
    cerr << "L = " << L << endl;

    /* Determine the number of quadrature rows (N_rows = N_theta/2 + 1) */
    int N_theta = get_N_theta(L, kappa, norm_r, norm_r0, eps);
    int N_rows = N_theta/2 + 1;
    int N_phiT = 2*L + 2;
    row.resize( N_rows );

    //cerr << "N_theta = " << N_theta << endl;

    Vec3 r0(norm_r0, 0, 0);
    vector< vector<complex> > TL = mod_transfer_fn(L,kappa,N_rows,N_phiT,r0);

    // Prune N_phi(theta) if possible
    for( int n = 0; n < N_rows; ++n ) {
      double theta = n * PI/(N_rows-1);

      if( n == 0 || n == N_rows-1 ) {   // The Poles

        double w0 = 2*PI/N_theta;
        int N_phi = 1;
        row[n] = QuadratureRow( theta, N_phi, w0 );

      } else if( n <= N_rows/2 ) {      // Positive z row

        // The weight of this row, 2x for symmetry
        double w0 = 2 * 2*PI/N_theta;

        // Determine N_phi for this row
        double knormrs = kappa * norm_r * sin( theta );

        // EBF METHOD (This is pretty good...)
        /*
          int maxL_ebf = ceil(knormrs + 1.8*pow(-log10(eps), 2.0/3.0)
          *pow(knormrs, 1.0/3.0));
          int N_phi_ebf = round(2*maxL_ebf + 1, 4);
          int N_phi = N_phi_ebf;
        */
        // Another Heuristic
        //int N_phi = round(2*L+1,4);

        // Compute the transfer function coefficients
        ifft( &TL[n][0], N_phiT );

        vector<double> absT(L+1);
        for( int m = 0; m <= L; ++m ) {
          absT[m] = abs( TL[n][m] );
          // Make sure this is correct and symmetric
          //assert
        }

        // Precompute Bessel J
        vector<double> absJ( ceil(2*maxL+1, 4) + 1 );
        for( int k = 0; k < (int) absJ.size(); ++k ) {
          absJ[k] = abs( bessel_J(k, knormrs) );
        }

        // Compute N_phi by searching for error
        int N_phi_max = ceil(2*maxL+8, 4);
        int N_phi;
        for( N_phi = 4; N_phi <= N_phi_max; N_phi += 4 ) {
          //for( N_phi = ceil(2*maxL+1, 4); N_phi > 4; N_phi -= 4 ) {

          // Get the truncation error sum_{|k| >= N_phi/2}
          double error_trunc = 0;
          for( int m = N_phi/2; m <= L; ++m ) {
            error_trunc += absT[m] * absJ[m];
          }

          // Get the aliasing error
          double error_alias = absT[0] * absJ[N_phi];
          for( int m = 1; m <= N_phi/2-1; ++m ) {
            error_alias += absT[m] * absJ[N_phi - m];
          }

          double error = 4*PI*PI * 2*( error_alias + error_trunc );
          //cerr << N_phi << ": " << error << "\t" << error_alias << "\t" << error_trunc << endl;

          if( error < eps ) break;
          //if( error > eps ) { N_phi += 4; break; }
        }

        if( N_phi >= N_phi_max ) {
          cerr << "WARNING: N_phi ceiling hit" << endl;
        }

        //cout << N_phi << ", " << N_phi_ebf << endl;
        row[n] = QuadratureRow( theta, N_phi, w0 );

      } else {                         // Negative z row

        // The weight of this row, 2x for symmetry
        double w0 = 2 * 2*PI/N_theta;
        int N_phi = row[N_rows-1-n].size();
        row[n] = QuadratureRow( theta, N_phi, w0 );

      }
    }

    makePointerAccess();
  }
  ~Quadrature_Uniform() {}

  inline int get_N_theta( int L, double kappa,
                          double norm_r, double norm_r0,
                          double eps )
  {
    int N_rows_max = L + 30;
    int N_theta_max = 2*N_rows_max - 2;

    // EBF METHOD
    /*
      int maxL_ebf = ell_ebf(kappa, norm_r, eps);
      int N_theta_ebf = 2*maxL_ebf + 2;
      cerr << "N_theta_ebf = " << N_theta_ebf << endl;
      return N_theta_ebf;
    */

    // Another Heuristic
    //return 2 * L + 2;

    // Get the transfer matrix for phi = 0, theta = [0,2PI)
    Vec3 r0(0,0,norm_r0);
    vector< vector<complex> > T = mod_transfer_fn(L,kappa,N_rows_max,2,r0);

    // Copy real-space values into a single vector t(theta)
    vector<complex> t(N_theta_max);
    t[0] = T[0][0];                        // North Pole
    for( int k = 1; k < N_rows_max; ++k ) {
      t[k] = T[k][0];
      t[N_theta_max-k] = T[k][1];
    }

    // Get the Fourier coefficients
    ifft( &t[0], N_theta_max );

    // Get the absolute value of the Fourier coefficients
    vector<double> absT( N_theta_max/2-1 );
    for( int k = 0; k < N_theta_max/2-1; ++k ) {
      absT[k] = abs(t[k]);
      // Make sure this is correct and symmetric
      //cout << k << ": " << abs( abs(t[k])-abs(t[N_theta_max-k]) ) << endl;
      //assert( k == 0 || abs( abs(t[k])-abs(t[N_theta_max-k]) ) < 1e-12 );
    }

    // Compute the bessel functions
    vector<double> absJ( N_theta_max+1 );
    for( int k = 0; k <= N_theta_max; ++k ) {
      absJ[k] = abs( bessel_J(k, kappa*norm_r) );
      //std::cerr << "|J_" << k << "(" << kappa*norm_r << ")| = " << absJ[k] << std::endl;
    }

    // Search for N_theta
    for( int N_theta = max(2*(L-10),0); N_theta < N_theta_max; N_theta += 2 ) {

      // Get the truncation error sum_{n >= N_theta/2} + sum_{n <= -N_theta/2}
      double error_trunc = 0;
      for( int n = N_theta/2; n < N_theta_max/2-1; ++n )
        error_trunc += absT[n] * absJ[n];

      // Get the aliasing error
      double error_alias = absT[0] * absJ[N_theta];
      for( int n = 1; n <= N_theta/2-1; ++n )
        error_alias += absT[n] * absJ[N_theta - n];

      double error = 4*PI*PI * 2*( error_alias + error_trunc );
      //cerr << N_theta << ": " << error << endl;

      if( error < eps ) return N_theta;
    }

    cerr << "Error: N_theta convergence failure" << endl;
    exit(1);
  }

};


// A standard spherical harmonic spherical quadrature
// Uniform quadrature points in phi
// Gauss-Legendre quadrature points in z(theta)
class Quadrature_S2 : public QuadratureHelm
{
 public:

  template <typename FUNC>
  Quadrature_S2( FUNC K, double boxSize, double eps )
      : QuadratureHelm(K)
  {
    // Characteristic lengths (worse case)
    double norm_r0 = 2 * boxSize;
    double norm_r = MLFMM_ALPHA * sqrt(3) * boxSize;

    // Determine the Gegenbauer truncation (L)
    L = get_Truncature(kappa, norm_r, norm_r0, eps);

    // Generate Gauss-Legendre quadrature of size N_rows in z(theta)
    int N_rows = L + 1;
    int N_phi  = 2*L + 2;

    row.resize(N_rows);

    vector<double> z(N_rows);
    vector<double> w(N_rows);
    getGaussLegendreQuad(N_rows, z, w);

    for( int n = 0; n < N_rows; ++n ) {
      double theta = acos( z[N_rows-1-n] );
      double w0 = w[N_rows-1-n];
      row[n] = QuadratureRow( theta, N_phi, w0 );
    }

    makePointerAccess();
  }
  ~Quadrature_S2() {}
};


typedef Quadrature_Uniform Quadrature;
//typedef QuadratureHelm Quadrature;


#endif
