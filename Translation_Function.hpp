#ifndef TRANSLATION_FUNCTION_H
#define TRANSLATION_FUNCTION_H

#include "Quadrature.hpp"
#include "NFunction.hpp"

class Translation_Function : public NFunction
{
 public:

  // Precomputes the Translation function for vector r
 Translation_Function(Quadrature* quad, const Vec3& r)
   : NFunction(quad)
  {
    //cout << "r = " << r << endl;

    Quadrature& q = *quad;

    double kappa = q.kappa;
    int K = q.size();
    for( int k = 0; k < K; ++k ) {
      const Vec3& s = q.getPoint(k).s;
      C[k] = exp( CI * kappa * r.dot(s) );
    }
  }

  // B += psi * A(r)
  inline void static add( const Vec3& r, const complex psi, NFunction& B )
  {
    //cout << "r = " << r << endl;

    Quadrature& quad = B.getQuad();

    double kappa = quad.kappa;
    int K = quad.size();
    for( int k = 0; k < K; ++k ) {
      const double krs = kappa * r.dot( quad.getPoint(k).s );
      B.C[k] += psi * complex( cos(krs), sin(krs) );
    }
  }

  // B = C .* A(r)
  inline void static times( const NFunction& C, const Vec3& r, NFunction& B )
  {
    //cout << "r = " << r << endl;

    Quadrature& quad = B.getQuad();

    double kappa = quad.kappa;
    int K = quad.size();
    for( int k = 0; k < K; ++k ) {
      const double krs = kappa * r.dot( quad.getPoint(k).s );
      B.C[k] = C.C[k] * complex( cos(krs), sin(krs) );
    }
  }

};

#endif
