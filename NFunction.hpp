#ifndef NFUNCTION_H
#define NFUNCTION_H

#include "General.hpp"
#include "Quadrature.hpp"

// A class to define a numerical function over a Quadrature

class NFunction
{
 public:
  Quadrature* quad;     // Pointer to quadrature of field values
  vector<complex> C;    // The field values for each point of the quadrature

  // Construct a numerical function over a given quadrature
  NFunction( Quadrature* q_ = NULL ) { setQuad(q_); }
  // Destructor
  ~NFunction() {
    // Don't delete quad, it doesn't belong to us
  }

  inline void setQuad( Quadrature* q )
  {
    quad = q;
    C = vector<complex>( quad == NULL ? 0 : quad->size() );
  }

  inline Quadrature& getQuad()
  {
    return *quad;
  }

  // This = 0
  inline void zero()
  {
    C.assign( C.size(), 0 );
  }

  // This += A
  inline void add( const NFunction& A )
  {
    assert( quad == A.quad );

    int K = C.size();
    //const vector<complex>& AC = A.C;

    for( int k = 0; k < K; ++k )
      C[k] += A.C[k];
  }

  // This = A*B
  inline void setProduct( const NFunction& A, const NFunction& B )
  {
    assert( quad == A.quad );
    assert( A.quad == B.quad );

    int K = C.size();
    for( int k = 0; k < K; ++k )
      C[k] = A.C[k] * B.C[k];
  }

  // This += A*B
  inline void addProduct( const NFunction& A, const NFunction& B )
  {
    assert( quad == A.quad );
    assert( A.quad == B.quad );

    int K = C.size();
    for( int k = 0; k < K; ++k )
      C[k] += A.C[k] * B.C[k];
  }

  inline complex integrate()
  {
    int K = C.size();
    Quadrature& q = *quad;

    complex total = 0;
    for( int k = 0; k < K; ++k ) {
      total += q.getPoint(k).w * C[k];
      std::cout << k << ":\t" << q.getPoint(k).w << std::endl;
    }

    return total;
  }
  friend ostream& operator<<(ostream& os, const NFunction& A) {
    return os << A.C;
  }
};

#endif
