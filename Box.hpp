#ifndef BOX_H
#define BOX_H

#include "General.hpp"
#include "Vec3.hpp"
#include "NFunction.hpp"

class Box
{
  NFunction M;             // Multipole (outgoing) expansion
  NFunction L;             // Local (incoming) expansion

 public:
  typedef vector<int> PointList;
  typedef PointList::iterator pointIter;

  int n;                   // The box number
  Vec3 center;             // The box center

  Box* parent;             // A pointer to this box's parent
  PointList pointIndex;    // If leaf, contains indices of points it contains

  Box() {}
  Box( int n_, const Vec3& center_ ) : n(n_), center(center_) {}
  ~Box() {}

  inline void addPointIndex( int index ) {
    pointIndex.push_back( index );
  }

  inline void makeMultipole( Quadrature* q ) {
    M = NFunction(q);
  }

  inline void makeLocal( Quadrature* q ) {
    L = NFunction(q);
  }

  inline NFunction& getMultipole() {
    return M;
  }

  inline NFunction& getLocal() {
    return L;
  }

 private:
  // Disable Copy and Assignment
  Box(const Box& S) { (void) S; }
  void operator=(const Box& S) { (void) S; }
};

#endif
