#ifndef TRANSFER_VECTOR_H
#define TRANSFER_VECTOR_H

#include "General.hpp"
#include "Box.hpp"
#include "Transfer_Function.hpp"

struct Close_Vector 
{
  Box* b1;
  Box* b2;
  
Close_Vector( Box* b1_, Box* b2_ ) : b1(b1_), b2(b2_) {}

  friend ostream& operator<<(ostream& os, const Close_Vector& t)
  {
    return (os << t.b1->n << "-" << t.b2->n);
  }
};


struct Trans_Vector 
{
  Box* b1;                  // Pointer to box 1 (from box)
  Box* b2;                  // Pointer to box 2 (to box)

  Transfer_Function* T;     // Pointer to Transfer Function to be applied

  int x, y, z;              // Representation of r_0 in number of boxes

  Trans_Vector() : b1(NULL), b2(NULL), T(NULL), x(0), y(0), z(0) {}
  
  Trans_Vector( Box* b1_, Box* b2_, Vec3 r0_ )
  : b1(b1_), b2(b2_), T(NULL), x(r0_.x), y(r0_.y), z(r0_.z) {}

  // Comparator for sorting
  inline bool operator<(const Trans_Vector& t) const {
    return (x < t.x 
	|| (x == t.x && (y < t.y
		     || (y == t.y && z < t.z))));
  }

  friend ostream& operator<<(ostream& os, const Trans_Vector& t)
  {
    return (os << t.b2->n << "->" << t.b1->n << ":"
	    << "\t(" << t.x << ", " << t.y << ", " << t.z << ")");
  }
};

#endif

