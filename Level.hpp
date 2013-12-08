#ifndef LEVEL_H
#define LEVEL_H

#include "Box.hpp"
#include "Transfer_Function.hpp"
#include "Transfer_Vector.hpp"
#include "Translation_Function.hpp"
#include "Reterpolator.hpp"

typedef list<Trans_Vector> TVList;
typedef TVList::iterator TVIter;

typedef list<Close_Vector> CVList;
typedef CVList::iterator CVIter;

typedef vector<Transfer_Function*> TFList;
typedef TFList::iterator TFIter;

typedef vector<Translation_Function*> TLList;
typedef TLList::iterator TLIter;

typedef list<Box*> BoxSet;
typedef BoxSet::iterator BoxIter;

template <int DIM>
class Level
{
  const static int BRANCH = 1 << DIM;   // The branching factor = 2^DIM

  double boxSize;          // Size of a box this levels contains

  BoxSet boxList;          // The list of nonempty boxes this level contains

  Quadrature* quad;        // Quadrature for numerical functions of this level
  NFunction scratchSpace;  // Scratch space of size quadrature for computation

  TVList transList;        // The list of transfer vectors
  CVList closeList;        // The list of close vectors (necessary?)

  TFList transferFunc;     // The list of transfer functions (necessary?)
  TLList translationFunc;  // The list of translation functions

  Interpolator* upReterp;
  Anterpolator* downReterp;

 public:

  // Constructors
  Level( double boxSize_ = 0 )
      : boxSize(boxSize_),
        quad(NULL), translationFunc(BRANCH),
        upReterp(NULL), downReterp(NULL) {}

  // Destructor
  ~Level()
  {
    // Delete Transfer Functions
    for( TFIter tfi = transferFunc.begin();
         tfi != transferFunc.end();
         ++tfi ) {
      Transfer_Function* T = *tfi;
      delete T;
    }
    // Delete Translation Functions
    for( TLIter tli = translationFunc.begin();
         tli != translationFunc.end();
         ++tli ) {
      Translation_Function* TL = *tli;
      delete TL;
    }
    // Delete Quadrature and Reterps
    delete quad;
    delete upReterp;
    delete downReterp;
  }

  template <typename FUNC>
  inline void defineQuadrature( FUNC K, double eps )
  {
    // Construct the Quadrature
    quad = new Quadrature( K, boxSize, eps );

    cerr << *quad << endl;

    // Construct some scratchspace for efficiency
    scratchSpace = NFunction( quad );

    // Initialize all the box multipole and locals with the quadrature
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ) {
      Box* b = *bi;
      b->makeMultipole( quad );
      b->makeLocal( quad );
    }
  }

  inline Quadrature* getQuad()
  {
    return quad;
  }

  inline NFunction& getScratch()
  {
    return scratchSpace;
  }

  inline void zeroFields()
  {
    // Zero the box multipole and locals
    for( BoxIter bi = boxbegin(); bi != boxend(); ++bi ) {
      Box* b = *bi;
      b->getMultipole().zero();
      b->getLocal().zero();
    }
  }

  inline void defineTransfers()
  {
    transList.sort();

    Trans_Vector LastTV;
    for( TVIter tvi = transList.begin(); tvi != transList.end(); ++tvi ) {
      Trans_Vector& tv = *tvi;

      if( LastTV.x != tv.x || LastTV.y != tv.y || LastTV.z != tv.z ) {
        // Compute a new Transfer Function
        Vec3 r0( tv.x * boxSize, tv.y * boxSize, tv.z * boxSize );
        std::cerr << "tv = " << tv << "    r0 = " << r0 << std::endl;
        tv.T = new Transfer_Function( quad, r0 );
      } else {
        // This Transfer Function is the same as the last one
        tv.T = LastTV.T;
      }

      LastTV = tv;
    }
  }

  // Define the translation functions from a lower level
  inline void defineTranslations()
  {
    // Depending on BRANCH, construct the translation vectors
    for( int k = 0; k < BRANCH; ++k ) {
      // Generate the kth translation vector
      // Could do this alot better with the NTree...
      Vec3 r;
      if( DIM == 1 ) {
        if( k == 0 ) r = Vec3( boxSize/4, 0, 0);
        if( k == 1 ) r = Vec3(-boxSize/4, 0, 0);
      }
      if( DIM == 2 ) {
        if( k == 0 ) r = Vec3( boxSize/4, boxSize/4, 0);
        if( k == 1 ) r = Vec3( boxSize/4,-boxSize/4, 0);
        if( k == 2 ) r = Vec3(-boxSize/4, boxSize/4, 0);
        if( k == 3 ) r = Vec3(-boxSize/4,-boxSize/4, 0);
      }
      if( DIM == 3 ) {
        if( k == 0 ) r = Vec3( boxSize/4, boxSize/4, boxSize/4);
        if( k == 1 ) r = Vec3( boxSize/4, boxSize/4,-boxSize/4);
        if( k == 2 ) r = Vec3( boxSize/4,-boxSize/4, boxSize/4);
        if( k == 3 ) r = Vec3( boxSize/4,-boxSize/4,-boxSize/4);
        if( k == 4 ) r = Vec3(-boxSize/4, boxSize/4, boxSize/4);
        if( k == 5 ) r = Vec3(-boxSize/4, boxSize/4,-boxSize/4);
        if( k == 6 ) r = Vec3(-boxSize/4,-boxSize/4, boxSize/4);
        if( k == 7 ) r = Vec3(-boxSize/4,-boxSize/4,-boxSize/4);
      }

      translationFunc[k] = new Translation_Function( quad, r );
    }
  }

  // Get the translation function from a child box b
  // to its parent box in this level
  inline Translation_Function& getTranslationUp( Box* b )
  {
    return *translationFunc[b->n & (BRANCH-1)];
  }

  // Get the translation function from a parent box to a child box b
  // in this level
  inline Translation_Function& getTranslationDown( Box* b )
  {
    return *translationFunc[(b->n & (BRANCH-1)) ^ (BRANCH-1)];
  }

  // Define an interpolator from the quadrature of this level
  // to the quadrature of level qB
  inline void defineInterp( Quadrature* qB )
  {
    upReterp = new Interpolator( quad, qB );
  }

  inline Interpolator& getInterp()
  {
    return *upReterp;
  }

  // Define an anterpolator from the quadrature of this level
  // to the quadrature of level qb
  inline void defineAnterp( Quadrature* qb )
  {
    downReterp = new Anterpolator( quad, qb );
  }

  inline Anterpolator& getAnterp()
  {
    return *downReterp;
  }

  inline double getBoxSize()
  {
    return boxSize;
  }

  inline void addBox( Box* b )
  {
    boxList.push_back(b);
  }

  inline void addTransfer( Box* b1, Box* b2 )
  {
    Vec3 r0 = b2->center - b1->center;
    r0 /= boxSize;
    r0.x = round(r0.x); r0.y = round(r0.y); r0.z = round(r0.z);
    //cerr << "Add Transfer " << b1->n << "->" << b2->n << ": " << r0 << endl;
    transList.push_back( Trans_Vector(b1, b2,  r0) );
    //cerr << "Add Transfer " << b2->n << "->" << b1->n << ": " << -r0 << endl;
    transList.push_back( Trans_Vector(b2, b1, -r0) );
  }

  inline void addClose( Box* b1, Box* b2 )
  {
    //cerr << "Add Close " << b1->n << "<->" << b2->n << endl;
    closeList.push_back( Close_Vector(b1, b2) );
  }

  inline void applyTransferFunctions()
  {
    for( TVIter tvi = transList.begin(); tvi != transList.end(); ++tvi ) {
      //cout << *tvi << endl;

      Transfer_Function* T = tvi->T;

      Box* b1 = tvi->b1;
      Box* b2 = tvi->b2;

      b2->getLocal().addProduct( *T, b1->getMultipole() );
    }
  }

  // Computes the direct product
  // omega[m] += sum_n  K(p[m] - p[n]) * psi[n]
  // for all pairs (m,n) that have been determined to be close
  template <typename FUNC>
  inline void applyClose( FUNC K,
                          const vector<complex>& psi,
                          const vector<Vec3>& p,
                          vector<complex>& omega )
  {
    for( CVIter cvi = closeList.begin(); cvi != closeList.end(); ++cvi ) {
      Box* b1 = cvi->b1;
      Box* b2 = cvi->b2;

      // For all pairs of points inside box b1 and b2
      const Box::pointIter nStart = b1->pointIndex.begin();
      const Box::pointIter nEnd   = b1->pointIndex.end();
      const Box::pointIter mStart = b2->pointIndex.begin();
      const Box::pointIter mEnd   = b2->pointIndex.end();

      for( Box::pointIter ni = nStart; ni != nEnd; ++ni ) {
        int n = *ni;
        const Vec3& pn = p[n];
        for( Box::pointIter mi = mStart; mi != mEnd; ++mi ) {
          int m = *mi;
          double r = pn.dist(p[m]);
          complex E = K(r);
          omega[m] += psi[n] * E;
          omega[n] += psi[m] * E;
        }
      }
    }
  }

  inline BoxIter boxbegin()
  {
    return boxList.begin();
  }

  inline BoxIter boxend()
  {
    return boxList.end();
  }
};

#endif
