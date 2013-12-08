#ifndef DIRECT_H
#define DIRECT_H

template <typename FUNC>
void Direct(FUNC K, const vector<Vec3>& p,
	    const vector<complex>& psi, vector<complex>& omega)
{
  int N = p.size();
  for( int m = 0; m < N; ++m ) {
    const Vec3& pm = p[m];
    omega[m] = 0;
    for( int n = 0; n != m; ++n ) {
      double r = pm.dist(p[n]);
      complex E = K(r);
      omega[m] += psi[n] * E;
      omega[n] += psi[m] * E;
    }
  }
}

#endif
