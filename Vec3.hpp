#ifndef VEC3_H
#define VEC3_H

struct Vec3
{
  double x, y, z;

  // Constructor
 Vec3(double x_ = 0, double y_ = 0, double z_ = 0)
     : x(x_), y(y_), z(z_) {
 }

  /**********************/
  // Operator Overloads //
  /**********************/

  inline bool operator==(const Vec3& v) const {
    return (x == v.x && y == v.y && z == v.z);
  }

  inline Vec3 operator+(const Vec3& v) const {
    return Vec3(x + v.x, y + v.y, z + v.z);
  }

  inline Vec3 operator+(double S) const {
    return Vec3(x + S, y + S, z + S);
  }

  inline void operator+=(const Vec3& v) {
    x += v.x; y += v.y; z += v.z;
  }

  inline void operator+=(double S) {
    x += S; y += S; z += S;
  }

  inline Vec3 operator-() const {
    return Vec3(-x, -y, -z);
  }

  inline Vec3 operator-(const Vec3& v) const {
    return Vec3(x - v.x,  y - v.y,  z - v.z);
  }

  inline Vec3 operator-(double S) const {
    return Vec3(x - S, y - S, z - S);
  }

  inline void operator-=(const Vec3& v) {
    x -= v.x; y -= v.y; z -= v.z;
  }

  inline void operator-=(double S) {
    x -= S; y -= S; z -= S;
  }

  inline Vec3 operator/(const Vec3& v) const {
    return Vec3(x / v.x,  y / v.y,  z / v.z);
  }

  inline Vec3 operator/(double S) const {
    return Vec3(x / S, y / S, z / S);
  }

  inline void operator/=(const Vec3& v) {
    x /= v.x; y /= v.y; z /= v.z;
  }

  inline void operator/=(double S) {
    x /= S; y /= S; z /= S;
  }

  inline Vec3 operator*(const Vec3& v) const {
    return Vec3(x * v.x,  y * v.y,  z * v.z);
  }

  inline Vec3 operator*(double S) const {
    return Vec3(x * S,  y * S,  z * S);
  }

  inline void operator*=(const Vec3& v) {
    x *= v.x; y *= v.y; z *= v.z;
  }

  inline void operator*=(double S) {
    x *= S; y *= S; z *= S;
  }

  /*************/
  // Functions //
  /*************/

  inline double dot(const Vec3 &v) const {
    return x*v.x + y*v.y + z*v.z;
  }

  inline Vec3 cross(const Vec3 &v) const {
    return Vec3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
  }

  inline double magSq() const {
    return x*x + y*y + z*z;
  }

  inline double mag() const {
    return sqrt(this->magSq());
  }

  inline double dist(const Vec3& v) const {
    return (v - *this).mag();
  }

  inline void normalize() {
    double mag = this->mag();
    if(mag != 0) (*this) *= 1.0/mag;
  }

  /**********/
  // Output //
  /**********/

  friend ostream& operator<<(ostream& os, const Vec3& v) {
    return (os << "(" << v.x << "," << v.y << "," << v.z << ")");
  }
};

#endif
