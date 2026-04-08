#pragma once

#include <cmath>
#include <ostream>

struct Pos1D {
  double x;  // координата X
  
  Pos1D& operator+=(double dx) {
    x += dx;
    return *this;
  }
  
  friend std::ostream &operator<<(std::ostream &s, const Pos1D &pos) {
    s << "{ " << pos.x << " }";
    return s;
  }
};