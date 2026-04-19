#pragma once

#include <cmath>
#include <ostream>

struct Pos2D {
  double x, y;
  
  Pos2D& operator+=(const Pos2D& other) {
    x += other.x;
    y += other.y;
    return *this;
  }
  
  friend std::ostream &operator<<(std::ostream &s, const Pos2D &pos) {
    s << "{ " << pos.x << ", " << pos.y << " }";
    return s;
  }
};

