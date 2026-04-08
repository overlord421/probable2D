#pragma once

#include <cmath>
#include <ostream>

struct Vec2 {
  double x, y;
  double dot(const Vec2 &other) const { return x * other.x + y * other.y; }
  double length() const { return std::sqrt(this->dot(*this)); }
  // В 2D cross product возвращает скаляр (псевдоскаляр)
  Vec2 cross_with_B(double Bz) const { return {y * Bz, - x * Bz}; }
  
  Vec2 sqrt() const { return {std::sqrt(x), std::sqrt(y)}; }
  
  Vec2 &operator+=(const Vec2 &other) {
    x += other.x;
    y += other.y;
    return *this;
  }
  friend Vec2 operator+(Vec2 a, const Vec2 &b) { return a += b; }
  
  Vec2 &operator-=(const Vec2 &other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }
  friend Vec2 operator-(Vec2 a, const Vec2 &b) { return a -= b; }
  
  Vec2 &operator*=(const Vec2 &other) {
    x *= other.x;
    y *= other.y;
    return *this;
  }
  friend Vec2 operator*(Vec2 a, const Vec2 &b) { return a *= b; }
  
  Vec2 &operator*=(double other) {
    x *= other;
    y *= other;
    return *this;
  }
  friend Vec2 operator*(Vec2 a, double b) { return a *= b; }
  friend Vec2 operator*(double b, Vec2 a) { return a *= b; }
  
  Vec2 &operator/=(double other) {
    x /= other;
    y /= other;
    return *this;
  }
  friend Vec2 operator/(Vec2 a, double b) { return a /= b; }
  
  friend std::ostream &operator<<(std::ostream &s, const Vec2 &v) {
    s << "{ " << v.x << ", " << v.y << " }";
    return s;
  }
};