
#ifndef __PET_BASE_VECTOR_H__
#define __PET_BASE_VECTOR_H__

class Vector
{
public:
  Vector() {};
  Vector(float x, float y, float z)
  { _data[0] = x; _data[1] = y; _data[2] = z; };
  ~Vector() {};

  float& operator[](size_t idx)
  { return _data[idx]; }

  float operator[](size_t idx) const
  { return _data[idx]; }


  Vector& operator+=(const Vector& rhs)
  {
    for(size_t i=0; i<3; ++i)
      _data[i] += rhs[i];
    return *this;
  }

  Vector& operator-=(const Vector& rhs)
  {
    for(size_t i=0; i<3; ++i)
      _data[i] -= rhs[i];
    return *this;
  }

  Vector& operator*=(const float& rhs) // compound assignment (does not need to be a member,
  {                           // but often is, to modify the private members)
    for(size_t i=0; i<3; ++i)
      _data[i] *= rhs;
    return *this; // return the result by reference
  }

  /**
   * @Note Passing first arg by value helps optimize chained a+b+c
   * alternatively, both parameters may be const references.
   */
  friend Vector operator+(Vector lhs, const Vector& rhs)
  { return lhs += rhs; } // reuse compound assignment and return the result by value

  friend Vector operator-(Vector lhs, const Vector& rhs)
  { return lhs -= rhs; } // reuse compound assignment and return the result by value

  friend Vector operator*(Vector lhs, const float& rhs)
  { return lhs *= rhs; } // reuse compound assignment and return the result by value

  friend Vector operator*(float lhs, const Vector& rhs)
  {
    Vector tmp;
    for(size_t i=0; i<3; ++i)
      tmp[i] = lhs * rhs[i];
    return tmp;
  } // reuse compound assignment and return the result by value

  friend Vector operator% (const Vector &v1, const Vector &v2)
  {
    return Vector(v1._data[1] * v2._data[2] - v1._data[2] * v2._data[1],
                  v1._data[2] * v2._data[0] - v1._data[0] * v2._data[2],
                  v1._data[0] * v2._data[1] - v1._data[1] * v2._data[0]);
  }

private:
  float _data[3];
};

#endif // __PET_BASE_VECTOR_H__
