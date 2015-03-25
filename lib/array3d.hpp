
#ifndef __PET_BASE_ARRAY3D_H__
#define __PET_BASE_ARRAY3D_H__

#include <cstddef>
#include <cstring>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#define DBG_ASSERT_BOUNDS(i,j,k) assert(i>=0 && j>=0 && k>=0 && i<_ni && j<_nj && k<_nk)

typedef std::size_t idx_t;

template<typename TYPE>
class Array3D
{
public:
  Array3D() {}
  Array3D(size_t ni, size_t nj, size_t nk)
  { initialize(ni, nj, nk); }

  ~Array3D()
  { free_mem(); }

  void initialize(size_t ni, size_t nj, size_t nk)
  { alloc_mem(ni, nj, nk); }

  TYPE& operator()(idx_t i, idx_t j, idx_t k)
  {
    DBG_ASSERT_BOUNDS(i,j,k);
    return _data[(i*_nj + j)*_nk + k];
  }

  const TYPE& operator()(idx_t i, idx_t j, idx_t k) const
  {
    DBG_ASSERT_BOUNDS(i,j,k);
    return _data[(i*_nj + j)*_nk + k];
  }

  size_t get_size_ni() const
  { return _ni; }
  size_t get_size_nj() const
  { return _nj; }
  size_t get_size_nk() const
  { return _nk; }

  TYPE* get_data_raw()
  { return _data; }

private:
  TYPE* _data = NULL;
  size_t _ni, _nj, _nk;

  //*************************************************************//

private:
  void free_mem()
  {
    if (_data)
      delete[] _data;
    _data = NULL;
    _ni = 0; _nj = 0; _nk = 0;
  }

  void alloc_mem(size_t ni, size_t nj, size_t nk)
  {
    free_mem();
    _data = new TYPE[ni*nj*nk];
    std::memset(_data, 0, ni*nj*nk * sizeof(TYPE));
    _ni = ni; _nj = nj; _nk=nk;
  }

};

#endif // __PET_BASE_ARRAY3D_H__

