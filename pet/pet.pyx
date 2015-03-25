
import numpy as np
cimport numpy as np
# from Cython.Compiler.Naming import retval_cname
from cpython cimport PyObject, Py_INCREF
# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
# cimport cython

np.import_array()

cdef extern from 'lib/vector.hpp': # namespace 'pet':
  cdef cppclass CVector "Vector":
    float& operator[](size_t idx)

cdef extern from 'lib/array3d.hpp': # namespace 'pet':
  cdef cppclass CArray3D "Array3D"[T]:
    CArray3D()
    T* get_data_raw()
    size_t get_size_ni() const
    size_t get_size_nj() const
    size_t get_size_nk() const

#  ctypedef CArray3D[float] CTArray3D

cdef extern from 'lib/patch.hpp':
  cdef cppclass CPatch "Patch":
    CPatch()
    CArray3D[CVector]* get_grid_E()
    CArray3D[CVector]* get_grid_B()
    CArray3D[CVector]* get_grid_DE(size_t idx)
    CArray3D[CVector]* get_grid_DB(size_t idx)
    CArray3D[CVector]* get_grid_AE(size_t idx)
    CArray3D[CVector]* get_grid_AB(size_t idx)
    void calc_grid_B()


cdef class Array3D:
  cdef CArray3D[float]* _this
  def __cinit__(self):
    pass

  # def __dealloc__(self):
  #   """ Frees the array. This is called by Python when all the
  #       references to the object are gone. """
  #   free(<void*>self.data_ptr)

  def __array__(self):
    """ Here we use the __array__ method, that is called when numpy
        tries to get an array from the object."""
    cdef np.npy_intp shape[3]
    shape[0] = <np.npy_intp> self._this.get_size_ni()
    shape[1] = <np.npy_intp> self._this.get_size_nj()
    shape[2] = <np.npy_intp> self._this.get_size_nk()
    # Create a 1D array, of length 'size'
    ndarray = np.PyArray_SimpleNewFromData(3, shape, np.NPY_FLOAT,
                                           self._this.get_data_raw())
    return ndarray

  def asarray(self, copy=False):
    cdef np.ndarray ndarray
    ndarray = np.array(self, copy=copy)
    ndarray.base = <PyObject*> self
    Py_INCREF(self)
    return ndarray



cdef class VectorField:
  cdef CArray3D[CVector]* _this
  def __cinit__(self):
    pass

  def __array__(self):
    cdef np.npy_intp shape[4]
    shape[0] = <np.npy_intp> self._this.get_size_ni()
    shape[1] = <np.npy_intp> self._this.get_size_nj()
    shape[2] = <np.npy_intp> self._this.get_size_nk()
    shape[3] = <np.npy_intp> 3
    ndarray = np.PyArray_SimpleNewFromData(4, shape, np.NPY_FLOAT,
                                           self._this.get_data_raw())
    return ndarray

  def asarray(self, copy=False):
    cdef np.ndarray ndarray
    ndarray = np.array(self, copy=copy)
    ndarray.base = <PyObject*> self
    Py_INCREF(self)
    return ndarray


cdef class Patch:
  cdef CPatch *_this
  def __cinit__(self):
    self._this = new CPatch()
  def __dealloc__(self):
    del self._this

  def calc_grid_B(self):
    self._this.calc_grid_B()

  def get_grid_E(self):
    array = VectorField()
    array._this = self._this.get_grid_E()
    return array

  def get_grid_B(self):
    array = VectorField()
    array._this = self._this.get_grid_B()
    return array

  def get_grid_DE(self, idx):
    array = VectorField()
    array._this = self._this.get_grid_DE(idx)
    return array

  def get_grid_DB(self, idx):
    array = VectorField()
    array._this = self._this.get_grid_DB(idx)
    return array

  def get_grid_AE(self, idx):
    array = VectorField()
    array._this = self._this.get_grid_AE(idx)
    return array

  def get_grid_AB(self, idx):
    array = VectorField()
    array._this = self._this.get_grid_AB(idx)
    return array

  property QE:
    def __get__(self):
      return self.get_grid_E().asarray()

  property QB:
    def __get__(self):
      return self.get_grid_B().asarray()

  property DE0:
    def __get__(self):
      return self.get_grid_DE(0).asarray()
  property DE1:
    def __get__(self):
      return self.get_grid_DE(1).asarray()
  property DE2:
    def __get__(self):
      return self.get_grid_DE(2).asarray()

  property DB0:
    def __get__(self):
      return self.get_grid_DB(0).asarray()
  property DB1:
    def __get__(self):
      return self.get_grid_DB(1).asarray()
  property DB2:
    def __get__(self):
      return self.get_grid_DB(2).asarray()


  property AE0:
    def __get__(self):
      return self.get_grid_AE(0).asarray()
  property AE1:
    def __get__(self):
      return self.get_grid_AE(1).asarray()
  property AE2:
    def __get__(self):
      return self.get_grid_AE(2).asarray()

  property AB0:
    def __get__(self):
      return self.get_grid_AB(0).asarray()
  property AB1:
    def __get__(self):
      return self.get_grid_AB(1).asarray()
  property AB2:
    def __get__(self):
      return self.get_grid_AB(2).asarray()
