
#ifndef __PET_BASE_PATCH_H__
#define __PET_BASE_PATCH_H__

#include "array3d.hpp"
#include "vector.hpp"

/// Number of cells in a patch (NC{I,J,K} >= 1)
#define NCI 10
#define NCJ 2
#define NCK 1

class Patch
{
public:
  typedef Array3D<float> TArray3D;
  typedef Array3D<Vector> TVectorField;
  //  typedef Array3D<Vector> TTensorField;

  Patch() :
    _QE(NCI+1, NCJ+1, NCK+1),
    _QB(NCI  , NCJ  , NCK  )
  {
    for(auto i=0; i<3; ++i)
    {
      size_t ipe = (i == 0) ? 1 : 0;
      size_t jpe = (i == 1) ? 1 : 0;
      size_t kpe = (i == 2) ? 1 : 0;

      _DE[i].initialize(NCI+1, NCJ+1, NCK+1);
      _DB[i].initialize(NCI+1, NCJ+1, NCK+1);

      _AE[i].initialize(NCI  , NCJ  , NCK  );
      _AB[i].initialize(NCI+ipe, NCJ+jpe, NCK+kpe);
    }
  };
  ~Patch() {};

  /**
   * This function interpolates B grid points $r_B$ from the E grid.
   * @note We assume that E grid $r_E$ is already set up.
   */
  void calc_grid_B();

  TVectorField* get_grid_E()
  { return &_QE; }
  TVectorField* get_grid_B()
  { return &_QB; }

  TVectorField* get_grid_DE(size_t idx)
  { return &_DE[idx]; }
  TVectorField* get_grid_DB(size_t idx)
  { return &_DB[idx]; }

  TVectorField* get_grid_AE(size_t idx)
  { return &_AE[idx]; }
  TVectorField* get_grid_AB(size_t idx)
  { return &_AB[idx]; }

private:
  TVectorField _QE, _QB; //< Grid points positions (respective center of div operation
  TVectorField _DE[3], _DB[3]; //< Grid tangent vectors
  TVectorField _AE[3], _AB[3]; //< Grid tangent vectors
};

#endif // __PET_BASE_PATCH_H__
