
#include "patch.hpp"
#include <algorithm> // std::min/max

void Patch::calc_grid_B()
{
  /// Calculate B-cell co-variant (cell edge) vectors.
  size_t im, jm, km, ip, jp, kp;
  const size_t nei = _QE.get_size_ni(),
    nej = _QE.get_size_nj(), nek = _QE.get_size_nk();
  const size_t nbi = _QB.get_size_ni(),
    nbj = _QB.get_size_nj(), nbk = _QB.get_size_nk();

  /// Calculate B grid points.
  ///
  /// @Note: we will need to repeat this process couple of times to make B-grid
  /// and E-grid as orthogonal as possible.
  for(size_t i=0; i<nbi; ++i)
    for(size_t j=0; j<nbj; ++j)
      for(size_t k=0; k<nbk; ++k)
        _QB(i,j,k) = (float)0.125 * (
          _QE(i  ,j  ,k  ) + _QE(i  ,j  ,k+1) +
          _QE(i  ,j+1,k  ) + _QE(i  ,j+1,k+1) +
          _QE(i+1,j  ,k  ) + _QE(i+1,j  ,k+1) +
          _QE(i+1,j+1,k  ) + _QE(i+1,j+1,k+1) );


  /// Calculate B-cell co-variant (__cell edge__) vectors.
  TVectorField &DB0 = _DB[0], &DB1 = _DB[1], &DB2 = _DB[2];
  for(size_t ii=0; ii<nei; ++ii)
  {
    im = std::min(ii, nei-2); ip = std::min(ii+1, nei-1);
    for(size_t jj=0; jj<nej; ++jj)
    {
      jm = std::min(jj, nej-2); jp = std::min(jj+1, nej-1);
      for(size_t kk=0; kk<nek; ++kk)
      {
        km = std::min(kk, nek-2); kp = std::min(kk+1, nek-1);
        DB0(ii,jj,kk) = _QE(ip,jj,kk) - _QE(im,jj,kk);
        DB1(ii,jj,kk) = _QE(ii,jp,kk) - _QE(ii,jm,kk);
        DB2(ii,jj,kk) = _QE(ii,jj,kp) - _QE(ii,jj,km);
      }
    }
  }

  /// Calculate B-cell contra-variant (__cell face__) vectors.
  TVectorField &DE0 = _DE[0], &DE1 = _DE[1], &DE2 = _DE[2];
  for(size_t ii=0; ii<nbi; ++ii)
  {
    //im = ii > 0 ? ii-1 : 0; ip = ii > 0 ? ii : 1;
    ip = ii+1;
    for(size_t jj=0; jj<nbj; ++jj)
    {
      //jm = jj > 0 ? jj-1 : 0; jp = jj > 0 ? jj : 1;
      jp = jj+1;
      for(size_t kk=0; kk<nbk; ++kk)
      {
        //km = kk > 0 ? kk-1 : 0; kp = kk > 0 ? kk : 1;
        kp = kk+1;

        if(ii == 0)
          DE0(ii,jj,kk) = (float)2.0 * (_QB(ii,jj,kk) -
            (float)0.25 * (_QE(ii,jj,kk) + _QE(ii,jj,kp) + _QE(ii,jp,kk) + _QE(ii,jp,kp)));
        if(ii < nbi-1)
          DE0(ip,jj,kk) = _QB(ip,jj,kk) - _QB(ii,jj,kk);
        else
          DE0(ip,jj,kk) = (float)2.0 *
            ((float)0.25 * (_QE(ip,jj,kk) + _QE(ip,jj,kp) + _QE(ip,jp,kk) + _QE(ip,jp,kp)) -
             _QB(ii,jj,kk));

        if(jj == 0)
          DE1(ii,jj,kk) = (float)2.0 * (_QB(ii,jj,kk) -
            (float)0.25 * (_QE(ii,jj,kk) + _QE(ii,jj,kp) + _QE(ip,jj,kk) + _QE(ip,jj,kp)));
        if(jj < nbj-1)
          DE1(ii,jp,kk) = _QB(ii,jp,kk) - _QB(ii,jj,kk);
        else
          DE1(ii,jp,kk) = (float)2.0 *
            ((float)0.25 * (_QE(ii,jp,kk) + _QE(ii,jp,kp) + _QE(ip,jp,kk) + _QE(ip,jp,kp)) -
             _QB(ii,jj,kk));

        if(kk == 0)
          DE2(ii,jj,kk) = (float)2.0 * (_QB(ii,jj,kk) -
            (float)0.25 * (_QE(ii,jj,kk) + _QE(ii,jp,kk) + _QE(ip,jj,kk) + _QE(ip,jp,kk)));
        if(kk < nbk-1)
          DE2(ii,jj,kp) = _QB(ii,jj,kp) - _QB(ii,jj,kk);
        else
          DE2(ii,jj,kp) = (float)2.0 *
            ((float)0.25 * (_QE(ii,jj,kp) + _QE(ii,jp,kp) + _QE(ip,jj,kp) + _QE(ip,jp,kp)) -
             _QB(ii,jj,kk));
      }
    }
  }

  /// Calculate B-cell face area normals.
  ///
  /// Note: there one less face in two out of three directions. For example,
  /// AB1 has NI+1 faces in "I" direction but only NJ and NK in J and K
  /// directions, respectively.
  TVectorField &AB0 = _AB[0], &AB1 = _AB[1], &AB2 = _AB[2];
  for(size_t ii=0; ii<nei; ++ii)
  {
    ip = ii+1;
    for(size_t jj=0; jj<nej; ++jj)
    {
      jp = jj+1;
      for(size_t kk=0; kk<nek; ++kk)
      {
        kp = kk+1;

        if (jj < nej-1 && kk < nek-1)
          AB0(ii,jj,kk) = (float)0.5 * (
            DB1(ii,jj,kk) % DB2(ii,jj,kk) + DB1(ii,jp,kk) % DB2(ii,jj,kp) );

        if (ii < nei-1 && kk < nek-1)
          AB1(ii,jj,kk) = (float)0.5 * (
            DB2(ii,jj,kk) % DB0(ii,jj,kk) + DB2(ip,jj,kk) % DB0(ii,jj,kp) );

        if (ii < nei-1 && jj < nej-1)
          AB2(ii,jj,kk) = (float)0.5 * (
            DB0(ii,jj,kk) % DB1(ii,jj,kk) + DB0(ii,jp,kk) % DB1(ip,jj,kk) );
      }
    }
  }

}
