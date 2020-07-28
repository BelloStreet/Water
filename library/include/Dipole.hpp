// cSpell:disable
#ifndef _DIPOLE_HPP_
#define _DIPOLE_HPP_

class Dipole {
public:
  /// A class that constructs various dipole operators in the length gauge
  Dipole();

  ~Dipole();

  /// Linear polarization $z=r\cos$
  void linear()

      /// Perpendicular polarization $x=r\cos(\phi)\sin(theta),
      /// y=r\sin(\phi)\sin(theta)
      void perpendicular()

      /// Circular polarization $x+iy$
      void circular()

      /// Elliptic polarization
      void elliptic()

          private:
}
#endif
