#ifndef _BPJ_FLUID_FLUIDSIM_H
#define _BPJ_FLUID_FLUIDSIM_H

#include <byteimage/byteimage.h>
#include <byteimage/matrix.h>

using byteimage::ByteImage;
using byteimage::Matrix;

class FluidSim {
protected:
  Matrix N[9], p, ux, uy;

  static constexpr unsigned char WALL      = 0x01;
  static constexpr unsigned char FIXED_VEL = 0x02;

  ByteImage wall;
  double w[9];

  inline static double sq(double d) {return d * d;}
  static double dot(int i, double ux, double uy);

  void setWeights();

public:
  static constexpr double EQ = 100.0;
  double omega;

  FluidSim();
  FluidSim(int nr, int nc);

  inline int rows() const {return N[0].rows();}
  inline int cols() const {return N[0].cols();}

  void step();

  void setWall(int r, int c);
  void setEq(int r, int c);
  void emitAt(int r, int c, double power = 24.0);
  void accelAt(int r, int c, double power = 0.2);
  void setWindTunnel(double power = 0.05);

  double pressureAt(int r, int c) const;
  double curlAt(int r, int c) const;
  double speedAt(int r, int c) const;
  inline double xVel(int r, int c) const {return ux.at(r, c);}
  inline double yVel(int r, int c) const {return uy.at(r, c);}

  bool isWall(int r, int c) const;
  bool isFixedVel(int r, int c) const;
};

#endif
