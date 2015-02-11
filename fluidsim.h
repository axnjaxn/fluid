#ifndef _BPJ_FLUID_FLUIDSIM_H
#define _BPJ_FLUID_FLUIDSIM_H

#include <byteimage/byteimage.h>
#include <byteimage/matrix.h>

class FluidSim {
protected:
  Matrix N[9], p, ux, uy;
  ByteImage wall;
  double w[9];

  inline static double sq(double d) {return d * d;}
  static double dot(int i, double ux, double uy);

  void setWeights();

public:
  static const double EQ = 100.0;
  double omega;

  FluidSim();
  FluidSim(int nr, int nc);

  inline int rows() const {return N[0].rows();}
  inline int cols() const {return N[0].cols();}
  
  void step();

  void setWall(int r, int c);
  void setEq(int r, int c);
  void emitAt(int r, int c, double power = 24.0);
  void accelAt(int r, int c, double power = 60.0);

  double pressureAt(int r, int c) const;
  double curlAt(int r, int c) const;
  double speedAt(int r, int c) const;
  unsigned char wallAt(int r, int c) const;
};

#endif
