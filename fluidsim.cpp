#include "fluidsim.h"
#include <cstring>
#include <cmath>

/*
 * A note on the axes
 *
 * +-----------> c, x
 * |
 * |     2
 * |   6   5
 * | 3   0   1
 * |   7   8
 * |     4
 * |
 * v
 *
 * r, y
 */

double FluidSim::dot(int i, double ux, double uy) {
  switch (i) {
  default:
  case 0: return 0.0;
  case 1: return ux;
  case 2: return -uy;
  case 3: return -ux;
  case 4: return uy;
  case 5: return ux - uy;
  case 6: return -ux -uy;
  case 7: return -ux + uy;
  case 8: return ux + uy;
  }
}

void FluidSim::setWeights() {
  w[0] = 4.0 / 9;
  w[1] = w[2] = w[3] = w[4] = 1.0 / 9;
  w[5] = w[6] = w[7] = w[8] = 1.0 / 36;
}

FluidSim::FluidSim() {
  omega = 1.0;
  setWeights();
}

FluidSim::FluidSim(int nr, int nc) {
  omega = 1.0;
  setWeights();

  wall = ByteImage(nr, nc);
  p = ux = uy = Matrix(nr, nc);

  for (int i = 0; i < 9; i++) N[i] = p;
  
  for (int r = 0; r < nr; r++)
    for (int c = 0; c < nc; c++) {
      setEq(r, c);
      p.at(r, c) = EQ;
    }
}

void FluidSim::step() {
  const int nr = rows(), nc = cols();

  /* Streaming step */

  //Buffer for horizontal movement
  double* d = new double [nc];

  //Right motion
  for (int r = 0; r < nr; r++) {
    memcpy(d, N[1].getArray() + r * nc, (nc - 1) * sizeof(double));
    memcpy(N[1].getArray() + r * nc + 1, d, (nc - 1) * sizeof(double));
  }

  //Upward motion
  for (int r = 0; r < nr - 1; r++)
    memcpy(N[2].getArray() + r * nc, N[2].getArray() + (r + 1) * nc, nc * sizeof(double));
    
  //Left motion
  for (int r = 0; r < nr; r++) {
    memcpy(d, N[3].getArray() + r * nc + 1, (nc - 1) * sizeof(double));
    memcpy(N[3].getArray() + r * nc, d, (nc - 1) * sizeof(double));
  }

  //Downward motion
  for (int r = nr - 1; r > 0; r--)
    memcpy(N[4].getArray() + r * nc, N[4].getArray() + (r - 1) * nc, nc * sizeof(double));

  //Up-right motion
  for (int r = 0; r < nr - 1; r++)
    memcpy(N[5].getArray() + r * nc + 1, N[5].getArray() + (r + 1) * nc, (nc - 1) * sizeof(double));

  //Up-left motion
  for (int r = 0; r < nr - 1; r++)
    memcpy(N[6].getArray() + r * nc, N[6].getArray() + (r + 1) * nc + 1, (nc - 1) * sizeof(double));

  //Down-left motion
  for (int r = nr - 1; r > 0; r--)
    memcpy(N[7].getArray() + r * nc, N[7].getArray() + (r - 1) * nc + 1, (nc - 1) * sizeof(double));

  //Down-right motion
  for (int r = nr - 1; r > 0; r--)
    memcpy(N[8].getArray() + r * nc + 1, N[8].getArray() + (r - 1) * nc, (nc - 1) * sizeof(double));

  delete [] d;

  //Check walls
  for (int r = 1; r < nr - 1; r++)
    for (int c = 1; c < nc - 1; c++)
      if (isWall(r, c)) {
	N[3].at(r, c - 1) += N[1].at(r, c);
	N[4].at(r + 1, c) += N[2].at(r, c);
	N[1].at(r, c + 1) += N[3].at(r, c);
	N[2].at(r - 1, c) += N[4].at(r, c);
	N[7].at(r + 1, c - 1) += N[5].at(r, c);
	N[8].at(r + 1, c + 1) += N[6].at(r, c);
	N[5].at(r - 1, c + 1) += N[7].at(r, c);
	N[6].at(r - 1, c - 1) += N[8].at(r, c);
	for (int i = 0; i < 9; i++)
	  N[i].at(r, c) = 0.0;
      }

  /* Reset edges */
  for (int r = 0; r < nr; r++) {
    setEq(r, 0);
    setEq(r, nc - 1);
  }
  for (int c = 0; c < nc; c++) {
    setEq(0, c);
    setEq(nr - 1, c);
  }

  /* Collision step */
    
  double Neq;
  for (int r = 0; r < nr; r++)
    for (int c = 0; c < nc; c++) {
      if (isWall(r, c)) {
	p.at(r, c) = ux.at(r, c) = uy.at(r, c) = 0.0;
	continue;
      }

      //Compute density as a sum of all molecules in a cell
      p.at(r, c) = 0.0;
      for (int i = 0; i < 9; i++) {
	if (N[i].at(r, c) < 0.0) N[i].at(r, c) = 0.0;
	p.at(r, c) += N[i].at(r, c);
      }
      if (p.at(r, c) <= 1.0)
	for (int i = 0; i < 9; i++) {
	  N[i].at(r, c) = w[i];
	  p.at(r, c) = 1.0;
	}
      
      //Compute flow velocity
      if (!isFixedVel(r, c)) {
	ux.at(r, c) = (N[1].at(r, c) + N[5].at(r, c) + N[8].at(r, c) - N[3].at(r, c) - N[6].at(r, c) - N[7].at(r, c)) / p.at(r, c);
	uy.at(r, c) = -(N[2].at(r, c) + N[5].at(r, c) + N[6].at(r, c) - N[4].at(r, c) - N[7].at(r, c) - N[8].at(r, c)) / p.at(r, c);
      }

      //Collision and relaxation
      for (int i = 0; i < 9; i++) {
	Neq = p.at(r, c) * w[i] * (1 + 3 * dot(i, ux.at(r, c), uy.at(r, c)) + 4.5 * sq(dot(i, ux.at(r, c), uy.at(r, c))) - 1.5 * (sq(ux.at(r, c)) + sq(uy.at(r, c))));
	N[i].at(r, c) += omega * (Neq - N[i].at(r, c));
      }
    }
}

void FluidSim::setWall(int r, int c) {
  if (r < 1 || r >= rows() - 1 || c < 1 || c >= cols() - 1) return;
  for (int i = 0; i < 9; i++) N[i].at(r, c) = 0;
  wall.at(r, c) = 0xFF;
}

void FluidSim::setEq(int r, int c) {
  for (int i = 0; i < 9; i++) 
    N[i].at(r, c) = EQ * w[i];
}

void FluidSim::emitAt(int r, int c, double power) {
  if (r < 0 || r >= rows() || c < 0 || c >= cols() || wall.at(r, c)) return;
  for (int i = 0; i < 9; i++) 
    N[i].at(r, c) = (EQ + power) * w[i];
}

void FluidSim::accelAt(int r, int c, double power) {
  if (r < 0 || r >= rows() || c < 0 || c >= cols() || wall.at(r, c)) return;
  setEq(r, c);
  ux.at(r, c) = power;
  uy.at(r, c) = 0.0;
  wall.at(r, c) = FIXED_VEL;
}

void FluidSim::setWindTunnel(double power) {
  for (int r = 0; r < rows(); r++) {
    accelAt(r, 0, power);
    accelAt(r, cols() - 1, power);
  }    
}

double FluidSim::pressureAt(int r, int c) const {return p.at(r, c);}

double FluidSim::curlAt(int r, int c) const {
  if (r <= 0 || r >= rows() - 1 || c <= 0 || c >= cols() - 1) return 0.0;
  return uy.at(r, c + 1) - uy.at(r, c - 1) - ux.at(r + 1, c) + ux.at(r - 1, c);
}
double FluidSim::speedAt(int r, int c) const {
  return sqrt(sq(ux.at(r, c)) + sq(uy.at(r, c)));
}

bool FluidSim::isWall(int r, int c) const {
  return wall.at(r, c) & WALL;
}

bool FluidSim::isFixedVel(int r, int c) const {
  return wall.at(r, c) & FIXED_VEL;
}
