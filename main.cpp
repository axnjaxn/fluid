#include <byteimage/byteimage_sdl2.h>
#include <byteimage/matrix.h>

class FluidSim {
protected:
  Matrix N[9], p, ux, uy;
  ByteImage wall;
  double w[9];

  inline static double sq(double d) {return d * d;}
  static double dot(int i, double ux, double uy) {
    switch (i) {
    default:
    case 0: return 0.0;
    case 1: return ux;
    case 2: return uy;
    case 3: return -ux;
    case 4: return -uy;
    case 5: return ux + uy;
    case 6: return uy - ux;
    case 7: return -ux - uy;
    case 8: return ux - uy;
    }
  }

  void setWeights() {
    w[0] = 4.0 / 9;
    w[1] = w[2] = w[3] = w[4] = 1.0 / 9;
    w[5] = w[6] = w[7] = w[8] = 1.0 / 36;
  }

public:
  static const double EQ = 100.0;
  double omega;

  FluidSim() {
    omega = 1.0;
    setWeights();
  }

  FluidSim(int nr, int nc) {
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

  inline int rows() const {return N[0].rows();}
  inline int cols() const {return N[0].cols();}
  
  void step() {
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
	if (wall.at(r, c)) {
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
	if (wall.at(r, c)) {
	  p.at(r, c) = ux.at(r, c) = uy.at(r, c) = 0.0;
	  continue;
	}

	//Compute density as a sum of all molecules in a cell
	p.at(r, c) = 0.0;
	for (int i = 0; i < 9; i++)
	  p.at(r, c) += N[i].at(r, c);

	//Compute flow velocity
	ux.at(r, c) = (N[1].at(r, c) + N[5].at(r, c) + N[8].at(r, c) - N[3].at(r, c) - N[6].at(r, c) - N[7].at(r, c)) / p.at(r, c);
	uy.at(r, c) = (N[2].at(r, c) + N[5].at(r, c) + N[6].at(r, c) - N[4].at(r, c) - N[7].at(r, c) - N[8].at(r, c)) / p.at(r, c);

	//Collision and relaxation
	for (int i = 0; i < 9; i++) {
	  Neq = p.at(r, c) * w[i] * (1 + 3 * dot(i, ux.at(r, c), uy.at(r, c)) + 4.5 * sq(dot(i, ux.at(r, c), uy.at(r, c))) - 1.5 * (sq(ux.at(r, c)) + sq(uy.at(r, c))));
	  N[i].at(r, c) += omega * (Neq - N[i].at(r, c));
	}
      }    

  }

  void setWall(int r, int c) {
    if (r < 1 || r >= rows() - 1 || c < 1 || c >= cols() - 1) return;
    for (int i = 0; i < 9; i++) N[i].at(r, c) = 0;
    wall.at(r, c) = 0xFF;
  }
  void setEq(int r, int c) {
    for (int i = 0; i < 9; i++) 
      N[i].at(r, c) = EQ * w[i];
  }

  void emitAt(int r, int c, double power = 24.0) {
    if (r < 0 || r >= rows() || c < 0 || c >= cols() || wall.at(r, c)) return;
    for (int i = 0; i < 9; i++) 
      N[i].at(r, c) = (EQ + power) * w[i];
  }
  void accelAt(int r, int c, double power = 50.0) {
    if (r < 0 || r >= rows() || c < 0 || c >= cols() || wall.at(r, c)) return;
    setEq(r, c);
    N[2].at(r, c) = (EQ + power) * w[2];
  }

  double pressureAt(int r, int c) const {return p.at(r, c);}
  double curlAt(int r, int c) const {
    if (r <= 0 || r >= rows() - 1 || c <= 0 || c >= cols() - 1) return 0.0;
    return uy.at(r, c + 1) - uy.at(r, c - 1) - ux.at(r + 1, c) + ux.at(r - 1, c);
  }
  unsigned char wallAt(int r, int c) const {
    return wall.at(r, c);
  }
};

#include <byteimage/render.h>
class FluidDisplay : public ByteImageDisplay {
protected:
  ByteImage canvas;
  FluidSim sim;
  int sc, radius;

  enum {
    EMIT,
    ACCEL,
    WALL
  } emitmode;

  enum {
    PRESSURE,
    CURL    
  } rendermode;

  int rate;
  bool emitting;
  int mx, my;

  void mapPressureColor(double v, ByteImage::BYTE& r, ByteImage::BYTE& g, ByteImage::BYTE& b) {
    r = g = b = 0;
    v -= sim.EQ;
    if (v > 0.0) r = ByteImage::clip(32.0 * v);
    else if (v < 0.0) g = b = ByteImage::clip(32.0 * -v);
  }
  void mapCurlColor(double v, ByteImage::BYTE& r, ByteImage::BYTE& g, ByteImage::BYTE& b) {
    r = g = b = 0;
    v *= 25.0;
    if (v > 0.0) r = ByteImage::clip(255.0 * v);
    else if (v < 0.0) g = b = ByteImage::clip(255.0 * -v);
  }

  void handleEvent(SDL_Event event) {
    if (event.type == SDL_MOUSEBUTTONDOWN) {
      emitting = 1;
      mx = event.button.x / sc;
      my = event.button.y / sc;
    }
    else if (event.type == SDL_MOUSEBUTTONUP) {
      emitting = 0;
    }
    else if (event.type == SDL_MOUSEMOTION) {
      mx = event.motion.x / sc;
      my = event.motion.y / sc;
    }
    else if (event.type == SDL_KEYDOWN) {
      switch (event.key.keysym.sym) {
      case SDLK_0:
	rate = 0;
	break;
      case SDLK_1:
	rate = 1;
	break;
      case SDLK_2:
	rate = 2;
	break;
      case SDLK_3:
	rate = 4;
	break;
      case SDLK_4:
	rate = 8;
	break;
      case SDLK_5:
	rate = 16;
	break;
      case SDLK_p:
	rendermode = PRESSURE;
	break;
      case SDLK_e:
	emitmode = EMIT;
	break;
      case SDLK_a:
	emitmode = ACCEL;
	break;
      case SDLK_w:
	emitmode = WALL;
	break;
      case SDLK_c:
	rendermode = CURL;
	break;
      case SDLK_o:
	printf("Set omega (current value: %.2lf)\n", sim.omega);
	scanf("%lf", &sim.omega);
	break;
      case SDLK_SPACE:
	sim.step();
	break;
      case SDLK_BACKSPACE:
	sim = FluidSim(sim.rows(), sim.cols());
	break;
      case SDLK_UP:
	printf("Radius: %d\n", ++radius);
	break;
      case SDLK_DOWN:
	if (radius > 1)
	  printf("Radius: %d\n", --radius);
	break;
      }
    }
    ByteImageDisplay::handleEvent(event);
  }

  void render() {
    ByteImage::BYTE R, G, B;
    for (int r = 0; r < sim.rows(); r++)
      for (int c = 0; c < sim.cols(); c++) {
	switch (rendermode) {
	case PRESSURE:
	  mapPressureColor(sim.pressureAt(r, c), R, G, B);
	  break;
	case CURL:
	  mapCurlColor(sim.curlAt(r, c), R, G, B);
	  break;
	}	
	if (sim.wallAt(r, c)) R = G = B = 255;
	DrawRect(canvas, c * sc, r * sc, sc, sc, R, G, B);
      }
    updateImage(canvas);
  }

  void emit() {
    if (!emitting) return;

    if (emitmode == EMIT) {
      for (int i = -radius; i <= radius; i++)
	for (int j = -radius; j <= radius; j++)
	  if (i * i + j * j <= radius * radius)
	    sim.emitAt(my + i, mx + j);
    }
    else if (emitmode == ACCEL) {
      for (int i = -radius; i <= radius; i++)
	for (int j = -radius; j <= radius; j++)
	  if (i * i + j * j <= radius * radius)
	    sim.accelAt(my + i, mx + j);
    }
    else if (emitmode == WALL) {
      for (int i = -radius; i <= radius; i++)
	for (int j = -radius; j <= radius; j++)
	  sim.setWall(my + i, mx + j);
    }
  }

  void update() {
    for (int i = 0; i < rate; i++) {
      emit();
      sim.step();
    }
    render();
    ByteImageDisplay::update();
  }

public:
  FluidDisplay(int w, int h, int sc) : ByteImageDisplay(h * sc, w * sc, "Fluid Simulation by Brian Jackson") {
    this->sc = sc;
    canvas = ByteImage(h * sc, w * sc, 3);
    updateImage(canvas);

    sim = FluidSim(h, w);

    emitmode = EMIT;
    rendermode = PRESSURE;
    emitting = 0;
    radius = 3;
    rate = 1;
  }
};

int main(int argc, char* argv[]) {
  FluidDisplay(200, 100, 4).main();
  return 0;
}
