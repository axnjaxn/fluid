#include <byteimage/byteimage_sdl2.h>
#include <byteimage/matrix.h>

class FluidSim {
protected:
  Matrix N[9], p, ux, uy;
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

  void setEq(int r, int c) {
    for (int i = 0; i < 9; i++) 
      N[i].at(r, c) = EQ * w[i];
  }

  void emitAt(int r, int c, double power = 16.0) {
    if (r < 0 || r >= rows() || c < 0 || c >= cols()) return;
    for (int i = 0; i < 9; i++) 
      N[i].at(r, c) = (EQ + power) * w[i];
  }

  double pressureAt(int r, int c) const {return p.at(r, c);}
};

#include <byteimage/render.h>
class FluidDisplay : public ByteImageDisplay {
protected:
  ByteImage canvas;
  FluidSim sim;
  int sc;

  int rate;
  bool emitting;
  int mx, my;

  void mapColor(double v, ByteImage::BYTE& r, ByteImage::BYTE& g, ByteImage::BYTE& b) {
    r = g = b = 0;
    v -= sim.EQ;
    if (v > 1.0) r = ByteImage::clip(32.0 * log2(v));
    else if (v < 1.0) g = b = ByteImage::clip(32.0 * log2(-v));
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
      case SDLK_w:
	printf("Set omega (current value: %.2lf)\n", sim.omega);
	scanf("%lf", &sim.omega);
	break;
      case SDLK_SPACE:
	sim.step();
	break;
      case SDLK_BACKSPACE:
	sim = FluidSim(sim.rows(), sim.cols());
	break;
      }
    }
    ByteImageDisplay::handleEvent(event);
  }

  void render() {
    ByteImage::BYTE R, G, B;
    for (int r = 0; r < sim.rows(); r++)
      for (int c = 0; c < sim.cols(); c++) {
	mapColor(sim.pressureAt(r, c), R, G, B);
	DrawRect(canvas, c * sc, r * sc, sc, sc, R, G, B);
      }
    updateImage(canvas);
  }

  void emit() {
    if (emitting)
      for (int i = -2; i <= 2; i++)
	for (int j = -2; j <= 2; j++)
	  sim.emitAt(my + i, mx + j);
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
    
    emitting = 0;
    rate = 1;
  }
};

int main(int argc, char* argv[]) {
  FluidDisplay(200, 100, 4).main();
  return 0;
}
