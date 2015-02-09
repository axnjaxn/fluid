#include <byteimage/byteimage_sdl2.h>
#include <byteimage/matrix.h>

class FluidSim {
protected:
  double T;
  Matrix N[9], p, ux, uy;

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

public:
  FluidSim() {
    T = 1.0;
  }
  FluidSim(int nr, int nc) {
    T = 1.0;

    N[0] = Matrix(nr, nc);
    for (int r = 0; r < nr; r++)
      for (int c = 0; c < nc; c++)
	N[0].at(r, c) = 1.0;

    for (int i = 1; i < 9; i++) N[i] = Matrix(nr, nc);

    p = ux = uy = Matrix(nr, nc);
  }

  inline int rows() const {return N[0].rows();}
  inline int cols() const {return N[0].cols();}
  
  void step() {
    const int nr = rows(), nc = cols();

    /* Collision step */
    
    const double w[9] = {4.0 / 9, 
			 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
			 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};
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
	  N[i].at(r, c) += (Neq - N[i].at(r, c)) / T;
	}
      }    

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
  }

  void emitAt(int r, int c, double level = 10.0) {
    const double w[9] = {4.0 / 9, 
			 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
			 1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

    if (p.at(r, c) > level) return;

    for (int i = 0; i < 9; i++) 
      N[i].at(r, c) = w[i] * level;
  }

  double pressureAt(int r, int c) const {return p.at(r, c);}
};

#include <byteimage/render.h>
class FluidDisplay : public ByteImageDisplay {
protected:
  ByteImage canvas;
  FluidSim sim;
  int rate;
  bool emitting;
  int mx, my;

  void mapColor(double v, ByteImage::BYTE& r, ByteImage::BYTE& g, ByteImage::BYTE& b) {
    if (v > 1.0) r = ByteImage::clip(64.0 * log2(v));
    else g = b = ByteImage::clip(64.0 * log2(1.0 / v));
  }

  void handleEvent(SDL_Event event) {
    if (event.type == SDL_MOUSEBUTTONDOWN) {
      emitting = 1;
      mx = event.button.x;
      my = event.button.y;
    }
    else if (event.type == SDL_MOUSEBUTTONUP) {
      emitting = 0;
    }
    else if (event.type == SDL_MOUSEMOTION) {
      mx = event.motion.x;
      my = event.motion.y;
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
      case SDLK_SPACE:
	sim.step();
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
	DrawRect(canvas, c, r, 1, 1, R, G, B);
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
  FluidDisplay(int w, int h) : ByteImageDisplay(h, w, "Fluid Simulation by Brian Jackson") {
    canvas = ByteImage(h, w, 3);
    updateImage(canvas);

    sim = FluidSim(h, w);
    
    emitting = 0;
    rate = 1;
  }
};

int main(int argc, char* argv[]) {
  FluidDisplay(200, 100).main();
  return 0;
}
