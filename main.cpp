#include "fluidsim.h"
#include <byteimage/byteimage_sdl2.h>
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
    CURL,
    SPEED
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
    v *= 50.0;
    if (v > 0.0) r = ByteImage::clip(255.0 * v);
    else if (v < 0.0) g = b = ByteImage::clip(255.0 * -v);
  }
  void mapSpeedColor(double v, ByteImage::BYTE& r, ByteImage::BYTE& g, ByteImage::BYTE& b) {
    r = g = b = 0;
    v *= 10.0;
    r = ByteImage::clip(255.0 * v);
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
      case SDLK_s:
	rendermode = SPEED;
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
	case SPEED:
	  mapSpeedColor(sim.speedAt(r, c), R, G, B);
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
