#include "fluidsim.h"
#include <byteimage/byteimage_sdl2.h>
#include <byteimage/render.h>
#include <vector>

class FluidDisplay : public ByteImageDisplay {
protected:
  ByteImage canvas;
  FluidSim sim;
  int sc, radius;
  std::vector<Matrix> trackers;

  enum {
    EMIT,
    ACCEL,
    WALL,
    TRACKER
  } emitmode;

  enum {
    PRESSURE,
    CURL,
    SPEED
  } rendermode;

  int rate;
  bool emitting, drawing, showgrid;
  int mx, my, nx, ny;

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
      if (event.button.button == SDL_BUTTON_LEFT)
	emitting = 1;
      else if (event.button.button == SDL_BUTTON_RIGHT)
	drawing = 1;
      mx = nx = event.button.x / sc;
      my = ny = event.button.y / sc;
    }
    else if (event.type == SDL_MOUSEBUTTONUP) {
      if (event.button.button == SDL_BUTTON_LEFT)
	emitting = 0;
      else if (event.button.button == SDL_BUTTON_RIGHT) {
	drawing = 0;
	Matrix d = makePoint(mx - nx, my - ny, 0);
	int n = length(d);
	for (int i = 0; i <= n; i++)
	  sim.setWall((int)(ny + d.at(1) * i / n), (int)(nx + d.at(0) * i / n));
      }
    }
    else if (event.type == SDL_MOUSEMOTION) {
      mx = event.motion.x / sc;
      my = event.motion.y / sc;
    }
    else if (event.type == SDL_KEYDOWN) {
      switch (event.key.keysym.sym) {
      case SDLK_0: rate = 0; break;
      case SDLK_1: rate = 1; break;
      case SDLK_2: rate = 2; break;
      case SDLK_3: rate = 4; break;
      case SDLK_4: rate = 6; break;
      case SDLK_5: rate = 8; break;
      case SDLK_p: rendermode = PRESSURE; break;
      case SDLK_e: emitmode = EMIT; break;
      case SDLK_a: emitmode = ACCEL; break;
      case SDLK_w: emitmode = WALL; break;
      case SDLK_r: emitmode = TRACKER; break;
      case SDLK_c: rendermode = CURL; break;
      case SDLK_s: rendermode = SPEED; break;
      case SDLK_x: showgrid = !showgrid; break;
      case SDLK_o:
	printf("Set omega (current value: %.2lf)\n", sim.omega);
	scanf("%lf", &sim.omega);
	break;
      case SDLK_SPACE:
	rate = 0;
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
      case SDLK_t:
	sim.setWindTunnel();
	break;
      }
    }
    ByteImageDisplay::handleEvent(event);
  }

  void render() {
    ByteImage::BYTE R, G, B;
    Matrix v, v1, color = makeColor(255, 255, 255);
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
	if (sim.isWall(r, c)) R = G = B = 255;
	DrawRect(canvas, c * sc, r * sc, sc, sc, R, G, B);

	if (showgrid) {
	  v = makePoint(c + 0.5, r + 0.5, 1.0 / sc);
	  v1 = v + 25.0 * makePoint(sim.xVel(r, c), -sim.yVel(r, c), 0.0);
	  DrawLine(canvas, v, v1, color);
	}
      }

    for (int i = 0; i < trackers.size(); i++) {
      DrawPoint(canvas, trackers[i], color, 3);
      DrawPoint(canvas, trackers[i], 0.5 * color, 3);
    }
    
    if (drawing)
      DrawLine(canvas, makePoint(nx + 0.5, ny + 0.5, 1.0 / sc), makePoint(mx + 0.5, my + 0.5, 1.0 / sc), makeColor(255, 255, 255));
    updateImage(canvas);
  }

  void emit() {
    if (!emitting) return;

    if (emitmode == TRACKER) {
      trackers.push_back(makePoint(mx, my, 1.0 / sc));
      return;
    }

    for (int i = -radius; i <= radius; i++)
      for (int j = -radius; j <= radius; j++)
	if (i * i + j * j <= radius * radius)
	  switch (emitmode) {
	  default:
	  case EMIT:
	    sim.emitAt(my + i, mx + j);
	    break;
	  case ACCEL:
	    sim.accelAt(my + i, mx + j);
	    break;
	  case WALL:
	    sim.setWall(my + i, mx + j);
	    break;
	  }
  }

  void moveTrackers() {
    double x, y;
    int r, c;
    
    for (int i = 0; i < trackers.size(); i++) {
      x = trackers[i].at(0);
      y = trackers[i].at(1);
      r = (int)(y + 0.5);
      c = (int)(x + 0.5);
      
      x += 10 * sim.xVel(r, c);
      y += 10 * sim.yVel(r, c);
      r = (int)(y + 0.5);
      c = (int)(x + 0.5);

      if (c < 1 || c >= sim.cols() - 1 || r < 1 || r >= sim.rows() - 1) {
	trackers.erase(trackers.begin() + i--);
	continue;
      }

      trackers[i] = makePoint(x, y, 1.0 / sc);
    }
  }

  void update() {
    for (int i = 0; i < rate; i++) {
      emit();
      sim.step();
      moveTrackers();
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
    emitting = drawing = showgrid = 0;
    radius = 3;
    rate = 1;
  }
};

int main(int argc, char* argv[]) {
  FluidDisplay(200, 100, 4).main();
  return 0;
}
