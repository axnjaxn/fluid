#include "fluidsim.h"
#include <byteimage/byteimage_sdl2.h>
#include <byteimage/video.h>
#include <byteimage/palette.h>
#include <byteimage/render.h>
#include <vector>

using namespace byteimage;

class FluidDisplay : public Display {
protected:
  ByteImage canvas;
  FluidSim sim;
  int sc, radius;
  std::vector<Pt2f> trackers;

  CachedPalette pressure_palette, curl_palette, speed_palette;

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
  
  void mapPressureColor(double v, Color& color) {
    color = pressure_palette.inRange(32.0 * (v - sim.EQ) / 255.0);
  }
  void mapCurlColor(double v, Color& color) {
    color = curl_palette.inRange(16.0 * v);
  }
  void mapSpeedColor(double v, Color& color) {
    color = speed_palette.inUnit(4.0 * v);
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
	Pt2f d(mx - nx, my - ny);
	int n = length(d);
	for (int i = 0; i <= n; i++)
	  sim.setWall((int)(ny + d.y * i / n), (int)(nx + d.x * i / n));
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
	trackers.clear();
	break;
      case SDLK_UP:
	printf("Radius: %d\n", ++radius);
	break;
      case SDLK_DOWN:
	if (radius > 1)
	  printf("Radius: %d\n", --radius);
	break;
      case SDLK_LEFT:
	if (sim.omega > 0.05) sim.omega -= 0.05;
	printf("Omega: %.2lf\n", sim.omega);
	break;
      case SDLK_RIGHT:
	if (sim.omega < 1.95) sim.omega += 0.05;
	printf("Omega: %.2lf\n", sim.omega);
	break;
      case SDLK_t:
	sim.setWindTunnel();
	break;
      case SDLK_RETURN:
	if (recording) stopRecording();
	else startRecording();
	break;
      }
    }
    Display::handleEvent(event);
  }

  void render() {
    Color cell_color;
    Pt2f v, v1;
    Color white(255, 255, 255), gray(128, 128, 128);
    for (int r = 0; r < sim.rows(); r++)
      for (int c = 0; c < sim.cols(); c++) {
	switch (rendermode) {
	case PRESSURE:
	  mapPressureColor(sim.pressureAt(r, c), cell_color);
	  break;
	case CURL:
	  mapCurlColor(sim.curlAt(r, c), cell_color);
	  break;
	case SPEED:
	  mapSpeedColor(sim.speedAt(r, c), cell_color);
	  break;
	}	
	if (sim.isWall(r, c)) cell_color = Color(255, 255, 255);
	DrawRect(canvas, c * sc, r * sc, sc, sc, cell_color.r, cell_color.g, cell_color.b);

	if (showgrid) {
	  v = Pt2f(c + 0.5, r + 0.5) * sc;
	  v1 = v + 25.0 * Pt2f(sim.xVel(r, c), -sim.yVel(r, c));
	  DrawLine(canvas, v, v1, white);
	}
      }

    for (int i = 0; i < trackers.size(); i++) {
      DrawPoint(canvas, trackers[i], white, 3);
      DrawPoint(canvas, trackers[i], gray, 1);
    }
    
    if (drawing)
      DrawLine(canvas, Pt2f(nx + 0.5, ny + 0.5) * sc, Pt2f(mx + 0.5, my + 0.5) * sc, white);
    updateImage(canvas);
  }

  void emit() {
    if (!emitting) return;

    if (emitmode == TRACKER) {
      trackers.push_back(Pt2f(mx, my) * sc);
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
      x = trackers[i].x;
      y = trackers[i].y;
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

      trackers[i] = Pt2f(x, y) * sc;
    }
  }

  void update() {
    for (int i = 0; i < rate; i++) {
      emit();
      sim.step();
      moveTrackers();
    }
    render();
    if (recording) recordFrame();
    Display::update();
    if (exitflag && recording) stopRecording(); 
  }

  VideoWriter* writer = nullptr;
  bool recording;

  void startRecording() {
    if (recording) return;
    recording = 1;

    char fn[256];
    sprintf(fn, "%d.avi", (int)time(NULL));

    printf("Writing to %s...\n", fn);
    writer = new VideoWriter(fn, canvas.nr, canvas.nc, 30);
  }

  void stopRecording() {
    if (!recording) return;
    recording = 0;

    delete writer;
    printf("Finished writing video.\n");
  }
  void recordFrame() {
    writer->write(canvas);
  }

public:
  FluidDisplay(int w, int h, int sc) : Display(h * sc, w * sc, "Fluid Simulation by Brian Jackson") {
    this->sc = sc;
    canvas = ByteImage(h * sc, w * sc, 3);
    updateImage(canvas);

    LinearPalette pal(3);
    pal[0] = Color(0, 255, 255);
    pal[1] = Color(0, 0, 0);
    pal[2] = Color(255, 0, 0);
    pressure_palette = pal.cache(256);
    pal[0] = Color(0, 0, 0);
    pal[1] = Color(255, 128, 64);
    pal[2] = Color(255, 255, 255);
    speed_palette = pal.cache(256);
    curl_palette = LinearPalette::jet().cache(256);
    
    sim = FluidSim(h, w);

    emitmode = EMIT;
    rendermode = PRESSURE;
    emitting = drawing = showgrid = 0;
    radius = 3;
    rate = 1;

    recording = 0;
  }
};

int main(int argc, char* argv[]) {
  FluidDisplay(200, 100, 4).main();
  return 0;
}
