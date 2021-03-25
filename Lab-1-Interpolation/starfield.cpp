// Introduction lab that covers:
// * C++
// * SDL
// * 2D graphics
// * Plotting pixels
// * Video memory
// * Color representation
// * Linear interpolation
// * glm::vec3 and std::vector

#include "SDL.h"
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include "SDL2auxiliary.h"

using namespace std;
using glm::vec3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
SDL2Aux *sdlAux;

// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw(vector<vec3>& stars);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

void Update(vector<vec3>& stars, float dt);

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main( int argc, char* argv[] )
{
	sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);

	vector<vec3> stars(1000);
	auto randNum = [](){ return float(rand()) / float(RAND_MAX); };

	for (int i = 0; i < stars.size(); i++){
		stars[i].x = -1 + 2 * randNum();
		stars[i].y = -1 + 2 * randNum();
		stars[i].z = randNum();
	}
	int t = SDL_GetTicks();
	while(!sdlAux->quitEvent())
	{
		// Calculate delta time
		int t2 = SDL_GetTicks();
		float dt = float(t2 - t);
		t = t2;

		Update(stars, dt);
		Draw(stars);
	}
	sdlAux->saveBMP("starfield.bmp");
	return 0;
}

void Draw(vector<vec3>& stars)
{
	sdlAux->clearPixels();
	// Focal length
	float f = SCREEN_HEIGHT / 2.0f;
	// Center points of camera projection
	float c_x = SCREEN_WIDTH / 2.0f;	
	float c_y = SCREEN_HEIGHT / 2.0f;	

	for (size_t i = 0; i < stars.size(); i++) {
		vec3 color = vec3(1,1,1) / (stars[i].z * stars[i].z);  

		float u = f * stars[i].x / stars[i].z + c_x;
		float v = f * stars[i].y / stars[i].z + c_y;
		sdlAux->putPixel(u, v, color);
	}
	sdlAux->render();
}

void Update(vector<vec3>& stars, float dt) {
	vec3 v(0, 0, 1);
	const float max_time = 1000;
	for (vec3& star : stars) {
		if (dt > max_time) dt = max_time;
		else if (dt <= 0) dt = 0;
		star.z -= v.z * dt / max_time;
		if (star.z <= 0)
			star.z += 1;
		if (star.z > 1)
			star.z = 1; 
	}	
}

void Interpolate (float a, float b, vector<float>& result) {
	// Inclusive both a and b
	const int size = result.size() - 1;
	for ( int i = 0; i <= size; i++ ) {
		float p = i / (float) size;
		result[i] = (1 - p) * a + p * b;
	}	
}

void Interpolate (vec3 a, vec3 b, vector<vec3>& result) {
	// Inclusive elements in a and b
	const int size = result.size() - 1;
	for ( int i = 0; i <= size; i++) {
		float p = i / (float) size;
		result[i].x = (1 - p) * a.x + p * b.x;
		result[i].y = (1 - p) * a.y + p * b.y;
		result[i].z = (1 - p) * a.z + p * b.z;
	}
}
