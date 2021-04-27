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
#include "SDL2Auxiliary.h"

using namespace std;
using glm::vec3;

// --------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
SDL2Aux *sdlAux;

// --------------------------------------------------------
// FUNCTION DECLARATIONS

void Draw(vector<vec3> leftSide, vector<vec3> rightSide);
void Interpolate(float a, float b, vector<float>& result);
void Interpolate(vec3 a, vec3 b, vector<vec3>& result);

// --------------------------------------------------------
// FUNCTION DEFINITIONS

int main( int argc, char* argv[] )
{
	sdlAux = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
	
	vec3 topLeft(1, 0, 0);
	vec3 topRight(0, 0, 1);
	vec3 bottomLeft(1, 1, 0);
	vec3 bottomRight(0, 1, 0);
	vector<vec3> leftSide (SCREEN_HEIGHT);
	vector<vec3> rightSide (SCREEN_HEIGHT);
	Interpolate(topLeft, bottomLeft, leftSide);
	Interpolate(topRight, bottomRight, rightSide);
	
	while (!sdlAux->quitEvent()) {
    	Draw(leftSide, rightSide);
  	}

	sdlAux->saveBMP("colors.bmp");
  	return 0;
}

void Draw(vector<vec3> leftSide, vector<vec3> rightSide)
{
	sdlAux->clearPixels();
	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		vector<vec3> rowColors(SCREEN_WIDTH);
		Interpolate(leftSide[y], rightSide[y], rowColors);
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			sdlAux->putPixel(x, y, rowColors[x]);
		}
	}

	sdlAux->render();
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
