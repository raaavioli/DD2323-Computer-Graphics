#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDL2Auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;

// ----------------------------------------------------------------------------
// STRUCTURES

struct Camera {
	vec3 position;
	float yaw, pitch;
	

	mat3 getRotation() const {
		return glm::mat3(
			glm::cos(this->yaw), 0, glm::sin(this->yaw),
			0.0, 1.0, 0.0,
			-glm::sin(this->yaw), 0.0, glm::cos(this->yaw) 
		) * glm::mat3(
			1.0, 0.0, 0.0,
			0.0, glm::cos(this->pitch), glm::sin(this->pitch),
			0.0, -glm::sin(this->pitch), glm::cos(this->pitch) 
		);
	}
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const int FOCAL_LENGTH = SCREEN_WIDTH / 2;
SDL2Aux* screen;
int t;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(Camera& camera);
void Draw(vector<Triangle> triangles, const Camera& camera);
void VertexShader(const vec3& v, ivec2& p, const Camera& camera);

int main (int argc, char* argv[]) {
	screen = new SDL2Aux (SCREEN_WIDTH, SCREEN_HEIGHT);
	Camera camera = {.position = vec3( 0, 0, -2.05 ), .yaw = 0.0};
	vector<Triangle> triangles;
	LoadTestModel(triangles);
	t = SDL_GetTicks();	// Set start value for timer.

	while(!screen->quitEvent())
	{
		Update(camera);
		Draw(triangles, camera);
	}

	screen->saveBMP("rasterization.bmp");
	return 0;
}

void Update (Camera& camera) {
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	mat3 R = camera.getRotation();
	vec3 right (R[0][0], R[0][1], R[0][2]);
	vec3 forward (R[2][0], R[2][1], R[2][2]);

	const Uint8* keystate = SDL_GetKeyboardState(NULL);

	if (dt > 1000) dt = 1000;
	if(keystate[SDL_SCANCODE_UP])
		camera.pitch += dt / 1000.f;

	if(keystate[SDL_SCANCODE_DOWN])
		camera.pitch -= dt / 1000.f;

	if(keystate[SDL_SCANCODE_RIGHT])
		camera.yaw -= dt / 1000.f;

	if(keystate[SDL_SCANCODE_LEFT])
		camera.yaw += dt / 1000.f;

	if(keystate[SDL_SCANCODE_RSHIFT])
		;

	if(keystate[SDL_SCANCODE_RCTRL])
		;

	if(keystate[SDL_SCANCODE_W])
		camera.position += forward * dt / 1000.f;

	if(keystate[SDL_SCANCODE_S])
		camera.position -= forward * dt / 1000.f;

	if(keystate[SDL_SCANCODE_D])
		camera.position += right * dt / 1000.f;

	if(keystate[SDL_SCANCODE_A])
		camera.position -= right * dt / 1000.f;

	if(keystate[SDL_SCANCODE_E])
		;

	if(keystate[SDL_SCANCODE_Q])
		;
}

void Draw (vector<Triangle> triangles, const Camera& camera) {
	screen->clearPixels();
	
	for(int i = 0; i<triangles.size(); i++ ) {
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		for(int v=0; v<3; ++v) {
			ivec2 projPos;
			VertexShader(vertices[v], projPos, camera);
			vec3 color(1,1,1);
			screen->putPixel(projPos.x, projPos.y, color);
		}
	}
	screen->render();
}


void VertexShader (const vec3& v, ivec2& p, const Camera& camera) {
	vec3 vertex(v);
	vertex -= camera.position;
	vertex = vertex * camera.getRotation();
	p.x = FOCAL_LENGTH * vertex.x / vertex.z + SCREEN_WIDTH / 2;
	p.y = FOCAL_LENGTH * vertex.y / vertex.z + SCREEN_HEIGHT / 2;
}