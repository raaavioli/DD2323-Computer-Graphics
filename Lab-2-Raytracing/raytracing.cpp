#include <iostream>
#include <algorithm>
#include <glm/glm.hpp>
#include "SDL.h"
#include "SDL2Auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const float PI = 3.141592654;
SDL2Aux* screen;
int current_time;

// ----------------------------------------------------------------------------
// STRUCTURES
struct Ray {
	vec3 s;
	vec3 d;
};

struct Light {
	vec3 position;
	vec3 color;
};

struct Intersection
{
      vec3 position;
      float distance;
      int triangleIndex;
};

struct Camera {
	vec3 position;
	float yaw;
	

	mat3 getRotation() {
		return glm::mat3(
			glm::cos(this->yaw), 0, glm::sin(this->yaw),
			0.0, 1.0, 0.0,
			-glm::sin(this->yaw), 0.0, glm::cos(this->yaw) 
		);
	}
};

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(Camera& camera, Light& light);
void Draw(vector<Triangle>& triangles, Camera& camera, const Light& light);
bool ClosestIntersection(const Ray& ray, const vector<Triangle>& triangles, Intersection& intersection);
vec3 DirectLight(const Light& light, const Intersection& i, const vector<Triangle>& triangles);

int main( int argc, char* argv[] )
{
	screen = new SDL2Aux(SCREEN_WIDTH, SCREEN_HEIGHT);
	current_time = SDL_GetTicks();	// Set start value for timer.

	Light light = {.position = vec3(0.0f, 0.0f, 0.0), .color = vec3(11.f, 11.f, 11.f)};
	Camera camera = {.position = vec3(0.0, 0.0, -2.0), .yaw = 0.0f};

	vector<Triangle> triangles;
	LoadTestModel(triangles);

	while(!screen->quitEvent())
	{
		Update(camera, light);
		Draw(triangles, camera, light);
	}

	screen->saveBMP("screenshot.bmp");
	return 0;
}

void Update(Camera& camera, Light& light)
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - current_time);
	//cout << "Render time: " << dt << " ms." << endl;
	current_time = t2;
	SDL_PumpEvents();
	const Uint8* keystate = SDL_GetKeyboardState(NULL);

	dt = dt > 1000 ? 1000 : dt;

	// Camera Rotation
	if (keystate[SDL_SCANCODE_LEFT]) 
		camera.yaw += dt / 1000.0f;
	if (keystate[SDL_SCANCODE_RIGHT]) 
		camera.yaw -= dt / 1000.0f;

	mat3 R = camera.getRotation();
	vec3 forward( R[2][0], R[2][1], R[2][2] );
	// Camera Position
	if (keystate[SDL_SCANCODE_UP] )
		camera.position += forward * dt / 1000.0f;
	if (keystate[SDL_SCANCODE_DOWN] )
		camera.position -= forward * dt / 1000.0f;

	if (keystate[SDL_SCANCODE_D])
		light.position.x += dt / 1000.0f;
	if (keystate[SDL_SCANCODE_A])
		light.position.x -= dt / 1000.0f;
	if (keystate[SDL_SCANCODE_E])
		light.position.y += dt / 1000.0f;
	if (keystate[SDL_SCANCODE_Q])
		light.position.y -= dt / 1000.0f;
	if (keystate[SDL_SCANCODE_W])
		light.position.z += dt / 1000.0f;
	if (keystate[SDL_SCANCODE_S])
		light.position.z -= dt / 1000.0f;
}

void Draw(vector<Triangle>& triangles, Camera& camera, const Light& light)
{
	screen->clearPixels();

	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			Intersection intersection;
			Ray ray = {
				.s = camera.position, 
				.d = camera.getRotation() * glm::normalize(vec3(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, SCREEN_WIDTH / 2))
			};
			if (ClosestIntersection(ray, triangles, intersection)) {
				vec3 indirect_light(0.3, 0.3, 0.3);
				vec3 color = DirectLight(light, intersection, triangles);
				color += indirect_light;
				color *= triangles[intersection.triangleIndex].color;
				screen->putPixel(x, y, color);
			} else {
				screen->putPixel(x, y, vec3(0.0, 0.0, 0.0));
			}
		}
	}

	screen->render();
}

bool ClosestIntersection(const Ray& ray, const vector<Triangle>& triangles, Intersection& intersection) {
	float m = std::numeric_limits<float>::max();
	intersection.distance = m;
	for (int i = 0; i < triangles.size(); i++) {
		Triangle triangle = triangles[i];
		vec3 v0 = triangle.v0;
		vec3 v1 = triangle.v1;
		vec3 v2 = triangle.v2;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = ray.s - v0;
		mat3 A( -ray.d, e1, e2 );
		vec3 x = glm::inverse(A) * b;
		
		float t = x.x;
		float u = x.y;
		float v = x.z;
		//cout << "t: " << t << ", u: " << u << ", v: " << v << endl;
		if (0 < t && t < intersection.distance && 
			0 <= u && 0 <= v && u + v <= 1) {
			intersection.distance = t;
			intersection.position = ray.s + ray.d * t;
			intersection.triangleIndex = i;
		}
	}

	return intersection.distance != m;
}

/**
 * Calculate direct light on intersection "i" given lightsource "light".
 */
vec3 DirectLight(const Light& light, const Intersection& i, const vector<Triangle>& triangles) {
	// D = B max(r̂ . n̂ , 0) = (P max (r̂ . n̂ , 0))/4πr^2
	vec3 light_dir = glm::normalize(light.position - i.position);
	vec3 normal = triangles[i.triangleIndex].normal;
	
	// A = 4πr^2;
	float r = glm::length(light_dir);
	float A = 4 * PI * r * r;

	// Check shadows
	Intersection shadow_i;
	// Adjust position to avoid self collision
	float bias = 1e-4;
	vec3 shadow_ray_start = i.position + light_dir * bias;
	const Ray shadow_ray = {.s = shadow_ray_start, .d = vec3(light.position - i.position)};
	ClosestIntersection(shadow_ray, triangles, shadow_i);
	vec3 shadow_factor(1.0, 1.0, 1.0);
	if (shadow_i.distance < 1) {
		shadow_factor = vec3(0., 0., 0.);
	}

	return shadow_factor * light.color * max(glm::dot(light_dir, normal), 0.0f) / A;
}