#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDL2Auxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::vec2;
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
void Draw(const vector<Triangle>& triangles, const Camera& camera);
void DrawLine (ivec2 a, ivec2 b, vec3 color);
void DrawPolygonLine (const vector<vec3>& vertices, const Camera& camera);
void DrawRows (const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels, const vec3& color);
void DrawPolygonFill (const vector<vec3>& vertices, const Camera& camera, const vec3& color);
void ComputePolygonRows (const vector<ivec2>& polygon, vector<ivec2>& left_pixels, vector<ivec2>& right_pixels);
void Interpolate (ivec2 a, ivec2 b, vector<ivec2>& result);
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
	//cout << "Render time: " << dt << " ms." << endl;

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

void Draw (const vector<Triangle>& triangles, const Camera& camera) {
	screen->clearPixels();
	
	cout << "pre draw" << endl;
	for(int i = 0; i<triangles.size(); i++ ) {
		if (i > 6) continue;
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		DrawPolygonFill (vertices, camera, triangles[i].color);
	}
	
	cout << "before render" << endl;
	screen->render();
	cout << "after render" << endl;
}

void DrawLine (ivec2 a, ivec2 b, vec3 color) {
	ivec2 delta = glm::abs (a - b);
	int pixels = glm::max (delta.x, delta.y) + 1;
	vector<ivec2> points (pixels);
	Interpolate (a, b, points);
	for (int i = 0; i < points.size(); i++) {
		screen->putPixel (points[i].x, points[i].y, color);
	}
}

void DrawPolygonLine (const vector<vec3>& vertices, const Camera& camera) {
	vector<ivec2> proj_vertices(vertices.size());
	for (int i = 0; i < vertices.size(); i++)
		VertexShader (vertices[i], proj_vertices[i], camera);

	int num_vertices = proj_vertices.size();
	for (int i = 0; i < num_vertices; i++) {
		int j = (i + 1) % num_vertices;
		vec3 color(1.0, 1.0, 1.0);
		DrawLine (proj_vertices[i], proj_vertices[j], color);
	}
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ) { 
	int N = result.size();
	vec2 step = vec2(b - a) / float (max (N-1,1));
	vec2 current (a);
	for( int i=0; i<N; ++i ) {
		result[i] = current;
		current += step;
	}
}

void ComputePolygonRows (const vector<ivec2>& polygon_pixels, 
	vector<ivec2>& left_pixels, vector<ivec2>& right_pixels) {
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int min_y = numeric_limits<int>::max();
	int max_y = -numeric_limits<int>::max();
	for (int i = 0; i < polygon_pixels.size (); i++) {
		if (polygon_pixels[i].y < min_y)
			min_y = polygon_pixels[i].y;
		if (polygon_pixels[i].y > max_y)
			max_y = polygon_pixels[i].y;
	}
	int polygon_rows = max_y - min_y + 1;

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	left_pixels.resize(polygon_rows);
	right_pixels.resize(polygon_rows);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < polygon_rows; i++) { 
			left_pixels[i].x = numeric_limits<int>::max();
			right_pixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i < polygon_pixels.size (); i++) {
		int j = (i + 1) % polygon_pixels.size ();

		vector<ivec2> edge_pixels (polygon_rows);
		Interpolate (polygon_pixels[i], polygon_pixels[j], edge_pixels);
		for (int e = 0; e < edge_pixels.size(); e++) {
			int x = edge_pixels[e].x;
			int y = edge_pixels[e].y;
			if (x < left_pixels[y - min_y].x) {
				left_pixels[y - min_y].x = x;
				left_pixels[y - min_y].y = y;
			}
			if (x > right_pixels[y - min_y].x) {
				right_pixels[y - min_y].x = x;
				right_pixels[y - min_y].y = y;
			}
		}
	}
}

void DrawRows (const vector<ivec2>& left_pixels, const vector<ivec2>& right_pixels, const vec3& color) {
	for (int i = 0; i < left_pixels.size(); i++) {
		for (int x = left_pixels[i].x; x <= right_pixels[i].x; x++) {
			int y = left_pixels[i].y;
			if (x >= 0 && x < SCREEN_WIDTH && y >= 0 && y < SCREEN_HEIGHT)
				screen->putPixel (x, left_pixels[i].y, color);
		}
	}
}

void DrawPolygonFill (const vector<vec3>& vertices, const Camera& camera, const vec3& color) {
	int V = vertices.size();
	vector<ivec2> vertexPixels( V );
	for (int i = 0; i < V; i++)
		VertexShader(vertices[i], vertexPixels[i], camera);
	vector<ivec2> left_pixels;
	vector<ivec2> right_pixels;
	ComputePolygonRows( vertexPixels, left_pixels, right_pixels );
	DrawRows (left_pixels, right_pixels, color);
}

void VertexShader (const vec3& v, ivec2& p, const Camera& camera) {
	vec3 vertex(v);
	vertex -= camera.position;
	vertex = vertex * camera.getRotation();
	p.x = FOCAL_LENGTH * vertex.x / vertex.z + SCREEN_WIDTH / 2;
	p.y = FOCAL_LENGTH * vertex.y / vertex.z + SCREEN_HEIGHT / 2;
}