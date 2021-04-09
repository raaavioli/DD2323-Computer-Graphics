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

struct Pixel {
    int x;
	int y;
	float z_inv;
	vec3 color;
};

struct Vertex {
	vec3 position;
	vec3 normal;
	vec3 color;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const int FOCAL_LENGTH = SCREEN_WIDTH / 2;
float depth_buffer[SCREEN_HEIGHT][SCREEN_WIDTH];
SDL2Aux* screen;
int t;


vec3 light_pos(0.0, -0.5, -0.7);
vec3 light_power = 1.0f*vec3( 1, 1, 1 );
vec3 indirect_light_power_per_area = 0.3f*vec3( 1, 1, 1 );

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update(Camera& camera);
void Draw(const vector<Triangle>& triangles, const Camera& camera);
void DrawLine (ivec2 a, ivec2 b, vec3 color);
void DrawPolygonLine (const vector<Vertex>& vertices, const Camera& camera);
void DrawRows (const vector<Pixel>& left_pixels, const vector<Pixel>& right_pixels);
void DrawPolygonFill (const vector<Vertex>& vertices, const Camera& camera);
void ComputePolygonRows (const vector<Pixel>& polygon, vector<Pixel>& left_pixels, vector<Pixel>& right_pixels);
void Interpolate (ivec2 a, ivec2 b, vector<ivec2>& result);
void Interpolate (Pixel a, Pixel b, vector<Pixel>& result);

void VertexShader (const Vertex& v, Pixel& p, const Camera& camera);
void FragmentShader (const Pixel& p);

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
	// cout << "Render time: " << dt << " ms." << endl;

	mat3 R = camera.getRotation();
	vec3 right (R[0][0], R[0][1], R[0][2]);
	vec3 forward (R[2][0], R[2][1], R[2][2]);

	const Uint8* keystate = SDL_GetKeyboardState(NULL);

	if (dt > 1000) dt = 1000;
	if(keystate[SDL_SCANCODE_UP]) camera.pitch += dt / 1000.f;
	if(keystate[SDL_SCANCODE_DOWN]) camera.pitch -= dt / 1000.f;
	if(keystate[SDL_SCANCODE_RIGHT]) camera.yaw -= dt / 1000.f;
	if(keystate[SDL_SCANCODE_LEFT]) camera.yaw += dt / 1000.f;
	if(keystate[SDL_SCANCODE_RSHIFT])
		;
	if(keystate[SDL_SCANCODE_RCTRL])
		;

	if(keystate[SDL_SCANCODE_W]) camera.position += forward * dt / 1000.f;
	if(keystate[SDL_SCANCODE_S]) camera.position -= forward * dt / 1000.f;
	if(keystate[SDL_SCANCODE_D]) camera.position += right * dt / 1000.f;
	if(keystate[SDL_SCANCODE_A]) camera.position -= right * dt / 1000.f;
	if(keystate[SDL_SCANCODE_T]) light_pos += forward * dt / 1000.0f;
	if(keystate[SDL_SCANCODE_G]) light_pos -= forward * dt / 1000.0f;
	if(keystate[SDL_SCANCODE_F]) light_pos -= right * dt / 1000.0f;
	if(keystate[SDL_SCANCODE_H]) light_pos += right * dt / 1000.0f;
}

void Draw (const vector<Triangle>& triangles, const Camera& camera) {
	screen->clearPixels();
	// Clear depth buffer
	for (int y = 0; y < SCREEN_HEIGHT; y++)
		for (int x = 0; x < SCREEN_WIDTH; x++)
			depth_buffer[y][x] = 0;
	
	for(int i = 0; i<triangles.size(); i++ ) {
		vector<Vertex> vertices(3);
		Vertex v0 = {.position = triangles[i].v0, .normal = triangles[i].normal, .color = triangles[i].color};
		Vertex v1 = {.position = triangles[i].v1, .normal = triangles[i].normal, .color = triangles[i].color};
		Vertex v2 = {.position = triangles[i].v2, .normal = triangles[i].normal, .color = triangles[i].color};
		vertices[0] = v0;
		vertices[1] = v1;
		vertices[2] = v2;
		DrawPolygonFill (vertices, camera);
	}
	
	screen->render();
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

void DrawPolygonLine (const vector<Vertex>& vertices, const Camera& camera) {
	vector<Pixel> proj_vertices(vertices.size());
	for (int i = 0; i < vertices.size(); i++)
		VertexShader (vertices[i], proj_vertices[i], camera);

	int num_vertices = proj_vertices.size();
	for (int i = 0; i < num_vertices; i++) {
		int j = (i + 1) % num_vertices;
		vec3 color(1.0, 1.0, 1.0);
		DrawLine (
			ivec2(proj_vertices[i].x, proj_vertices[i].y), 
			ivec2(proj_vertices[j].x, proj_vertices[j].y), 
			color
		);
	}
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result ) { 
	vec2 start(a);
	vec2 end(b);
	float len = result.size() - 1;
	for (int i = 0; i <= len; i++) {
		float p = i / len;
		result[i] = glm::round((1 - p) * start + p * end);
	}
}

void Interpolate (Pixel a, Pixel b, vector<Pixel>& result ) { 
	float len = result.size() - 1;
	for (int i = 0; i <= len; i++) {
		float p = i / len;
		Pixel res = {
			.x = static_cast<int>(glm::round((1 - p) * a.x + p * b.x)),
			.y = static_cast<int>(glm::round((1 - p) * a.y + p * b.y)),
			.z_inv = (1 - p) * a.z_inv + p * b.z_inv,
			.color = (1 - p) * a.color + p * b.color,
		};
		result[i] = res;
	}
}

void ComputePolygonRows (const vector<Pixel>& polygon_pixels, 
	vector<Pixel>& left_pixels, vector<Pixel>& right_pixels) {
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

		vector<Pixel> edge_pixels (polygon_rows);
		Interpolate (polygon_pixels[i], polygon_pixels[j], edge_pixels);
		for (int e = 0; e < edge_pixels.size(); e++) {
			int x = edge_pixels[e].x;
			int y = edge_pixels[e].y;
			if (x < left_pixels[y - min_y].x)
				left_pixels[y - min_y] = edge_pixels[e];
			if (x > right_pixels[y - min_y].x)
				right_pixels[y - min_y] = edge_pixels[e];
		}
	}
}

void DrawRows (const vector<Pixel>& left_pixels, const vector<Pixel>& right_pixels) {
	for (int row = 0; row < left_pixels.size(); row++) {
		Pixel left = left_pixels[row];
		Pixel right = right_pixels[row];
		vector<Pixel> row_pixels(right.x - left.x + 1);
		Interpolate (left, right, row_pixels);
		for (int i = 0; i < row_pixels.size(); i++) {
			FragmentShader (row_pixels[i]);
		}
	}
}

void DrawPolygonFill (const vector<Vertex>& vertices, const Camera& camera) {
	int V = vertices.size();
	vector<Pixel> vertexPixels( V );
	for (int i = 0; i < V; i++)
		VertexShader(vertices[i], vertexPixels[i], camera);
	vector<Pixel> left_pixels;
	vector<Pixel> right_pixels;
	ComputePolygonRows( vertexPixels, left_pixels, right_pixels );
	DrawRows (left_pixels, right_pixels);
}

void VertexShader (const Vertex& v, Pixel& p, const Camera& camera) {
	vec3 position(v.position);
	vec3 light_position(light_pos);
	vec3 light_dir = glm::normalize(light_position - position);
	const mat3 rotation = camera.getRotation ();
	position -= camera.position;
	position = position * rotation;
	light_position -= camera.position;
	p.z_inv = 1 / position.z;
	p.x = FOCAL_LENGTH * position.x / position.z + SCREEN_WIDTH / 2;
	p.y = FOCAL_LENGTH * position.y / position.z + SCREEN_HEIGHT / 2;
	
	vec3 normal = glm::normalize(v.normal);
	float cosv = glm::max(glm::dot(light_dir, normal), 0.0f);
	p.color = v.color * (cosv * light_power + indirect_light_power_per_area);
}

void FragmentShader (const Pixel& p) {
	int x = p.x;
	int y = p.y;
	if(x >= 0 && x < SCREEN_WIDTH && y >= 0 && y < SCREEN_HEIGHT && p.z_inv > depth_buffer[y][x]) {
		depth_buffer[y][x] = p.z_inv;
		screen->putPixel(x, y, p.color); 
	}
}