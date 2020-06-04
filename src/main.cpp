#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <SDL.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

int screen_width = 1280;
int screen_height = 720;
bool running = true;

SDL_Window* window = nullptr;

void app_init();
void app_update(float delta);
void app_render();
void app_release();

int main(int argc, char** argv) {

	SDL_Init(SDL_INIT_EVERYTHING);

	window = SDL_CreateWindow(
		"CPU Raytracer", 
		SDL_WINDOWPOS_UNDEFINED, 
		SDL_WINDOWPOS_UNDEFINED, 
		screen_width, 
		screen_height, 
		SDL_WINDOW_SHOWN);

	SDL_Event e;

	
	app_init();

	uint32_t pre = SDL_GetTicks();
	uint32_t curr = 0;
	float delta = 0.0f;

	// create images
	while (running) {
		curr = SDL_GetTicks();
		delta = (curr - pre) / 1000.0f;
		pre = curr;

		while (SDL_PollEvent(&e)) {
			if (e.type == SDL_QUIT) {
				running = false;
			}
		}

		app_update(delta);
		app_render();

		SDL_UpdateWindowSurface(window);
	}


	app_release();

	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}


struct Sphere {
	glm::vec3 position;
	float radius;
	glm::vec3 color;
	float shiny;
	float reflective;
	float refractive;

	Sphere(glm::vec3 p, float r, glm::vec3 color, float shiny, float reflective, float refractive) {
		this->position = p;
		this->radius = r;
		this->color = color;
		this->shiny = shiny;
		this->reflective = reflective;
		this->refractive = refractive;
	}
};

struct Ray {
	glm::vec3 position;
	glm::vec3 direction;
};

struct Camera {
	glm::vec3 pos;
	glm::vec3 forward;
	glm::vec3 right;
	glm::vec3 up;
	float width;
	float height;

} camera;

struct Light {
	float intensity;
	glm::vec3 color;
	glm::vec3 position;

	Light(float intensity, glm::vec3 color, glm::vec3 position) : intensity(intensity), color(color), position(position) {}

};

void ray_sphere_intersection(Ray* ray, Sphere* sphere, float& t1, float& t2);

void fb_init();
void fb_set(int x, int y, glm::vec3 color);
void fb_clear(glm::vec3 color);
void fb_present();
void fb_release();

void camera_setup();
void camera_update(float delta);

Ray camera_makeRay(glm::vec2 p);

std::vector<Sphere> spheres = {
	Sphere(glm::vec3(-8, 0, 0), 1, glm::vec3(0.5f, 0.5f, 0.5f), 500.0f, 0.2f, 0.0f),
	Sphere(glm::vec3(0, 0, -8), 1, glm::vec3(0.0f, 0.5f, 0.0f), 500.0f, 0.3f, 0.0f),
	Sphere(glm::vec3(8, 0, 0), 1, glm::vec3(0.0f, 0.0f, 0.5f), 500.0f, 0.4f, 0.0f),
	Sphere(glm::vec3(0, 0, 8), 1, glm::vec3(0.5f, 0.0f, 0.0f), 500.0f, 0.5f, 0.0f),
	Sphere(glm::vec3(0, -5001, 0), 5000, glm::vec3(1.0f, 1.0f, 0.0f), 500.0f, 0.5f, 0.0f)
};

std::vector<Light> lights = {
	Light(0.6, glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(0, 1, 0))
};

bool isMultiThreaded = true;

void app_init() {
	fb_init();
	camera_setup();
}

void closestIntersection(Ray* ray, float zmin, float zmax, Sphere** sphere, float& t) {
	for (int i = 0; i < spheres.size(); i++) {
		float t1, t2;

		ray_sphere_intersection(ray, &spheres[i], t1, t2);

		if ((t1 >= zmin && t1 <= zmax) && t1 < t) {
			t = t1;
			*sphere = &spheres[i];
		}

		if ((t2 >= zmin && t2 <= zmax) && t2 < t) {
			t = t2;
			*sphere = &spheres[i];
		}
	}
}

glm::vec3 computeLighting(glm::vec3 P, glm::vec3 N, glm::vec3 V, float shiny, glm::vec3 albedo, float zmin, float zmax) {
	glm::vec3 light = glm::vec3(0.0f);

	for (int i = 0; i < lights.size(); i++) {
		glm::vec3 L = glm::normalize(lights[i].position - P);
		glm::vec3 H = glm::normalize(L + V);
		float t_max = INFINITY;

		Sphere* shadowSphere = nullptr;
		float shadowT = INFINITY;

		Ray shadowRay;
		shadowRay.position = P;
		shadowRay.direction = lights[i].position - P;

		closestIntersection(&shadowRay, 0.001f, t_max, &shadowSphere, shadowT);
		
		if(shadowSphere != nullptr) {
			continue;
		}
		

		// Diffuse Lambersion
		float diffuse = glm::dot(N, L);

		// Specular Blinn-Phong
		float spec = powf(glm::dot(N, H), shiny);

		light += (diffuse * albedo * lights[i].color + spec * lights[i].color) * lights[i].intensity;
	}

	light += glm::vec3(0.1f) * albedo;

	return light;
}

glm::vec3 reflect(glm::vec3 R, glm::vec3 N) {
	return 2.0f * N * glm::dot(N, R) - R;
}

glm::vec3 refract(glm::vec3 I, glm::vec3 N, float eta) {
	float k = 1.0f - eta * eta * (1.0f - glm::dot(N, I) * glm::dot(N, I));

	if (k < 0.0f) {
		return glm::vec3(0.0f);
	}
	else {
		return eta * I - (eta * glm::dot(N, I) + sqrtf(k)) * N;
	}
}

glm::vec3 refract_raytrace(Ray* ray, float zmin, float zmax, float depth) {
	Sphere* sphere = nullptr;
	float t = INFINITY;

	closestIntersection(ray, zmin, zmax, &sphere, t);

	if (sphere == nullptr) {
		return glm::vec3(0.0f);
	}

	glm::vec3 P = ray->position + ray->direction * t;
	glm::vec3 N = P - sphere->position;
	N = glm::normalize(N);

	glm::vec3 localColor = computeLighting(P, N, -ray->direction, sphere->shiny, sphere->color, zmin, zmax);

	// Refract
	float r = sphere->refractive;

	
	if (depth <= 0 || r <= 0) {
		return localColor;
	}

	glm::vec3 R = refract(ray->direction, -N, r);
	Ray refractRay;
	refractRay.direction = R;
	refractRay.position = P;

	return localColor * 0.5f + refract_raytrace(&refractRay, zmin, zmax, depth - 1) * 0.5f;

}

glm::vec3 raytrace(Ray* ray, float zmin, float zmax, float depth) {
	Sphere* s = nullptr;
	float closest_t = INFINITY;

	closestIntersection(ray, zmin, zmax, &s, closest_t);

	if (s == nullptr) {
		return glm::vec3(0.0f);
	}

	glm::vec3 P = ray->position + ray->direction * closest_t;
	glm::vec3 N = P - s->position;
	N = glm::normalize(N);

	glm::vec3 localColor = computeLighting(P, N, -ray->direction, s->shiny, s->color, zmin, zmax);
	
	if (s->refractive > 0) {
		float r = s->refractive;
		glm::vec3 R = refract(ray->direction, -N, r);
		Ray refractRay;
		refractRay.direction = R;
		refractRay.position = P;

		localColor = refract_raytrace(&refractRay, zmin, zmax, depth) * 0.5f + localColor * 0.5f;
	}
	
	// Reflections
	float r = s->reflective;
	if (depth <= 0 || r <= 0) {
		return localColor;
	}

	glm::vec3 R = reflect(-ray->direction, N);

	Ray reflectRay;
	reflectRay.position = P;
	reflectRay.direction = R;

	glm::vec3 reflectedColor = raytrace(&reflectRay, zmin, zmax, depth - 1);

	glm::vec3 color = localColor * (1.0f - r) + reflectedColor * r;

	return color;
}

struct ThreadData {
	int x;
	int y;
	int width;
	int height;
};

struct ThreadHolder {
	int id = 0;
	SDL_Thread* thread = nullptr;
	ThreadData data;
};
int raytracer_thread(void* data) {
	ThreadData* d = (ThreadData*)data;

	for (int y = d->y; y < d->height; y++) {
		for (int x = d->x; x < d->width; x++) {
			glm::vec2 sc = glm::vec2(
				(2.0f * x) / (float)screen_width - 1.0f,
				(2.0f * y) / (float)screen_height - 1.0f
			);

			Ray ray = camera_makeRay(sc);

			glm::vec3 color = raytrace(&ray, 1.0f, 1024.0f, 5);

			fb_set(x, y, color);
		}
	}

	return 0;
}

void app_update(float delta) {
	std::cout << delta << std::endl;
	camera_update(delta);
}

void app_render() {
	fb_clear(glm::vec3(0.0f, 0.0f, 0.0f));

	uint32_t pre = SDL_GetTicks();
	
	if (isMultiThreaded) {
		uint32_t threadX = 16;
		uint32_t threadY = 9;

		uint32_t sizeX = screen_width / threadX;
		uint32_t sizeY = screen_height / threadY;

		std::vector<ThreadHolder> threadHolders(threadX * threadY);

		for (int y = 0; y < threadY; y++) {
			for (int x = 0; x < threadX; x++) {
				threadHolders[y * threadX + x].id = y * threadX + x;

				threadHolders[y * threadX + x].data.x = x * sizeX;
				threadHolders[y * threadX + x].data.width = x * sizeX + sizeX;
				threadHolders[y * threadX + x].data.y = y * sizeY;
				threadHolders[y * threadX + x].data.height = y * sizeY + sizeY;

				std::stringstream ss;

				ss << "RayTraced_Thread_" << threadHolders[y * threadX + x].id << std::endl;

				threadHolders[y * threadX + x].thread = SDL_CreateThread(raytracer_thread, ss.str().c_str(), &threadHolders[y * threadX + x].data);
			}
		}

		for (int i = 0; i < threadHolders.size(); i++) {
			SDL_WaitThread(threadHolders[i].thread, nullptr);
		}

		threadHolders.clear();
	}
	else {

		for (int y = 0; y < screen_height; y++) {
			for (int x = 0; x < screen_width; x++) {

				glm::vec2 sc = glm::vec2(
					(2.0f * x) / (float)screen_width - 1.0f,
					(2.0f * y) / (float)screen_height - 1.0f
				);

				Ray ray = camera_makeRay(sc);

				glm::vec3 color = raytrace(&ray, 1.0f, 1024.0f, 10);

				fb_set(x, y, color);
			}
		}
	}

	float time = (SDL_GetTicks() - pre) / 1000.0f;
	std::cout << time << "S" << std::endl;
	fb_present();
}

void app_release() {
	fb_release();
}


SDL_Surface* screen = nullptr;
SDL_Color* pixels = nullptr;

void fb_init() {
	screen = SDL_CreateRGBSurfaceWithFormat(0, screen_width, screen_height, 32, SDL_PIXELFORMAT_RGBA32);
	pixels = (SDL_Color*)screen->pixels;
}

void fb_set(int x, int y, glm::vec3 color) {

	if (x < 0 || y < 0 || x > screen_width - 1 || y > screen_height - 1) {
		return;
	}

	SDL_Color c;
	c.r = glm::clamp(color.r, 0.0f, 1.0f) * 255;
	c.g = glm::clamp(color.g, 0.0f, 1.0f) * 255;
	c.b = glm::clamp(color.b, 0.0f, 1.0f) * 255;
	c.a = 1.0f * 255;

	pixels[y * screen_width + x] = c;

}

void fb_clear(glm::vec3 clearColor) {
	for (int y = 0; y < screen_height; y++) {
		for (int x = 0; x < screen_width; x++) {
			fb_set(x, y, clearColor);
		}
	}
}

void fb_present() {
	SDL_BlitSurface(screen, nullptr, SDL_GetWindowSurface(window), nullptr);
}

void fb_release() {
	screen = nullptr;
	pixels = nullptr;
}

void ray_sphere_intersection(Ray* ray, Sphere* sphere, float& t1, float& t2) {
	glm::vec3 v = ray->position - sphere->position;

	float k1 = glm::dot(ray->direction, ray->direction);
	float k2 = 2 * glm::dot(v, ray->direction);
	float k3 = glm::dot(v, v) - sphere->radius * sphere->radius;

	float d = k2 * k2 - 4 * k1 * k3;

	if (d < 0.0f) {
		t1 = INFINITY, t2 = INFINITY;
	}

	t1 = (-k2 + sqrt(d)) / (2 * k1);
	t2 = (-k2 - sqrt(d)) / (2 * k1);
}

float yaw = 0.0f;
float pitch = 0.0f;

void camera_setup() {
	float fov = 60.0f;
	float aspect = (float)screen_width / (float)screen_height;

	camera.pos = glm::vec3(0.0f, 0.0f, 0.0f);


	camera.forward = glm::normalize(glm::vec3(0, 0, -1));
	camera.right = glm::normalize(glm::cross(camera.forward, glm::vec3(0.0f, 1.0f, 0.0f)));
	camera.up = glm::cross(camera.forward, camera.right);

	camera.height = tan(fov);
	camera.width = camera.height * aspect;
}

Ray camera_makeRay(glm::vec2 point) {
	glm::vec3 d = camera.forward + point.x * camera.width * camera.right + point.y * camera.height * camera.up;
	Ray ray;
	ray.position = camera.pos;
	ray.direction = glm::normalize(d);
	return ray;
}

void camera_update(float delta) {
	const uint8_t* keys = SDL_GetKeyboardState(nullptr);

	if (keys[SDL_SCANCODE_LEFT]) {
		yaw -= 32.0f * delta;
	}

	if (keys[SDL_SCANCODE_RIGHT]) {
		yaw += 32.0f * delta;
	}

	if (yaw < -360.0f) {
		yaw += 360.0f;
	}

	if (yaw > 360.0f) {
		yaw -= 360.0f;
	}

	if (keys[SDL_SCANCODE_UP]) {
		pitch += 32.0f * delta;
	}

	if (keys[SDL_SCANCODE_DOWN]) {
		pitch -= 32.0f * delta;
	}


	if (pitch < -90.0f) {
		pitch = -90.0f;
	}

	if (pitch > 90.0f) {
		pitch = 90.0f;
	}

	glm::vec3 direction = glm::vec3(
		glm::cos(glm::radians(yaw)) * glm::cos(glm::radians(pitch)),
		glm::sin(glm::radians(pitch)),
		glm::sin(glm::radians(yaw)) * glm::cos(glm::radians(pitch))
	);

	camera.forward = glm::normalize(direction);
	camera.right = glm::normalize(glm::cross(camera.forward, glm::vec3(0.0f, 1.0f, 0.0f)));
	camera.up = glm::cross(camera.forward, camera.right);


	glm::vec3 forward = glm::vec3(
		camera.forward.x,
		0.0f,
		camera.forward.z);

	if (keys[SDL_SCANCODE_W]) {
		camera.pos += 2.5f * delta * forward;
	}

	if (keys[SDL_SCANCODE_S]) {
		camera.pos -= 2.5f * delta * forward;
	}

	if (keys[SDL_SCANCODE_A]) {
		camera.pos -= 2.5f * delta * camera.right;
	}

	if (keys[SDL_SCANCODE_D]) {
		camera.pos += 2.5f * delta * camera.right;
	}

	if (keys[SDL_SCANCODE_SPACE]) {
		camera.pos.y += 2.5f * delta;
	}

	if (keys[SDL_SCANCODE_LSHIFT]) {
		camera.pos.y -= 2.5f * delta;
	}
}