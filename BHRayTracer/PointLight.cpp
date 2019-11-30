#include "lights.h"

Vec3f GetSampleAlongNormal(const Vec3f& N, float R);
float Rnd01();
#define PI 3.14159265
#define One_Over_FourPI 0.0795775
cy::Color PointLight::Illuminate(Vec3f const &p, Vec3f const &N) const
{
	Vec3f centerDir = position - p;
	float r = centerDir.Length();
	float rr = r * r;
	if (size > 0) {
		return Shadow(Ray(p, GetSampleAlongNormal(centerDir, size)), 1) *intensity / rr;
	}
	else
		return Shadow(Ray(p, centerDir), 1) * intensity / rr;
}

Ray PointLight::RandomPhoton() const
{
	Vec3f axisZ = Vec3f(0, 0, 1);
	Vec3f axisX = Vec3f(1, 0, 0);
	Vec3f axisY = Vec3f(0, 1, 0);

	float phi = Rnd01() * 2 * PI;
	float theta = ACosSafe(1 - 2 * Rnd01());
	Ray ray = Ray();

	ray.dir = sin(theta) * (axisX * cos(phi) + axisY * sin(phi)) + axisZ * cos(theta);
	ray.p = position;

	return ray;
}
