#include "lights.h"
#define GlossyShadow

#ifdef GlossyShadow
#define InitialSampleCount 1
#define MaxSampleCount 1
#endif // GlossyShadow
Vec3f GetSampleAlongNormal(const Vec3f& N, float R);
float Rnd01();
#define PI 3.14159265

cy::Color PointLight::Illuminate(Vec3f const &p, Vec3f const &N) const
{
	if (size > 0) {
		// Adaptive Sample
		int shadowRayHitCount = 0;
		Vec3f centerDir = position - p;
		for (int i = 0; i < InitialSampleCount; ++i)
		{
			shadowRayHitCount += Shadow(Ray(p, GetSampleAlongNormal(centerDir, size)), 1);
		}
		if (shadowRayHitCount == 0) return Color(0, 0, 0);
		else if (shadowRayHitCount == InitialSampleCount) return intensity;
		else {
			// need more samples
			for (int i = 0; i < MaxSampleCount - InitialSampleCount; ++i)
			{
				shadowRayHitCount += Shadow(Ray(p, GetSampleAlongNormal(centerDir, size)), 1);
			}
			float percentage = (float)shadowRayHitCount / MaxSampleCount;
			return percentage * intensity;
		}
	}
	else
	return Shadow(Ray(p, position - p), 1) * intensity;
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
