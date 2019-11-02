#include "lights.h"
cy::Color PointLight::Illuminate(Vec3f const &p, Vec3f const &N) const
{
	return Shadow(Ray(p, position - p), 1) * intensity;
}
