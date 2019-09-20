#include "objects.h"

bool Plane::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	// in obj space
	float rayPz = ray.p.z;
	float rayDz = ray.dir.z;

	return false;
}
