#include "objects.h"

bool Plane::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	// in obj space
	float rayPz = ray.p.z;
	float rayDz = ray.dir.z;

	float t = -rayPz / rayDz;
	// hit the opposite face, or not the closest one
	if (t < 0 || t > hInfo.z) return false;
	// x is the hit point in the unit plane's plane
	Vec3f x = ray.p + t * ray.dir;
	if (x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1) {
		return false;
	}
	hInfo.p = x;
	hInfo.N = Vec3f(0, 0, 1);
	hInfo.front = true;
	hInfo.z = t;
	return true;
}
