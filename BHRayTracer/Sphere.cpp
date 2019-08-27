#include "objects.h"
#include "cyVector.h"
bool Sphere::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	// dot(p-c, p-c) = (x-cx)^2 + (y-cy)^2 + (z-cz)^2 = r*r
	// p = origin + dist * dir;
	Vec3f dir = ray.dir;
	dir.Normalize();

	Vec3f oc = ray.p - c;

	float A = dir.Dot(dir);
	float B = 2 * dir.Dot(oc);
	float C = oc.Dot(oc) - r * r;

	float DD = B * B - 4 * A*C;
	if (DD > 0) {
		float t1 = (-B + sqrt(DD)) / (2 * A);
		float t2 = (-B - sqrt(DD)) / (2 * A);
		hInfo.z = t2;
		hInfo.front = true;

		return true;
	}
	
	return false;
}