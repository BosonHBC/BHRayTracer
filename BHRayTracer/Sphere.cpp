#include "objects.h"
#include "cyVector.h"
bool Sphere::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	// dot(p, p) = (x)^2 + (y)^2 + (z)^2 = 1
	// p = origin + dist * dir;
	Vec3f dir = ray.dir;
	Vec3f oc = ray.p;

	float A = dir.Dot(dir);
	float B = 2 * dir.Dot(oc);
	float C = oc.Dot(oc) - 1;

	float DD = B * B - 4 * A*C;
	if (DD > 0) {
		float t1 = (-B + sqrt(DD)) / (2 * A);
		float t2 = (-B - sqrt(DD)) / (2 * A);
		hInfo.z = Min(t1, t2);
		hInfo.front = true;

		return true;
	}
	
	return false;
}