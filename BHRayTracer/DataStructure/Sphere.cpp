#include "objects.h"
#include "cyVector.h"
#define HIT_FRONT 0
#define HIT_BACK 1
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
		float t = hitSide == HIT_FRONT ? Min(t1, t2) : Max(t1, t2);


		if (t < 0) {
			if (hInfo.z > Max(t1, t2)) return false;
			hInfo.z = Max(t1, t2);
		}
		else {
			if (hInfo.z < Min(t1, t2)) return false;
			hInfo.z = Min(t1, t2);
		}
		hInfo.p = oc + hInfo.z * dir;
		hInfo.N = hInfo.p;// - Vec3f(0, 0, 0);

		hInfo.front = true;
		return true;
	}

	return false;
}
