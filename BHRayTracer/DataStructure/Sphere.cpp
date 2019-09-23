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
		float t = BIGFLOAT;
		if (t1 < 0 && t2 < 0) {
			return false;
		}
		else if (t1 *t2 <= 0) {
			t = cy::Max(t1, t2);
		}
		else if (t1 > 0 && t2 > 0) {
			t = Min(t1, t2);
		}

		if (hInfo.z < t) return false;

		hInfo.z = t;

		hInfo.p = oc + hInfo.z * dir;
		hInfo.N = hInfo.p; // - Vec3f(0, 0, 0);

		hInfo.front = true;
		return true;
	}

	return false;
}
