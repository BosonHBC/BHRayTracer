#include "objects.h"
#include "cyVector.h"
#include <math.h>
#define PI 3.14159265

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
		// if this hit is longer than current closet hit, return false
		if (hInfo.z < t) return false;
		
		// Set hit info
		hInfo.z = t;
		hInfo.p = oc + hInfo.z * dir;
		hInfo.N = hInfo.p; // - Vec3f(0, 0, 0);
		hInfo.front = true;

		// Set uv Info
		Vec3f uvw;
		Vec3f d = hInfo.N.GetNormalized();
		uvw.x = 0.5f + atan2(d.y, d.x) / (2 * PI);
		uvw.y = 0.5f - asin(d.z) / (PI);
		hInfo.uvw = uvw;

		// Ray differential
		Vec3f duvw[2];
		duvw[0] = Vec3f(0, 0, 0);
		duvw[1] = Vec3f(0, 0, 0);
		hInfo.duvw[0] = duvw[0];
		hInfo.duvw[1] = duvw[1];

		return true;
	}

	return false;
}
