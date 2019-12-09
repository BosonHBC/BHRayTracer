#include "objects.h"
#include "cyVector.h"
#include <math.h>
#define PI 3.14159265
// cos(85 degree) = 0.08715574274
#define PerpendicularFaceDeterminance 0.087f

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
		// A must > 0 because dir.Dot(dir) > 0
		// so t1 is bigger than t2 if DD > 0
		float t1 = (-B + sqrt(DD)) / (2 * A);
		float t2 = (-B - sqrt(DD)) / (2 * A);
		float t = BIGFLOAT;
		bool _hitFront = true;
		if (t1 < 0 && t2 < 0) {
			// from the ray direction, there is no intersection, all intersections happen on the opposite direction of the sphere
			return false;
		}
		else if (t1 *t2 <= 0) {
			// origin is in the center of the sphere, the intersection is back face, so if the hitSide requires front face, it should not be counted as intersection
			if (hitSide == HIT_FRONT) return false;
			t = t1;
			_hitFront = false;
		}
		else if (t1 > 0 && t2 > 0) {
			// origin is outside the sphere
			if (hitSide == HIT_FRONT || hitSide == HIT_FRONT_AND_BACK) {
				t = t2;
				_hitFront = true;
			}
			else if (hitSide == HIT_BACK) {
				t = t1;
				_hitFront = false;
			}
		}
		// if this hit is longer than current closet hit, return false
		if (hInfo.z < t || t <= 0) return false;

		// Set hit info
		hInfo.z = t;
		hInfo.p = oc + hInfo.z * dir;
		hInfo.N = hInfo.p; // - Vec3f(0, 0, 0);
		hInfo.front = _hitFront;

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
