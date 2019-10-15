#include "objects.h"

extern Vec3f dd_x;
extern Vec3f dd_y;


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

	// Set hit info
	hInfo.p = x;
	hInfo.N = Vec3f(0, 0, 1);
	hInfo.front = true;
	hInfo.z = t;
	// Set uv Info
	Vec3f uvw;
	uvw.x = (1 + hInfo.p.x) / 2.f;
	uvw.y = (1 + hInfo.p.y) / 2.f;
	hInfo.uvw = uvw;
	// Ray Differential
	Vec3f duvw[2];
	duvw[0] = Vec3f(0, 0, 0);
	duvw[1] = Vec3f(0, 0, 0);

	// dx = dd
	{
		Vec3f d = ray.dir.GetNormalized();
		float _t = (t * ray.dir).Length();
		Vec3f dDx = (d.Dot(d) * dd_x - d.Dot(dd_x) *	d) / pow(d.Dot(d), 1.5f);
		Vec3f dDy = (d.Dot(d) * dd_y - d.Dot(dd_y) *	d) / pow(d.Dot(d), 1.5f);

		float dtx = -(0 + _t * dDx.Dot(hInfo.N) / d.Dot(hInfo.N));
		float dty = -(0 + _t * dDy.Dot(hInfo.N) / d.Dot(hInfo.N));
																					  
		// delta hit point on plane
		Vec3f dXx = 0 +_t* dDx + dtx * d;
		Vec3f dXy = 0 +_t* dDy + dty * d;

		duvw[0] = dXx / 2.f;
		duvw[1] = dXy / 2.f;
	}

	hInfo.duvw[0] = duvw[0];
	hInfo.duvw[1] = duvw[1];

	return true;
}
