#include "objects.h"

extern Vec3f dd_x;
extern Vec3f dd_y;

#define RAY_DIFFERENTIAL
#define PASTED_FACES_BIAS 0.0001f
bool Plane::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	// in obj space

	float rayPz = ray.p.z;
	float rayDz = ray.dir.z;
	// this direction is parallel to this plane
	if (rayDz == 0.0f) return false;
	float t = -rayPz / rayDz;
	// hit the opposite face, or not the closest one
	if (t <= 0 || t > hInfo.z) return false;
	// x is the hit point in the unit plane's plane
	Vec3f x = ray.p + t * ray.dir;
	if (x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1) {
		return false;
	}
	Vec3f faceNormal = Vec3f(0, 0, 1);
	// Back face check
	bool _hitFront = ((-ray.dir).Dot(faceNormal) > 0);
	{
		if (!_hitFront && hitSide == HIT_FRONT) // back face
			return false;
		else if (_hitFront && hitSide == HIT_BACK) // front face
			return false;
	}

	// Set hit info
	hInfo.p = x;
	hInfo.N = faceNormal;
	hInfo.z = t;
	hInfo.front = _hitFront;
	
	// Set uv Info
	Vec3f uvw = Vec3f(0, 0, 0);
	uvw.x = (1 + hInfo.p.x) / 2.f;
	uvw.y = (1 + hInfo.p.y) / 2.f;
	hInfo.uvw = uvw;

	// Ray Differential
	Vec3f duvw[2];
	duvw[0] = Vec3f(0, 0, 0);
	duvw[1] = Vec3f(0, 0, 0);

#ifdef RAY_DIFFERENTIAL
	// d_: prefix of dervitive;
// Direction: the direction of original ray
	{
		Vec3f normalized_dir = ray.dir.GetNormalized();
		float scaled_t = (t * ray.dir).Length();
		Vec3f d_Direction_X = (normalized_dir.Dot(normalized_dir) * dd_x - normalized_dir.Dot(dd_x) *	normalized_dir) / pow(normalized_dir.Dot(normalized_dir), 1.5f);
		Vec3f d_Direction_Y = (normalized_dir.Dot(normalized_dir) * dd_y - normalized_dir.Dot(dd_y) *	normalized_dir) / pow(normalized_dir.Dot(normalized_dir), 1.5f);

		float d_t_x = -(0 + scaled_t * d_Direction_X.Dot(hInfo.N) / normalized_dir.Dot(hInfo.N));
		float d_t_y = -(0 + scaled_t * d_Direction_Y.Dot(hInfo.N) / normalized_dir.Dot(hInfo.N));

		// delta hit point on plane
		Vec3f d_HitPont_x = 0 + scaled_t * d_Direction_X + d_t_x * normalized_dir;
		Vec3f d_HitPont_y = 0 + scaled_t * d_Direction_Y + d_t_y * normalized_dir;

		duvw[0] = d_HitPont_x / 2.f;
		duvw[1] = d_HitPont_y / 2.f;
	}
#endif // RAY_DIFFERENTIAL


	hInfo.duvw[0] = duvw[0];
	hInfo.duvw[1] = duvw[1];

	return true;
}
