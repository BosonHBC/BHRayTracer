#include "objects.h"

bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	bool bHit;
	for (unsigned int i = 0; i < nf; i++)
	{
		if (IntersectTriangle(ray, hInfo, hitSide, i)) {
			bHit = true;
		}
	}
	return bHit;
}


bool TriObj::IntersectTriangle(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const
{
	TriFace face = F(faceID);
	unsigned int v0i = face.v[0];
	unsigned int v1i = face.v[1];
	unsigned int v2i = face.v[2];

	Vec3f v0 = V(v0i);
	Vec3f v1 = V(v1i);
	Vec3f v2 = V(v2i);

	Vec3f vN = (v1 - v0).Cross(v2 - v0);

	float t;

	t = (vN.Dot(v0) - vN.Dot(ray.p)) / (vN.Dot(ray.dir));
	// on the back side or the distance is larger than the current shortest one
	if (t < 0 || t > hInfo.z) return false;

	// vX is the hit point on the XY plane
	Vec3f vX = ray.p + t * ray.dir;

	Vec3f absVN = vN.Abs();

	Vec2f v0_2d, v1_2d, v2_2d, vX_2d;
	if (absVN.x >= absVN.y && absVN.x >= absVN.z)
	{
		// X plane as projection
		v0_2d.x = v0.y;
		v0_2d.y = v0.z;

		v1_2d.x = v1.y;
		v1_2d.y = v1.z;

		v2_2d.x = v2.y;
		v2_2d.y = v2.z;

		vX_2d.x = vX.y;
		vX_2d.y = vX.z;

	}
	else 	if (absVN.y >= absVN.x && absVN.y >= absVN.z)
	{
		// Y plane as projection
		v0_2d.x = v0.x;
		v0_2d.y = v0.z;

		v1_2d.x = v1.x;
		v1_2d.y = v1.z;

		v2_2d.x = v2.x;
		v2_2d.y = v2.z;

		vX_2d.x = vX.x;
		vX_2d.y = vX.z;
	}
	else 	if (absVN.z >= absVN.y && absVN.z >= absVN.x)
	{
		// Z plane as projection
		v0_2d.x = v0.x;
		v0_2d.y = v0.y;

		v1_2d.x = v1.x;
		v1_2d.y = v1.y;

		v2_2d.x = v2.x;
		v2_2d.y = v2.y;

		vX_2d.x = vX.x;
		vX_2d.y = vX.y;
	}

	float a0_2d = (v1_2d - vX_2d).Cross(v2_2d - vX_2d) / 2.f;

	float a1_2d = (v2_2d - vX_2d).Cross(v0_2d - vX_2d) / 2.f;

	float a2_2d = (v0_2d - vX_2d).Cross(v1_2d - vX_2d) / 2.f;
	if ((a2_2d < 0 || a1_2d < 0 || a0_2d < 0) && !(a0_2d < 0 && a1_2d < 0 && a2_2d < 0))return false;

	// until here, this triangle is intersecting with the ray
	// Can calculate the barycentric coordinates coordinate of this point and return the normal

	float a_2d = a0_2d + a1_2d + a2_2d;
	Vec3f bc = Vec3f(a0_2d / a_2d, a1_2d / a_2d, a2_2d / a_2d);
	Vec3f hitPointNormal = GetNormal(faceID, bc);
	hInfo.z = t;
	hInfo.N = hitPointNormal;
	hInfo.p = vX;
	hInfo.front = true;
	return true;
}


bool TriObj::TraceBVHNode(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID) const
{
	if (bvh.IsLeafNode(nodeID)) {
		

		//IntersectTriangle(ray, hInfo, hitSide, )
	}
	return false;
}
