#include "objects.h"
#define USE_BVH
#ifdef USE_BVH
//#define SHOW_BOUNDINGBOX
#define  USE_BVH_SHADOW
#endif // USE_BVH

#ifndef Bias
#define  Bias 0.001f
#endif
#ifndef PI
#define PI 3.14159265
#endif

bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
#ifdef USE_BVH
	auto rootNodeId = bvh.GetRootNodeID();
	Box rootBox = Box(bvh.GetNodeBounds(rootNodeId));
	float tmin = -1;
	if (rootBox.IntersectRay(ray, hInfo.z, tmin)) {

		return TraceBVHNode(ray, hInfo, hitSide, rootNodeId);
	}
	return false;
#else
	// Do not use BVH
	bool bHit = false;
	for (unsigned int i = 0; i < nf; i++)
	{
		if (IntersectTriangle(ray, hInfo, hitSide, i)) {
			bHit = true;
		}
	}
	return bHit;
#endif // USE_BVH
}

bool TriObj::ShadowRecursive(Ray const&ray, float t_max)
{
#ifdef USE_BVH_SHADOW
	auto rootNodeId = bvh.GetRootNodeID();
	Box rootBox = Box(bvh.GetNodeBounds(rootNodeId));
	float tmin = -1;
	float t_min = BIGFLOAT;
	bool bHit = false;
	if (rootBox.IntersectRay(ray, BIGFLOAT, tmin)) {
		TraceBVHShadow(ray, t_min, bHit, rootNodeId);
	}
	if(bHit && t_min > Bias && t_min < t_max)
		return true;
	return false;
#else
	for (int i = 0; i < NF(); i++)
	{
		HitInfo hinfo;
		if (IntersectTriangle(ray, hInfo, 0, i)) {
			if (hInfo.z < t_max && hInfo.z > Bias)
				return true;
		}
	}
	return false;
#endif
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
	// Set hit info
	hInfo.z = t;
	hInfo.N = hitPointNormal;
	hInfo.p = vX;
	hInfo.front = true;
	// Set uv info
	Vec3f uvw;
	Vec3f tc = GetTexCoord(faceID, bc);
	hInfo.uvw = tc;

	return true;
}


bool TriObj::TraceBVHNode(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID) const
{
	// if it is a leaf, intersect triangles
	if (bvh.IsLeafNode(nodeID)) {
		// if it is a leaf
		auto elements = bvh.GetNodeElements(nodeID);
		auto elementCount = bvh.GetNodeElementCount(nodeID);
		bool bHit = false;

#ifdef SHOW_BOUNDINGBOX
		bHit = GetBoundBox().IntersectRayWithHitInfo(ray, hInfo);
#else
		for (size_t i = 0; i < elementCount; ++i)
		{
			if (IntersectTriangle(ray, hInfo, hitSide, elements[i]))
			{
				bHit = true;
			}
		}
#endif
		return bHit;
	}
	else {
		unsigned int child1Id = 0;
		unsigned int child2Id = 0;
		bvh.GetChildNodes(nodeID, child1Id, child2Id);
		Box child1 = Box(bvh.GetNodeBounds(child1Id));
		Box child2 = Box(bvh.GetNodeBounds(child2Id));
		float tmin1 = BIGFLOAT;
		float tmin2 = BIGFLOAT;
		bool bChild1Hit = child1.IntersectRay(ray, hInfo.z, tmin1);
		bool bChild2Hit = child2.IntersectRay(ray, hInfo.z, tmin2);

		if (!bChild1Hit && !bChild2Hit) return false;

		if (tmin1 < tmin2) {
			// Trace left child first
			if (TraceBVHNode(ray, hInfo, hitSide, child1Id)) {
				// Store the left hit info for tracing right child
				HitInfo tempHitInfo = hInfo;
				// It trace something on left child
				// If the hit result is even smaller than the right child's bounding box, do not need to trace the right child
				// but in the following case, hit point is inside right child bounding box, needs to trace right child
				if (tempHitInfo.z > tmin2)
				{
					if (TraceBVHNode(ray, tempHitInfo, hitSide, child2Id)) {
						// Right child gets a closer result, use right child's hit info
						hInfo = tempHitInfo;
					}
				}
				// since left child already trace something, so it has to return true
				return true;
			}
			else {
				// trace right child
				return TraceBVHNode(ray, hInfo, hitSide, child2Id);
			}
		}
		else {
			// Trace right child first
			if (TraceBVHNode(ray, hInfo, hitSide, child2Id)) {
				HitInfo tempHitInfo = hInfo;

				if (tempHitInfo.z > tmin1)
				{
					if (TraceBVHNode(ray, tempHitInfo, hitSide, child1Id)) {
						hInfo = tempHitInfo;
					}
				}
				return true;
			}
			else {
				// trace right child
				return TraceBVHNode(ray, hInfo, hitSide, child1Id);
			}
		}
	}
	return false;
}

bool TriObj::TraceBVHShadow(const Ray &ray, float& t_min, bool& hitOnce, unsigned int nodeID) const
{
	if (hitOnce) return hitOnce;
	if (bvh.IsLeafNode(nodeID)) {
		// if it is a leaf
		auto elements = bvh.GetNodeElements(nodeID);
		auto elementCount = bvh.GetNodeElementCount(nodeID);

		HitInfo hInfo;
		for (size_t i = 0; i < elementCount; ++i)
		{
			if (IntersectTriangle(ray, hInfo, 0, elements[i]))
			{
				hitOnce = true;
				t_min = hInfo.z;
			}
		}

		return hitOnce;
	}
	else {
		unsigned int child1Id = 0;
		unsigned int child2Id = 0;
		bvh.GetChildNodes(nodeID, child1Id, child2Id);
		Box child1 = Box(bvh.GetNodeBounds(child1Id));
		Box child2 = Box(bvh.GetNodeBounds(child2Id));
		float tmin1 = BIGFLOAT;
		float tmin2 = BIGFLOAT;
		bool bChild1Hit = child1.IntersectRay(ray, BIGFLOAT, tmin1);
		bool bChild2Hit = child2.IntersectRay(ray, BIGFLOAT, tmin2);

		if (!bChild1Hit && !bChild2Hit) return false;
		return TraceBVHShadow(ray, t_min, hitOnce, child1Id) | TraceBVHShadow(ray, t_min, hitOnce, child2Id);
	}
	return false;
}
