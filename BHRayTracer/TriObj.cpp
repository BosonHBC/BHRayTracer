#include "objects.h"
#define USE_BVH
#ifdef USE_BVH
//#define SHOW_BOUNDINGBOX
#endif // USE_BVH

bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
#ifdef USE_BVH
	auto rootNodeId = bvh.GetRootNodeID();
	Box rootBox = Box(bvh.GetNodeBounds(rootNodeId));
	float tmin = -1;
	if (rootBox.IntersectRay(ray, BIGFLOAT, tmin)) {
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

		if (tmin1 < tmin2) {
			if (bChild1Hit) {
				if (TraceBVHNode(ray, hInfo, hitSide, child1Id)) return true;
				if (bChild2Hit) {
					return TraceBVHNode(ray, hInfo, hitSide, child2Id);
				}
			}
		}
		else if (tmin1 > tmin2) {
			if (bChild2Hit) {
				if(TraceBVHNode(ray, hInfo, hitSide, child2Id)) return true;
				if (bChild1Hit) {
					return TraceBVHNode(ray, hInfo, hitSide, child1Id);
				}
			}
		}
		else {
			if (bvh.IsLeafNode(bvh.GetFirstChildNode(nodeID)) && bvh.IsLeafNode(bvh.GetSecondChildNode(nodeID))) {
				auto firstElements = bvh.GetNodeElements(child1Id);
				auto secondElements = bvh.GetNodeElements(child2Id);
			}
			else {
				if (bChild2Hit) {
					if (TraceBVHNode(ray, hInfo, hitSide, child2Id)) return true;
					if (bChild1Hit) {
						return TraceBVHNode(ray, hInfo, hitSide, child1Id);
					}
				}
			}
		}

	}
	return false;
}
