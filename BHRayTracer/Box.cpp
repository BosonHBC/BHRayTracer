#include "scene.h"

bool Box::IntersectRay(Ray const &r, float t_max, float& t_min) const
{
	// Reference points in box
	Vec3f v0 = Box::Corner(0); // in XY plane1 and ZX plane1 and YZ plane1
	Vec3f v7 = Box::Corner(7); // in  XY plane2 and ZX plane2 and YZ plane2

	// for XY plane
	Vec3f vXYN = Vec3f(0, 0, 1); // Up
	float nDotrr = vXYN.Dot(r.dir);
	float nDotrp = vXYN.Dot(r.p);
	float tz1 = (vXYN.Dot(v0) - nDotrp) / nDotrr; // XY plane1
	float tz2 = (vXYN.Dot(v7) - nDotrp) / nDotrr; // XY plane2

	// for ZX plane
	Vec3f vZXN = Vec3f(0, 1, 0); // Right
	nDotrr = vZXN.Dot(r.dir);
	nDotrp = vZXN.Dot(r.p);
	float ty1 = (vZXN.Dot(v0) - nDotrp) / nDotrr; // ZX plane1
	float ty2 = (vZXN.Dot(v7) - nDotrp) / nDotrr; // ZX plane2

	// for YZ pane
	Vec3f vYZN = Vec3f(1, 0, 0);  // forward
	nDotrr = vYZN.Dot(r.dir);
	nDotrp = vYZN.Dot(r.p);
	float tx1 = (vYZN.Dot(v0) - nDotrp) / nDotrr; // YZ plane1
	float tx2 = (vYZN.Dot(v7) - nDotrp) / nDotrr; // YZ plane2

	float tMin = Max(
		Max(
			Min(tx1, tx2), Min(ty1, ty2)
		), Min(tz1, tz2)
	);
	float tMax = Min(
		Min(
			Max(tx1, tx2), Max(ty1, ty2)
		), Max(tz1, tz2)
	);

	if (tMin <= tMax && tMin < t_max) {
		t_min = tMin;
		return true;
	}
	return false;
}

bool Box::IntersectRayWithHitInfo(Ray const &r, HitInfo& hInfo) const
{
	// Reference points in box
	Vec3f v0 = Box::Corner(0); // in XY plane1 and ZX plane1 and YZ plane1
	Vec3f v7 = Box::Corner(7); // in  XY plane2 and ZX plane2 and YZ plane2
	// for XY plane
	Vec3f vXYN = Vec3f(0, 0, 1); // Up
	float nDotrr = vXYN.Dot(r.dir);
	float nDotrp = vXYN.Dot(r.p);
	float tz1 = (vXYN.Dot(v0) - nDotrp) / nDotrr; // XY plane1
	float tz2 = (vXYN.Dot(v7) - nDotrp) / nDotrr; // XY plane2

	// for ZX plane
	Vec3f vZXN = Vec3f(0, 1, 0); // Right
	nDotrr = vZXN.Dot(r.dir);
	nDotrp = vZXN.Dot(r.p);
	float ty1 = (vZXN.Dot(v0) - nDotrp) / nDotrr; // ZX plane1
	float ty2 = (vZXN.Dot(v7) - nDotrp) / nDotrr; // ZX plane2

	// for YZ pane
	Vec3f vYZN = Vec3f(1, 0, 0);  // forward
	nDotrr = vYZN.Dot(r.dir);
	nDotrp = vYZN.Dot(r.p);
	float tx1 = (vYZN.Dot(v0) - nDotrp) / nDotrr; // YZ plane1
	float tx2 = (vYZN.Dot(v7) - nDotrp) / nDotrr; // YZ plane2

	float tMin = Max(
		Max(
			Min(tx1, tx2), Min(ty1, ty2)
		), Min(tz1, tz2)
	);
	float tMax = Min(
		Min(
			Max(tx1, tx2), Max(ty1, ty2)
		), Max(tz1, tz2)
	);

	if (tMin <= tMax) {
		hInfo.z = tMin;
		hInfo.p = r.p + tMin * r.dir;
		hInfo.N = hInfo.p;
		return true;
	}
	return false;
}
