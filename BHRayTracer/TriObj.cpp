#include "objects.h"

bool TriObj::IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide /*= HIT_FRONT*/) const
{
	return false;
}


bool TriObj::IntersectTriangle(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const
{
	return false;
}