
//-------------------------------------------------------------------------------
///
/// \file       objects.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    6.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifndef _OBJECTS_H_INCLUDED_
#define _OBJECTS_H_INCLUDED_

#include "scene.h"
#include "cyTriMesh.h"
#include "cyBVH.h"

//-------------------------------------------------------------------------------

class Sphere : public Object
{
public:
	virtual bool IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(-1, -1, -1, 1, 1, 1); }
	virtual void ViewportDisplay(const Material *mtl) const;
};

extern Sphere theSphere;

//-------------------------------------------------------------------------------

class Plane : public Object
{
public:
	virtual bool IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(-1, -1, 0, 1, 1, 0); }
	virtual void ViewportDisplay(const Material *mtl) const;
};

extern Plane thePlane;

//-------------------------------------------------------------------------------

class TriObj : public Object, public cyTriMesh
{
public:
	virtual bool IntersectRay(Ray const &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(GetBoundMin(), GetBoundMax()); }
	virtual void ViewportDisplay(const Material *mtl) const;

	bool Load(char const *filename)
	{
		bvh.Clear();
		if (!LoadFromFileObj(filename)) return false;
		if (!HasNormals()) ComputeNormals();
		ComputeBoundingBox();
		bvh.SetMesh(this, 4);
		return true;
	}
	bool IntersectTriangle(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const;
	bool ShadowRecursive(Ray const&ray, HitInfo &hInfo, float t_max);
private:
	cyBVHTriMesh bvh;
	bool TraceBVHNode(Ray const &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID) const;
};

//-------------------------------------------------------------------------------

#endif