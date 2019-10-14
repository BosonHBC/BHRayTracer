
//-------------------------------------------------------------------------------
///
/// \file       objects.h 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    7.0
/// \date       October 6, 2015
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
	virtual bool IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(-1, -1, -1, 1, 1, 1); }
	virtual void ViewportDisplay(const Material *mtl) const;
};

extern Sphere theSphere;

//-------------------------------------------------------------------------------

class Plane : public Object
{
public:
	virtual bool IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(-1, -1, 0, 1, 1, 0); }
	virtual void ViewportDisplay(const Material *mtl) const;
};

extern Plane thePlane;

//-------------------------------------------------------------------------------

class TriObj : public Object, public cyTriMesh
{
public:
	virtual bool IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide = HIT_FRONT) const;
	virtual Box GetBoundBox() const { return Box(GetBoundMin(), GetBoundMax()); }
	virtual void ViewportDisplay(const Material *mtl) const;

	bool Load(const char *filename, bool loadMtl)
	{
		bvh.Clear();
		if (!LoadFromFileObj(filename, loadMtl)) return false;
		if (!HasNormals()) ComputeNormals();
		ComputeBoundingBox();
		bvh.SetMesh(this, 4);
		return true;
	}
	bool IntersectTriangle(const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int faceID) const;
	bool ShadowRecursive(Ray const&ray, float t_max);
private:
	cyBVHTriMesh bvh;
	bool TraceBVHNode(const Ray &ray, HitInfo &hInfo, int hitSide, unsigned int nodeID) const;
	bool TraceBVHShadow(const Ray &ray, float& t_min,bool& hitOnce, unsigned int nodeID)const;
};



//-------------------------------------------------------------------------------

#endif
