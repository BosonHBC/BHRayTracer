#include "lights.h"
#include "scene.h"
#include "objects.h"

#define ShadowBias 0.00001f

extern Node rootNode;
bool ShadowRayRecursive(Node* root, const Ray& ray, float t_max);

float GenLight::Shadow(Ray ray, float t_max /*= BIGFLOAT*/)
{
	return ShadowRayRecursive(&rootNode, ray, t_max) ? 0.f : 1.f;
}

bool ShadowRayRecursive(Node* root, const Ray& ray, float t_max) {

	Ray transformedRay = root->ToNodeCoords(ray);
	for (int i = 0; i < root->GetNumChild(); i++)
	{
		if (ShadowRayRecursive(root->GetChild(i), transformedRay, t_max)) return true;
	}
	if (root->GetNodeObj() != nullptr) {
		// it is a sphere
		if (Sphere* ptr = dynamic_cast<Sphere*>(root->GetNodeObj())) {
			// transform ray to child coordinate
			Vec3f dir = transformedRay.dir;
			Vec3f oc = transformedRay.p;
			float A = dir.Dot(dir);
			float B = 2 * dir.Dot(oc);
			float C = oc.Dot(oc) - 1;

			float DD = B * B - 4 * A*C;
			if (DD > 0) {
				float t1 = (-B + sqrt(DD)) / (2 * A);
				float t2 = (-B - sqrt(DD)) / (2 * A);
				float t = Min(t1, t2);
				if (t < 0) return false;
				if (t < t_max && t > ShadowBias)
				{
					return true;
				}
			}
		}
		else if (Plane* ptr = dynamic_cast<Plane*>(root->GetNodeObj()))
		{
			// in obj space
			float rayPz = transformedRay.p.z;
			float rayDz = transformedRay.dir.z;

			float t = -rayPz / rayDz;
			// hit the opposite face, or not the closest one
			if (t < 0) return false;
			// x is the hit point in the unit plane's plane
			Vec3f x = ray.p + t * ray.dir;
			if (x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1) {
				return false;
			}
			if (t < t_max && t > ShadowBias)
			{
				return true;
			}
			return false;
		}
		else if (TriObj* ptr = dynamic_cast<TriObj*>(root->GetNodeObj())) {
			return  ptr->ShadowRecursive(transformedRay, t_max);
		}
	}
	return false;
}