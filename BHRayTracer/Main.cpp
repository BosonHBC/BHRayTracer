
#include "scene.h"
#include "objects.h"
#include "cyVector.h"
#include "cyColor.h"
#include <math.h>
#include "lights.h"
//#include <omp.h>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
Plane thePlane;
MaterialList materials;
LightList lights;
ObjFileList objList;

#define PI 3.14159265
#define  Bias 0.00001f
#define REFLECTION_BOUNCE 3


int LoadScene(char const *filename);

void recursive(Node* root, Ray ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/) {
	if (root->GetNumChild() <= 0) return;
	for (int i = 0; i < root->GetNumChild(); i++)
	{
		Ray transformedRay = root->GetChild(i)->ToNodeCoords(ray);
		if (root->GetChild(i)->GetNodeObj() != nullptr) {

			// transform ray to child coordinate
			if (root->GetChild(i)->GetNodeObj()->IntersectRay(transformedRay, outHit, hitSide))
			{
				outHit.node = root->GetChild(i);
				_bHit = true;
				root->GetChild(i)->FromNodeCoords(outHit);
			}

		}
		recursive(root->GetChild(i), transformedRay, outHit, _bHit, hitSide);
	}
	for (int i = 0; i < root->GetNumChild(); i++) {
		if (root->GetChild(i) == outHit.node) {
			root->FromNodeCoords(outHit);
			break;
		}
	}
}

void ShowViewport();
void BeginRender() {
	Vec3f rayStart = camera.pos;
	float aor = camera.imgWidth / (float)camera.imgHeight;
	float tan_h_pov = tan(camera.fov / 2 * PI / 180.0);
	float l = 1;
	float h = 2 * l * tan_h_pov;
	float w = aor * h;

	Vec3f camZAxis = -camera.dir;
	Vec3f camYAxis = camera.up;
	Vec3f camXAxis = camYAxis.Cross(camZAxis);
	Vec3f topLeft = rayStart - camZAxis * l + camYAxis * h / 2 - camXAxis * w / 2;


	renderImage.ResetNumRenderedPixels();
	//#pragma omp parallel for
	for (int i = 0; i < camera.imgWidth; ++i)
	{
		for (int j = 0; j < camera.imgHeight; ++j)
		{
			Vec3f pixelPos = topLeft + (i + 1 / 2) * w / camera.imgWidth * camXAxis - (j + 1 / 2) * h / camera.imgHeight * camYAxis;
			Ray tRay;
			tRay.p = rayStart;
			tRay.dir = pixelPos - rayStart;
			// For this ray, if it hits or not
			bool bHit = false;
			HitInfo outHit;
			recursive(&rootNode, tRay, outHit, bHit, 0);
			if (bHit) {
				Ray worldRay;
				worldRay.dir = pixelPos - rayStart;
				worldRay.p = rayStart;
				Color outColor = outHit.node->GetMaterial()->Shade(worldRay, outHit, lights, REFLECTION_BOUNCE);

				renderImage.GetPixels()[j*camera.imgWidth + i] = Color24(outColor);
				renderImage.GetZBuffer()[j*camera.imgWidth + i] = outHit.z;
			}
			else {
				renderImage.GetPixels()[j*camera.imgWidth + i].Set(0, 0, 0);
				renderImage.GetZBuffer()[j*camera.imgWidth + i] = BIGFLOAT;
			}
			renderImage.IncrementNumRenderPixel(1);
			//printf("Percent: %f\n", renderImage.GetNumRenderedPixels() / (float)(renderImage.GetWidth() * renderImage.GetHeight()));
		}
	}

	renderImage.ComputeZBufferImage();
	//renderImage.SaveZImage("Resource/Result/prj5.png");
	renderImage.SaveImage("Resource/Result/prj6.png");
}
void StopRender() {

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
				if (t < t_max && t > Bias)
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
			if (t < t_max && t > Bias)
			{
				return true;
			}
			return false;
		}
		else if (TriObj* ptr = dynamic_cast<TriObj*>(root->GetNodeObj())) {
			HitInfo hit;
			return ptr->ShadowRecursive(transformedRay, hit, t_max);
		}
	}
	return false;
}
float GenLight::Shadow(Ray ray, float t_max /*= BIGFLOAT*/)
{
	return ShadowRayRecursive(&rootNode, ray, t_max) ? 0.f : 1.f;
}



int main() {

	//	omp_set_num_threads(16);
	const char* filename = "Resource/Data/proj5.xml";
	LoadScene(filename);

	printf("Render image width: %d\n", renderImage.GetWidth());
	printf("Render image height: %d", renderImage.GetHeight());
	ShowViewport();

	return 0;
}