
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
// Lights & Materials
MaterialList materials;
LightList lights;
// Triangular meshes
ObjFileList objList;
// Textures
TexturedColor background;
TexturedColor environment;
TextureList textureList;

Vec3f topLeft;
Vec3f dd_x;
Vec3f dd_y;
#define PI 3.14159265
#define Bias 0.001f
#define REFLECTION_BOUNCE 3

//#define USE_MSAA
#ifdef USE_MSAA
#define MSAA_AdaptiveThreshold 0.2f
#define  MSAA_InitialRayCount 4
#endif // USE_MSAA

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

// Trace only single ray
Color TraceRaySingle(Ray &ray, HitInfo &outHit, int i, int j)
{
	Vec3f pixelPos = topLeft + (i + 1 / 2) * dd_x - (j + 1 / 2) * dd_y;
	ray.dir = pixelPos - ray.p;
	// For this ray, if it hits or not
	bool bHit = false;
	recursive(&rootNode, ray, outHit, bHit, 0);
	if (bHit) {
		// Shade the hit object 
		return outHit.node->GetMaterial()->Shade(ray, outHit, lights, REFLECTION_BOUNCE);
	}
	else {
		// Shade Background color
		Vec3f bguvw = Vec3f((float)i / camera.imgWidth, (float)j / camera.imgHeight, 0.0f);
		return background.Sample(bguvw);
	}
}
// Multi-Sampling
Color TraceRayMultiple(Ray &ray, HitInfo &outHit, int i, int j) { return Color::Black(); }

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

	topLeft = rayStart - camZAxis * l + camYAxis * h / 2 - camXAxis * w / 2;

	dd_x = camXAxis * w / camera.imgWidth;
	dd_y = camYAxis * h / camera.imgHeight;

	renderImage.ResetNumRenderedPixels();
	//#pragma omp parallel for
	for (int i = 0; i < camera.imgWidth; ++i)
	{
		for (int j = 0; j < camera.imgHeight; ++j)
		{
			Color outColor = Color::Black();
			Ray tRay;
			tRay.p = rayStart;
			HitInfo outHit;
			
#ifdef USE_MSAA
			outColor = TraceRayMultiple(tRay, outHit, i, j);
#else
			outColor = TraceRaySingle(tRay, outHit, i, j);
#endif // USE_MSAA

			// Set out color
			renderImage.GetPixels()[j*camera.imgWidth + i] = Color24(outColor);
			renderImage.GetZBuffer()[j*camera.imgWidth + i] = outHit.z;
			renderImage.IncrementNumRenderPixel(1);
			//printf("Percent: %f\n", renderImage.GetNumRenderedPixels() / (float)(renderImage.GetWidth() * renderImage.GetHeight()));
		}
	}

	renderImage.ComputeZBufferImage();
	renderImage.SaveImage("Resource/Result/prj7.png");
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
			return ptr->ShadowRecursive(transformedRay, t_max);
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
	const char* filename = "Resource/Data/proj7.xml";
	LoadScene(filename);

	printf("Render image width: %d\n", renderImage.GetWidth());
	printf("Render image height: %d", renderImage.GetHeight());
	ShowViewport();

	return 0;
}