
#include "scene.h"
#include "objects.h"
#include "lights.h"
#include "cyVector.h"
#include "cyColor.h"
#include <math.h>
#include "lights.h"
#include <omp.h>
#include <algorithm>
#include "cyPhotonMap.h"

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
Vec3f camZAxis;
Vec3f camYAxis;
Vec3f camXAxis;

#define PI 3.14159265

#define INTERNAL_REFLECTION_BOUNCE 15

int LoadScene(char const *filename);
void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);

void ShowViewport();
void SaveImages();

//-------------------
/** Build PhotonMap */
#define MAX_PhotonCount 10000
#define PhotonBounceCount 3
PhotonMap photonMap;
bool BuildPhotonMap(PhotonMap& o_photonMap, const LightList& i_lights, Node* i_root);
//-------------------
/** Calculate all light intensity*/
float allLightIntensity;
void CalculateLightsIntensity() {
	sort(lights.begin(), lights.end(), CompareGenLight);

	for (int i = 0; i < lights.size(); ++i)
	{
		allLightIntensity += ((Light*)lights[i])->GetIntensity();
	}
}
//---------------
// Use Components

#define USE_PathTracing
#define USE_GamaCorrection

#define GIBounceCount 3
//---------------
Vec3f RandomPositionInPixel(Vec3f i_center, float i_pixelLength) {
	Vec3f result = i_center;
	const Vec3f unit_dx = dd_x.GetNormalized();
	const Vec3f unit_dy = dd_y.GetNormalized();
	result += unit_dx * (((double)rand() / (RAND_MAX)) * 2 - 1)* i_pixelLength / 2;
	result += unit_dy * (((double)rand() / (RAND_MAX)) * 2 - 1)* i_pixelLength / 2;
	return result;
}
#ifdef USE_PathTracing
#define PT_SampleCount 16

Color PathTracing(int i_i, int i_j) {
	Color outColor = Color::Black();
	Vec3f pixelCenter = topLeft + (i_i + 1 / 2) * dd_x - (i_j + 1 / 2) * dd_y;

	/** Square Pixel Only*/
	const float pixelLen = dd_x.Length();
	Color colorSum = Color::Black();
	for (int i = 0; i < PT_SampleCount; ++i)
	{
		/** Jitter sampling in a pixel */
		Ray ray = Ray(camera.pos, RandomPositionInPixel(pixelCenter, pixelLen) - camera.pos);

		bool bHit = false;
		HitInfo tHitInfo = HitInfo();

		recursive(&rootNode, ray, tHitInfo, bHit, HIT_FRONT);
		if (bHit) {
			// Shade the hit object 
			colorSum += tHitInfo.node->GetMaterial()->Shade(ray, tHitInfo, lights, INTERNAL_REFLECTION_BOUNCE, GIBounceCount);

		}
		else {
			// Shade Background color
			Vec3f bguvw = Vec3f((float)i_i / camera.imgWidth, (float)i_j / camera.imgHeight, 0.0f);
			colorSum += background.Sample(bguvw);
		}
	}
	outColor = colorSum / PT_SampleCount;
	return outColor;
}

#endif // USE_PathTracing

bool CompareGenLight(Light* l1, Light* l2) {
	return (l1->GetIntensity() < l2->GetIntensity());
}

void BeginRender() {
	float aor = camera.imgWidth / (float)camera.imgHeight;
	float tan_h_pov = tan(camera.fov / 2 * PI / 180.0);
	float l = camera.focaldist;
	float h = 2 * l * tan_h_pov;
	float w = aor * h;

	camZAxis = -camera.dir;
	camYAxis = camera.up;
	camXAxis = camYAxis.Cross(camZAxis);

	topLeft = camera.pos - camZAxis * l + camYAxis * h / 2 - camXAxis * w / 2;

	dd_x = camXAxis * w / camera.imgWidth;
	dd_y = camYAxis * h / camera.imgHeight;
	renderImage.ResetNumRenderedPixels();

	CalculateLightsIntensity();

	BuildPhotonMap();

#pragma omp parallel for
	for (int i = 0; i < camera.imgWidth; ++i)
	{
		for (int j = 0; j < camera.imgHeight; ++j)
		{
			// Initialize the data
			{
				renderImage.GetSampleCount()[j*camera.imgWidth + i] = 0;
			}
			Color outColor = Color::Black();

			outColor = PathTracing(i, j);

#ifdef USE_GamaCorrection
			Color afterGamaCorrection = Color::Black();
			const float inverseGama = 1 / 2.2f;
			afterGamaCorrection.r = pow(outColor.r, inverseGama);
			afterGamaCorrection.g = pow(outColor.g, inverseGama);
			afterGamaCorrection.b = pow(outColor.b, inverseGama);
			outColor = afterGamaCorrection;
#endif // USE_GamaCorrection
			// Set out color
			renderImage.GetPixels()[j*camera.imgWidth + i] = Color24(outColor);
			//renderImage.GetZBuffer()[j*camera.imgWidth + i] = outHit.z;
			renderImage.IncrementNumRenderPixel(1);
		}
	}
	SaveImages();
}
void StopRender() {

}


bool BuildPhotonMap(PhotonMap& o_photonMap, const LightList& i_lights, Node* i_root)
{
	for (auto it = i_lights.begin(); it != i_lights.end(); ++ it)
	{
		// if it is a point light
		if (PointLight* ptr = reinterpret_cast<PointLight*>(*it)) {

			for (int i = 0; i < MAX_PhotonCount; ++i)
			{
				Ray ray = ptr->RandomPhoton();

			}

		}
	}
}


void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/) {
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
void SaveImages() {
	//renderImage.ComputeZBufferImage();
	renderImage.SaveImage("Resource/Result/proj13.png");
}
int main() {

	omp_set_num_threads(16);
	const char* filename = "Resource/Data/proj13.xml";
	LoadScene(filename);

	printf("Render image width: %d\n", renderImage.GetWidth());
	printf("Render image height: %d", renderImage.GetHeight());
	ShowViewport();

	return 0;
}