
#include "scene.h"
#include "objects.h"
#include "cyVector.h"
#include "cyColor.h"
#include <math.h>
#include "lights.h"
#include <omp.h>

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


#define REFLECTION_BOUNCE 5

int LoadScene(char const *filename);
void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
void ShowViewport();
void SaveImages();
Color TraceRaySingle(HitInfo &outHit, int i, int j);
Color TraceRayMultiple(int i, int j);

//---------------
// Use Components

//#define USE_MSAA
//#define ENABLE_DEPTH_OF_VIEW
#define USE_GamaCorrection

#define GIBounceCount 1
//---------------


#ifdef USE_MSAA
#define MSAA_AdaptiveThreshold_Sqr 0.03f * 0.03f
#define  MSAA_RayCountPerSlot 4
#define  MSAA_AddadptiveStopLevel 3
Vec3f RandomPositionInPixel(Vec3f i_center, float i_pixelLength) {
	Vec3f result = i_center;
	result += dd_x.GetNormalized() *(((double)rand() / (RAND_MAX)) * 2 - 1)* i_pixelLength / 2;
	result += dd_y.GetNormalized() *(((double)rand() / (RAND_MAX)) * 2 - 1)* i_pixelLength / 2;
	return result;
}
Color JitteredAddaptiveSampling(Vec3f pixelCenter, int level, int i_i, int i_j) {
	level++;
	renderImage.GetSampleCount()[i_j*camera.imgWidth + i_i] += MSAA_RayCountPerSlot;
	Color overallColor = Color::Black();
	// square of variance
	float varianceSqr[3] = { 1 ,1, 1 };
	// average of r,g,b channel
	float average[3];

	Color subPixelColor[] = { Color::Black() ,Color::Black() ,Color::Black() ,Color::Black() };
	Color subPixelColor_Sum = Color::Black();
	// Trace 4 rays for 4 sub-pixels
	for (int i = 0; i < 4; ++i)
	{
		Ray tRay;
		tRay.p = camera.pos;

		switch (i)
		{
		case 0:
			tRay.dir = RandomPositionInPixel(pixelCenter - dd_x / (float)(level * MSAA_RayCountPerSlot) - dd_y / (float)(level * MSAA_RayCountPerSlot), dd_x.Length() * 2 / (float)(level * MSAA_RayCountPerSlot)) - tRay.p;
			break;
		case 1:
			tRay.dir = RandomPositionInPixel(pixelCenter + dd_x / (float)(level * MSAA_RayCountPerSlot) - dd_y / (float)(level * MSAA_RayCountPerSlot), dd_x.Length() * 2 / (float)(level * MSAA_RayCountPerSlot)) - tRay.p;
			break;
		case 2:
			tRay.dir = RandomPositionInPixel(pixelCenter + dd_x / (float)(level * MSAA_RayCountPerSlot) + dd_y / (float)(level * MSAA_RayCountPerSlot), dd_x.Length() * 2 / (float)(level * MSAA_RayCountPerSlot)) - tRay.p;
			break;
		case 3:
			tRay.dir = RandomPositionInPixel(pixelCenter - dd_x / (float)(level * MSAA_RayCountPerSlot) + dd_y / (float)(level * MSAA_RayCountPerSlot), dd_x.Length() * 2 / (float)(level * MSAA_RayCountPerSlot)) - tRay.p;
			break;
		default:
			break;
		}

		bool bHit = false;
		HitInfo tHitInfo = HitInfo();

		recursive(&rootNode, tRay, tHitInfo, bHit, HIT_FRONT);
		if (bHit) {
			// Shade the hit object 
			subPixelColor[i] = tHitInfo.node->GetMaterial()->Shade(tRay, tHitInfo, lights, REFLECTION_BOUNCE, GIBounceCount);

		}
		else {
			// Shade Background color
			Vec3f bguvw = Vec3f((float)i_i / camera.imgWidth, (float)i_j / camera.imgHeight, 0.0f);
			subPixelColor[i] = background.Sample(bguvw);
		}
		subPixelColor_Sum += subPixelColor[i];
	}
	// Get the average color of 4 sub-pixels
	average[0] = subPixelColor_Sum.r / (MSAA_RayCountPerSlot);
	average[1] = subPixelColor_Sum.g / (MSAA_RayCountPerSlot);
	average[2] = subPixelColor_Sum.b / (MSAA_RayCountPerSlot);

	float variance_Sum[3] = { 0,0,0 };
	for (int i = 0; i < 4; i++)
	{
		variance_Sum[0] += (subPixelColor[i].r - average[0]) * (subPixelColor[i].r - average[0]);
		variance_Sum[1] += (subPixelColor[i].g - average[1]) * (subPixelColor[i].g - average[1]);
		variance_Sum[2] += (subPixelColor[i].b - average[2]) * (subPixelColor[i].b - average[2]);
	}
	varianceSqr[0] = variance_Sum[0] / MSAA_RayCountPerSlot;
	varianceSqr[1] = variance_Sum[1] / MSAA_RayCountPerSlot;
	varianceSqr[2] = variance_Sum[2] / MSAA_RayCountPerSlot;

	float widthPerSubPixel = level * MSAA_RayCountPerSlot;
	// top left
	if (
		(varianceSqr[0] > MSAA_AdaptiveThreshold_Sqr ||
			varianceSqr[1] > MSAA_AdaptiveThreshold_Sqr ||
			varianceSqr[2] > MSAA_AdaptiveThreshold_Sqr)
		&& level < MSAA_AddadptiveStopLevel
		) {
		// return the mean color of all sub pixel
		return (
			JitteredAddaptiveSampling(pixelCenter - dd_x / widthPerSubPixel - dd_y / widthPerSubPixel, level, i_i, i_j) +
			JitteredAddaptiveSampling(pixelCenter + dd_x / widthPerSubPixel - dd_y / widthPerSubPixel, level, i_i, i_j) +
			JitteredAddaptiveSampling(pixelCenter - dd_x / widthPerSubPixel + dd_y / widthPerSubPixel, level, i_i, i_j) +
			JitteredAddaptiveSampling(pixelCenter + dd_x / widthPerSubPixel + dd_y / widthPerSubPixel, level, i_i, i_j)
			) / MSAA_RayCountPerSlot;
	}

	return subPixelColor_Sum / MSAA_RayCountPerSlot;
}
#endif // USE_MSAA

#ifdef ENABLE_DEPTH_OF_VIEW
#define DoV_SampleCount 128
Vec3f GetSampleInAperture(const Camera& cam) {
	// Center of aperture
	Vec3f O = cam.pos;
	// Radius of aperture
	float R = cam.dof;
	// Non-uniform distribution
	float r = ((double)rand() / (RAND_MAX));
	// Uniform distribution
	r = sqrt(r) * R;
	float theta = ((double)rand() / (RAND_MAX)) * 2 * PI;
	// Random point in a circle
	float x = r * cos(theta);
	float y = r * sin(theta);

	return O + camXAxis * x + camYAxis * y;
}
Color DOVSampling(Vec3f pixelCenter, int i_i, int i_j) {
	Color colorSum = Color::Black();
	for (int i = 0; i < DoV_SampleCount; ++i)
	{
		Ray ray = Ray();
		ray.p = GetSampleInAperture(camera);
		ray.dir = pixelCenter - ray.p;

		bool bHit = false;
		HitInfo tHitInfo = HitInfo();

		recursive(&rootNode, ray, tHitInfo, bHit, HIT_FRONT);
		if (bHit) {
			// Shade the hit object 
			colorSum += tHitInfo.node->GetMaterial()->Shade(ray, tHitInfo, lights, REFLECTION_BOUNCE);

		}
		else {
			// Shade Background color
			Vec3f bguvw = Vec3f((float)i_i / camera.imgWidth, (float)i_j / camera.imgHeight, 0.0f);
			colorSum += background.Sample(bguvw);
		}
	}
	return colorSum / DoV_SampleCount;
}
#endif // ENABLE_DEPTH_OF_VIEW



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
#ifdef USE_MSAA | ENABLE_DEPTH_OF_VIEW
			outColor = TraceRayMultiple(i, j);
#else
			HitInfo outHit;
			outColor = TraceRaySingle(outHit, i, j);
#endif // USE_MSAA
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


// Trace only single ray
Color TraceRaySingle(HitInfo &outHit, int i, int j)
{
	Vec3f pixelCenter = topLeft + (i + 1 / 2) * dd_x - (j + 1 / 2) * dd_y;
	Ray ray;
	ray.p = camera.pos;
	ray.dir = pixelCenter - ray.p;
	// For this ray, if it hits or not
	bool bHit = false;
	recursive(&rootNode, ray, outHit, bHit, HIT_FRONT);
	if (bHit) {
		// Shade the hit object 
		return outHit.node->GetMaterial()->Shade(ray, outHit, lights, REFLECTION_BOUNCE, GIBounceCount);
	}
	else {
		// Shade Background color
		Vec3f bguvw = Vec3f((float)i / camera.imgWidth, (float)j / camera.imgHeight, 0.0f);
		return background.Sample(bguvw);
	}
}
// Multi-Sampling
Color TraceRayMultiple(int i, int j)
{
	Vec3f pixelCenter = topLeft + (i + 1 / 2) * dd_x - (j + 1 / 2) * dd_y;
	Color overallColor = Color::Black();
#ifdef ENABLE_DEPTH_OF_VIEW
	overallColor = DOVSampling(pixelCenter, i, j);
#else
#ifdef USE_MSAA
	overallColor = JitteredAddaptiveSampling(pixelCenter, 0, i, j);
#endif // USE_MSAA
#endif // ENABLE_DEPTH_OF_VIEW

	return overallColor;
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
	renderImage.SaveImage("Resource/Result/proj11_2.png");
	renderImage.ComputeSampleCountImage();
	renderImage.SaveSampleCountImage("Resource/Result/proj11_2_sample.png");
}
int main() {

	omp_set_num_threads(16);
	const char* filename = "Resource/Data/proj11_2.xml";
	LoadScene(filename);

	printf("Render image width: %d\n", renderImage.GetWidth());
	printf("Render image height: %d", renderImage.GetHeight());
	ShowViewport();

	return 0;
}