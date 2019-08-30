
#include "scene.h"
#include "objects.h"
#include "cyVector.h"
#include "cyColor.h"
#include <math.h>

Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
MaterialList materials;
LightList lights;
#define PI 3.14159265
int LoadScene(char const *filename);

void recursive(int _i, int _j,Node* root, Ray ray, bool &_bHit) {
	if (root->GetNumChild() <= 0) return;
	for (int i = 0; i < root->GetNumChild(); i++)
	{
		if (root->GetChild(i)->GetNodeObj() != nullptr) {
			
			Ray transformedRay = root->GetChild(i)->ToNodeCoords(ray);
			recursive(_i, _j,root->GetChild(i), transformedRay, _bHit);

			HitInfo outHit;
			outHit.node = root->GetChild(i);

			if (root->GetChild(i)->GetNodeObj()->IntersectRay(transformedRay, outHit, 1))
			{
				renderImage.GetPixels()[_j*camera.imgWidth + _i].Set(255, 255, 255);
				renderImage.GetZBuffer()[_j*camera.imgWidth + _i] = outHit.z;

				_bHit = true;
			}
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
	Ray tRay;
	tRay.p = rayStart;

	renderImage.ResetNumRenderedPixels();
	for (size_t i = 0; i < camera.imgWidth; i++)
	{
		for (size_t j = 0; j < camera.imgHeight; j++)
		{
			Vec3f pixelPos = topLeft + (i + 1 / 2) * w / camera.imgWidth * camXAxis - (j + 1 / 2) * h / camera.imgHeight * camYAxis;
			tRay.dir = pixelPos - rayStart;
			// For this ray, if it hits or not
			bool bHit = false;
			recursive(i, j, &rootNode, tRay, bHit);
			// if it is not hit, write as black and the z buffer is big float
			if (!bHit) {
				renderImage.GetPixels()[j*camera.imgWidth + i].Set(0, 0, 0);
				renderImage.GetZBuffer()[j*camera.imgWidth + i] = BIGFLOAT;
			}
			renderImage.IncrementNumRenderPixel(1);
		}
	}
	renderImage.ComputeZBufferImage();
	renderImage.SaveImage("Resource/Result/proj2.png");
	renderImage.SaveZImage("Resource/Result/proj2_z.png");
}
void StopRender() {

}


int main() {
	const char* filename = "Resource/Data/proj2.xml";
	LoadScene(filename);

	printf("Render image width: %d\n", renderImage.GetWidth());
	printf("Render image height: %d", renderImage.GetHeight());

	ShowViewport();
	return 0;
}