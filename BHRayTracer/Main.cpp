
#include "scene.h"
#include "objects.h"
Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
int LoadScene(char const *filename);
void ShowViewport();
void BeginRender() {

}
void StopRender() {

}

int main() {
	const char* filename = "data/proj1.xml";
	LoadScene(filename);

	ShowViewport();
	return 0;
}