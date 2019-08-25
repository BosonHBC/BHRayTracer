
#include "scene.h"
#include "objects.h"
Node rootNode;
Camera camera;
RenderImage renderImage;
Sphere theSphere;
int LoadScene(char const *filename);

void BeginRender() {

}
void StopRender() {

}

int main() {
	const char* filename = "testscene";
	LoadScene(filename);

	return 0;
}