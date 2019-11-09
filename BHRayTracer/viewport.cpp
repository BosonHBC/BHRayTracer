
//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp 
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    11.0
/// \date       August 21, 2019
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "scene.h"
#include "objects.h"
#include "lights.h"
#include "materials.h"
#include "texture.h"
#include <stdlib.h>
#include <time.h>

#ifdef USE_GLUT
# ifdef __APPLE__
#  include <GLUT/glut.h>
# else
#  include <GL/glut.h>
# endif
#else
# include <GL/freeglut.h>
#endif

//-------------------------------------------------------------------------------

void BeginRender(); // Called to start rendering (renderer must run in a separate thread)
void StopRender();  // Called to end rendering (if it is not already finished)

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern LightList lights;
extern TexturedColor background;

//-------------------------------------------------------------------------------

enum Mode {
	MODE_READY,         // Ready to render
	MODE_RENDERING,     // Rendering the image
	MODE_RENDER_DONE    // Rendering is finished
};

enum ViewMode
{
	VIEWMODE_OPENGL,
	VIEWMODE_IMAGE,
	VIEWMODE_Z,
	VIEWMODE_SAMPLECOUNT,
	VIEWMODE_IRRADCOMP,
};

enum MouseMode {
	MOUSEMODE_NONE,
	MOUSEMODE_DEBUG,
	MOUSEMODE_ROTATE,
};

static Mode     mode = MODE_READY;       // Rendering mode
static ViewMode viewMode = VIEWMODE_OPENGL;  // Display mode
static MouseMode mouseMode = MOUSEMODE_NONE;   // Mouse mode
static int      startTime;                      // Start time of rendering
static int mouseX = 0, mouseY = 0;
static float viewAngle1 = 0, viewAngle2 = 0;
static GLuint viewTexture;

static int      dofDrawCount = 0;
static Color    *dofImage = nullptr;
static Color24 *dofBuffer = nullptr;

#define MAX_DOF_DRAW    32

//-------------------------------------------------------------------------------

void GlutDisplay();
void GlutReshape(int w, int h);
void GlutIdle();
void GlutKeyboard(unsigned char key, int x, int y);
void GlutMouse(int button, int state, int x, int y);
void GlutMotion(int x, int y);

//-------------------------------------------------------------------------------

void ShowViewport()
{
	int argc = 1;
	char argstr[] = "raytrace";
	char *argv = argstr;
	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	if (glutGet(GLUT_SCREEN_WIDTH) > 0 && glutGet(GLUT_SCREEN_HEIGHT) > 0) {
		glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - camera.imgWidth) / 2, (glutGet(GLUT_SCREEN_HEIGHT) - camera.imgHeight) / 2);
	}
	else glutInitWindowPosition(50, 50);
	glutInitWindowSize(camera.imgWidth, camera.imgHeight);

	glutCreateWindow("Ray Tracer - CS 6620");
	glutDisplayFunc(GlutDisplay);
	glutReshapeFunc(GlutReshape);
	glutIdleFunc(GlutIdle);
	glutKeyboardFunc(GlutKeyboard);
	glutMouseFunc(GlutMouse);
	glutMotionFunc(GlutMotion);

	Color bg = background.GetColor();
	glClearColor(bg.r, bg.g, bg.b, 0);

	glPointSize(3.0);
	glEnable(GL_CULL_FACE);

	float zero[] = { 0,0,0,0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, zero);

	glEnable(GL_NORMALIZE);

	glLineWidth(2);

	if (camera.dof > 0) {
		dofBuffer = new Color24[camera.imgWidth * camera.imgHeight];
		dofImage = new Color[camera.imgWidth * camera.imgHeight];
		memset(dofImage, 0, camera.imgWidth * camera.imgHeight * sizeof(Color));
	}

	glGenTextures(1, &viewTexture);
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	glutMainLoop();
}

//-------------------------------------------------------------------------------

void GlutReshape(int w, int h)
{
	if (w != camera.imgWidth || h != camera.imgHeight) {
		glutReshapeWindow(camera.imgWidth, camera.imgHeight);
	}
	else {
		glViewport(0, 0, w, h);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		float r = (float)w / float(h);
		gluPerspective(camera.fov, r, 0.02, 1000.0);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
}

//-------------------------------------------------------------------------------

void DrawNode(Node *node)
{
	glPushMatrix();

	const Material *mtl = node->GetMaterial();
	if (mtl) mtl->SetViewportMaterial();

	Matrix3f tm = node->GetTransform();
	Vec3f p = node->GetPosition();
	float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
	glMultMatrixf(m);

	Object *obj = node->GetNodeObj();
	if (obj) obj->ViewportDisplay(mtl);

	for (int i = 0; i < node->GetNumChild(); i++) {
		DrawNode(node->GetChild(i));
	}

	glPopMatrix();
}

//-------------------------------------------------------------------------------

void DrawScene()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	const TextureMap *bgMap = background.GetTexture();
	if (bgMap) {
		glDepthMask(GL_FALSE);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		Color c = background.GetColor();
		glColor3f(c.r, c.g, c.b);
		if (bgMap->SetViewportTexture()) {
			glEnable(GL_TEXTURE_2D);
			glMatrixMode(GL_TEXTURE);
			Matrix3f tm = bgMap->GetInverseTransform();
			Vec3f p = tm * bgMap->GetPosition();
			float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
			glLoadMatrixf(m);
			glMatrixMode(GL_MODELVIEW);
		}
		else {
			glDisable(GL_TEXTURE_2D);
		}
		glBegin(GL_QUADS);
		glTexCoord2f(0, 1);
		glVertex2f(-1, -1);
		glTexCoord2f(1, 1);
		glVertex2f(1, -1);
		glTexCoord2f(1, 0);
		glVertex2f(1, 1);
		glTexCoord2f(0, 0);
		glVertex2f(-1, 1);
		glEnd();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glDepthMask(GL_TRUE);

		glDisable(GL_TEXTURE_2D);
	}

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glPushMatrix();
	Vec3f p = camera.pos;
	Vec3f t = camera.pos + camera.dir*camera.focaldist;
	Vec3f u = camera.up;
	if (camera.dof > 0) {
		Vec3f v = camera.dir ^ camera.up;
		float r = sqrtf(float(rand()) / RAND_MAX)*camera.dof;
		float a = float(M_PI) * 2.0f * float(rand()) / RAND_MAX;
		p += r * cosf(a)*v + r * sinf(a)*u;
	}
	gluLookAt(p.x, p.y, p.z, t.x, t.y, t.z, u.x, u.y, u.z);

	glRotatef(viewAngle1, 1, 0, 0);
	glRotatef(viewAngle2, 0, 0, 1);

	if (lights.size() > 0) {
		for (unsigned int i = 0; i < lights.size(); i++) {
			lights[i]->SetViewportLight(i);
		}
	}
	else {
		float white[] = { 1,1,1,1 };
		float black[] = { 0,0,0,0 };
		Vec4f p(camera.pos, 1);
		glEnable(GL_LIGHT0);
		glLightfv(GL_LIGHT0, GL_AMBIENT, black);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
		glLightfv(GL_LIGHT0, GL_SPECULAR, white);
		glLightfv(GL_LIGHT0, GL_POSITION, &p.x);
	}

	DrawNode(&rootNode);

	glPopMatrix();

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------------

void DrawImage(void *data, GLenum type, GLenum format)
{
	glBindTexture(GL_TEXTURE_2D, viewTexture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, renderImage.GetWidth(), renderImage.GetHeight(), 0, format, type, data);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);

	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glColor3f(1, 1, 1);
	glBegin(GL_QUADS);
	glTexCoord2f(0, 1);
	glVertex2f(-1, -1);
	glTexCoord2f(1, 1);
	glVertex2f(1, -1);
	glTexCoord2f(1, 0);
	glVertex2f(1, 1);
	glTexCoord2f(0, 0);
	glVertex2f(-1, 1);
	glEnd();

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	glDisable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------------

void DrawProgressBar(float done)
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glBegin(GL_LINES);
	glColor3f(1, 1, 1);
	glVertex2f(-1, -1);
	glVertex2f(done * 2 - 1, -1);
	glColor3f(0, 0, 0);
	glVertex2f(done * 2 - 1, -1);
	glVertex2f(1, -1);
	glEnd();

	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

//-------------------------------------------------------------------------------

void DrawRenderProgressBar()
{
	int rp = renderImage.GetNumRenderedPixels();
	int np = renderImage.GetWidth() * renderImage.GetHeight();
	if (rp >= np) return;
	float done = (float)rp / (float)np;
	DrawProgressBar(done);
}

//-------------------------------------------------------------------------------

void GlutDisplay()
{
	switch (viewMode) {
	case VIEWMODE_OPENGL:
		if (dofImage) {
			if (dofDrawCount < MAX_DOF_DRAW) {
				DrawScene();
				glReadPixels(0, 0, camera.imgWidth, camera.imgHeight, GL_RGB, GL_UNSIGNED_BYTE, dofBuffer);
				for (int i = 0, y = 0; y < camera.imgHeight; y++) {
					int j = (camera.imgHeight - y - 1)*camera.imgWidth;
					for (int x = 0; x < camera.imgWidth; x++, i++, j++) {
						dofImage[i] = (dofImage[i] * float(dofDrawCount) + dofBuffer[j].ToColor()) / float(dofDrawCount + 1);
					}
				}
				dofDrawCount++;
			}
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
			DrawImage(dofImage, GL_FLOAT, GL_RGB);
			if (dofDrawCount < MAX_DOF_DRAW) {
				DrawProgressBar(float(dofDrawCount) / MAX_DOF_DRAW);
				glutPostRedisplay();
			}
		}
		else {
			DrawScene();
		}
		break;
	case VIEWMODE_IMAGE:
		DrawImage(renderImage.GetPixels(), GL_UNSIGNED_BYTE, GL_RGB);
		DrawRenderProgressBar();
		break;
	case VIEWMODE_Z:
		if (!renderImage.GetZBufferImage()) renderImage.ComputeZBufferImage();
		DrawImage(renderImage.GetZBufferImage(), GL_UNSIGNED_BYTE, GL_LUMINANCE);
		break;
	case VIEWMODE_SAMPLECOUNT:
		if (!renderImage.GetSampleCountImage()) renderImage.ComputeSampleCountImage();
		DrawImage(renderImage.GetSampleCountImage(), GL_UNSIGNED_BYTE, GL_LUMINANCE);
		break;
	case VIEWMODE_IRRADCOMP:
		if (renderImage.GetIrradianceComputationImage()) {
			DrawImage(renderImage.GetIrradianceComputationImage(), GL_UNSIGNED_BYTE, GL_LUMINANCE);
		}
		break;
	}

	glutSwapBuffers();
}

//-------------------------------------------------------------------------------

void GlutIdle()
{
	static int lastRenderedPixels = 0;
	if (mode == MODE_RENDERING) {
		int nrp = renderImage.GetNumRenderedPixels();
		if (lastRenderedPixels != nrp) {
			lastRenderedPixels = nrp;
			if (renderImage.IsRenderDone()) {
				mode = MODE_RENDER_DONE;
				int endTime = (int)time(nullptr);
				int t = endTime - startTime;
				int h = t / 3600;
				int m = (t % 3600) / 60;
				int s = t % 60;
				printf("\nRender time is %d:%02d:%02d.\n", h, m, s);
			}
			glutPostRedisplay();
		}
		if (viewMode == VIEWMODE_IRRADCOMP) glutPostRedisplay();
	}
}

//-------------------------------------------------------------------------------

void GlutKeyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 27:    // ESC
		exit(0);
		break;
	case ' ':
		switch (mode) {
		case MODE_READY:
			mode = MODE_RENDERING;
			viewMode = VIEWMODE_IMAGE;
			if (dofImage) {
				Color24 *p = renderImage.GetPixels();
				for (int i = 0; i < camera.imgWidth*camera.imgHeight; i++) p[i] = Color24(dofImage[i]);
			}
			else {
				DrawScene();
				glReadPixels(0, 0, renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels());
				{
					Color24 *c = renderImage.GetPixels();
					for (int y0 = 0, y1 = renderImage.GetHeight() - 1; y0 < y1; y0++, y1--) {
						int i0 = y0 * renderImage.GetWidth();
						int i1 = y1 * renderImage.GetWidth();
						for (int x = 0; x < renderImage.GetWidth(); x++, i0++, i1++) {
							Color24 t = c[i0]; c[i0] = c[i1]; c[i1] = t;
						}
					}
				}
			}
			startTime = (int)time(nullptr);
			BeginRender();
			break;
		case MODE_RENDERING:
			mode = MODE_READY;
			StopRender();
			glutPostRedisplay();
			break;
		case MODE_RENDER_DONE:
			mode = MODE_READY;
			viewMode = VIEWMODE_OPENGL;
			glutPostRedisplay();
			break;
		}
		break;
	case '1':
		viewAngle1 = viewAngle2 = 0;
		viewMode = VIEWMODE_OPENGL;
		glutPostRedisplay();
		break;
	case '2':
		viewMode = VIEWMODE_IMAGE;
		glutPostRedisplay();
		break;
	case '3':
		viewMode = VIEWMODE_Z;
		glutPostRedisplay();
		break;
	case '4':
		viewMode = VIEWMODE_SAMPLECOUNT;
		glutPostRedisplay();
		break;
	case '5':
		viewMode = VIEWMODE_IRRADCOMP;
		glutPostRedisplay();
		break;
	}
}

//-------------------------------------------------------------------------------

void PrintPixelData(int x, int y)
{
	if (x < renderImage.GetWidth() && y < renderImage.GetHeight()) {
		Color24 *colors = renderImage.GetPixels();
		float *zbuffer = renderImage.GetZBuffer();
		int i = (renderImage.GetHeight() - y - 1) *renderImage.GetWidth() + x;
		printf("Pixel [ %d, %d ] Color24: %d, %d, %d   Z: %f\n", x, y, colors[i].r, colors[i].g, colors[i].b, zbuffer[i]);
	}
	else {
		printf("-- Invalid pixel (%d,%d) --\n", x, y);
	}
}

//-------------------------------------------------------------------------------

void GlutMouse(int button, int state, int x, int y)
{
	if (state == GLUT_UP) {
		mouseMode = MOUSEMODE_NONE;
	}
	else {
		switch (button) {
		case GLUT_LEFT_BUTTON:
			mouseMode = MOUSEMODE_DEBUG;
			PrintPixelData(x, y);
			break;
		case GLUT_RIGHT_BUTTON:
			mouseMode = MOUSEMODE_ROTATE;
			mouseX = x;
			mouseY = y;
			break;
		}
	}
}

//-------------------------------------------------------------------------------

void GlutMotion(int x, int y)
{
	switch (mouseMode) {
	case MOUSEMODE_DEBUG:
		PrintPixelData(x, y);
		break;
	case GLUT_RIGHT_BUTTON:
		viewAngle1 -= 0.2f * (mouseY - y);
		viewAngle2 -= 0.2f * (mouseX - x);
		mouseX = x;
		mouseY = y;
		glutPostRedisplay();
		break;
	}
}

//-------------------------------------------------------------------------------
// Viewport Methods for various classes
//-------------------------------------------------------------------------------
void Sphere::ViewportDisplay(const Material *mtl) const
{
	static GLUquadric *q = nullptr;
	if (q == nullptr) {
		q = gluNewQuadric();
		gluQuadricTexture(q, true);
	}
	gluSphere(q, 1, 50, 50);
}
void Plane::ViewportDisplay(const Material *mtl) const
{
	const int resolution = 32;
	float xyInc = 2.0f / resolution;
	float uvInc = 1.0f / resolution;
	glPushMatrix();
	glNormal3f(0, 0, 1);
	glBegin(GL_QUADS);
	float y1 = -1, y2 = xyInc - 1, v1 = 0, v2 = uvInc;
	for (int y = 0; y < resolution; y++) {
		float x1 = -1, x2 = xyInc - 1, u1 = 0, u2 = uvInc;
		for (int x = 0; x < resolution; x++) {
			glTexCoord2f(u1, v1);
			glVertex3f(x1, y1, 0);
			glTexCoord2f(u2, v1);
			glVertex3f(x2, y1, 0);
			glTexCoord2f(u2, v2);
			glVertex3f(x2, y2, 0);
			glTexCoord2f(u1, v2);
			glVertex3f(x1, y2, 0);
			x1 = x2; x2 += xyInc; u1 = u2; u2 += uvInc;
		}
		y1 = y2; y2 += xyInc; v1 = v2; v2 += uvInc;
	}
	glEnd();
	glPopMatrix();
}
void TriObj::ViewportDisplay(const Material *mtl) const
{
	unsigned int nextMtlID = 0;
	unsigned int nextMtlSwith = NF();
	if (mtl && NM() > 0) {
		mtl->SetViewportMaterial(0);
		nextMtlSwith = GetMaterialFaceCount(0);
		nextMtlID = 1;
	}
	glBegin(GL_TRIANGLES);
	for (unsigned int i = 0; i < NF(); i++) {
		while (i >= nextMtlSwith) {
			if (nextMtlID >= NM()) nextMtlSwith = NF();
			else {
				glEnd();
				nextMtlSwith += GetMaterialFaceCount(nextMtlID);
				mtl->SetViewportMaterial(nextMtlID);
				nextMtlID++;
				glBegin(GL_TRIANGLES);
			}
		}
		for (int j = 0; j < 3; j++) {
			if (HasTextureVertices()) glTexCoord3fv(&VT(FT(i).v[j]).x);
			if (HasNormals()) glNormal3fv(&VN(FN(i).v[j]).x);
			glVertex3fv(&V(F(i).v[j]).x);
		}
	}
	glEnd();
}
void MtlBlinn::SetViewportMaterial(int subMtlID) const
{
	ColorA c;
	c = ColorA(diffuse.GetColor());
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r);
	c = ColorA(specular.GetColor());
	glMaterialfv(GL_FRONT, GL_SPECULAR, &c.r);
	glMaterialf(GL_FRONT, GL_SHININESS, glossiness*1.5f);
	c = ColorA(emission.GetColor());
	glMaterialfv(GL_FRONT, GL_EMISSION, &c.r);
	const TextureMap *dm = diffuse.GetTexture();
	if (dm && dm->SetViewportTexture()) {
		glEnable(GL_TEXTURE_2D);
		glMatrixMode(GL_TEXTURE);
		Matrix3f tm = dm->GetInverseTransform();
		Vec3f p = tm * dm->GetPosition();
		float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, -p.x,-p.y,-p.z,1 };
		glLoadMatrixf(m);
		glMatrixMode(GL_MODELVIEW);
	}
	else {
		glDisable(GL_TEXTURE_2D);
	}
}
void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Vec4f pos) const
{
	glEnable(GL_LIGHT0 + lightID);
	glLightfv(GL_LIGHT0 + lightID, GL_AMBIENT, &ambient.r);
	glLightfv(GL_LIGHT0 + lightID, GL_DIFFUSE, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r);
	glLightfv(GL_LIGHT0 + lightID, GL_POSITION, &pos.x);
}
bool TextureFile::SetViewportTexture() const
{
	if (viewportTextureID == 0) {
		glGenTextures(1, &viewportTextureID);
		glBindTexture(GL_TEXTURE_2D, viewportTextureID);
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, &data[0].r);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	glBindTexture(GL_TEXTURE_2D, viewportTextureID);
	return true;
}
bool TextureChecker::SetViewportTexture() const
{
	if (viewportTextureID == 0) {
		const int texSize = 256;
		glGenTextures(1, &viewportTextureID);
		glBindTexture(GL_TEXTURE_2D, viewportTextureID);
		Color24 c[2] = { Color24(color1), Color24(color2) };
		Color24 *tex = new Color24[texSize*texSize];
		for (int i = 0; i < texSize*texSize; i++) {
			int ix = (i%texSize) < 128 ? 0 : 1;
			if (i / 256 >= 128) ix = 1 - ix;
			tex[i] = c[ix];
		}
		gluBuild2DMipmaps(GL_TEXTURE_2D, 3, texSize, texSize, GL_RGB, GL_UNSIGNED_BYTE, &tex[0].r);
		delete[] tex;
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	glBindTexture(GL_TEXTURE_2D, viewportTextureID);
	return true;
}
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------