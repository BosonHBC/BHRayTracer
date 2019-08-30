#include "materials.h"
#include "scene.h"
#include <math.h>
Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights) const {
	// fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
	Color outColor = Color::Black();
	Vec3f vN = hInfo.N.GetNormalized();

	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	Color ambientColor = Color::Black();
	for (auto it = lights.begin(); it!= lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {

			Vec3f trLight = (*it)->Direction(hInfo.p);
			trLight = hInfo.node->VectorTransformTo(trLight);
			Vec3f vL = (-trLight).GetNormalized();
			float cosTheta = vL.Dot(vN);
			if (cosTheta <= 0) {
				// from back side
				break;
			}
			Vec3f vH = (vL + vV).GetNormalized();
			Color bdrf = diffuse * cosTheta + specular * pow(vH.Dot(vN), glossiness);
			outColor += bdrf * (*it)->Illuminate(hInfo.p, vN);
		}
		else {
			// it is ambient
			ambientColor = diffuse *(*it)->Illuminate(hInfo.p, vN);
		}
	}
	//outColor.ClampMax();

	outColor += ambientColor;
	return outColor;
}
