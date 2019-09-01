#include "materials.h"
#include "scene.h"
#include <math.h>
#include "lights.h"
Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &_lights) const {
	// fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
	Color outColor = Color::Black();
	Vec3f vN = hInfo.N.GetNormalized();

	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	Color ambientColor = Color::Black();
	for (auto it = _lights.begin(); it != _lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			// transform light
			Vec3f vL = (-(*it)->Direction(hInfo.p)).GetNormalized();
			float cosTheta = vL.Dot(vN);
			if (cosTheta < 0) {
				// from back side
				continue;
			}
			Vec3f vH = (vL + vV).GetNormalized();
			Color bdrf = diffuse  + specular * pow(vH.Dot(vN), glossiness);
			outColor += bdrf * (*it)->Illuminate(hInfo.p, vN)* cosTheta;
		}
		else {
			// it is ambient
			ambientColor = diffuse * (*it)->Illuminate(hInfo.p, vN);
		}
	}
	//outColor.ClampMax();

	outColor += ambientColor;
	return outColor;
}
