#include "materials.h"
#include "scene.h"
#include <math.h>
Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights) const {
	// fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
	Color outColor(0,0,0);
	Vec3f vN = hInfo.N.GetNormalized();

	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	for (auto it = lights.begin(); it!= lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			Vec3f vL = (-(*it)->Direction(hInfo.p)).GetNormalized();
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
			outColor += (*it)->Illuminate(hInfo.p, vN);
		}
	}
	//outColor.ClampMax();

	return outColor;
}
