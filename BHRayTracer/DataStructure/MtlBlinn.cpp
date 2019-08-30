#include "materials.h"
#include "scene.h"
Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights) const {
	// fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
	Color outColor(0,0,0);
	Vec3f vN = hInfo.N;
	vN.Normalize();
	Vec3f vV = -ray.dir + hInfo.p;
	vV.Normalize();
	for (auto it = lights.begin(); it!= lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			Vec3f vL = -(*it)->Direction(hInfo.p);
			float cosTheta = vL.Dot(vN);
			if (cosTheta <= 0) {
				// from back side
				break;
			}
			Vec3f vH = (vL + vV);
			vH.Normalize();
			Color bdrf = diffuse * cosTheta + specular * vH.Dot(vN);
			outColor += bdrf;// +(*it)->Illuminate(hInfo.p, vN);
		}
		else {
			// it is ambient
			outColor += (*it)->Illuminate(hInfo.p, vN);
		}
	}


	return outColor;
}
