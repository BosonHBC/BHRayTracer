#include "materials.h"
#include "scene.h"
#include <math.h>
#include "lights.h"
extern Node rootNode;
void recursive(Node* root, Ray ray, HitInfo & outHit, bool &_bHit);

Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const {
	Color ambientColor = Color::Black();
	Color outColor = Color::Black();
	Vec3f vN = hInfo.N.GetNormalized();
	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	// Phi is dot product of View and Normal
	float cosPhi1 = vN.Dot(vV);

	for (auto it = lights.begin(); it != lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			// transform light
			Vec3f vL = (-(*it)->Direction(hInfo.p)).GetNormalized();
			// Theta is  dot product of Normal and light
			float cosTheta = vL.Dot(vN);

			// Reflection color
			{
				float R0 = pow((1 - ior / 1 + ior), 2);
				float RPhi = R0 + (1 - R0)* pow((1 - cosPhi1), 5);
				Color FresnelReflectionFactor = refraction * RPhi;
				Color reflectionFactor = reflection + FresnelReflectionFactor; //add Fresnel Reflection later
				if (!reflectionFactor.IsBlack()) {
					Ray reflectionRay;
					reflectionRay.p = hInfo.p;
					reflectionRay.dir = 2 * cosPhi1 * vN - vV;
					HitInfo reflHInfo;
					bool bReflectionHit;
					recursive(&rootNode, reflectionRay, reflHInfo, bReflectionHit);
					if (bReflectionHit && reflHInfo.node != nullptr && bounceCount > 0) {
						bounceCount--;
						outColor += reflectionFactor * reflHInfo.node->GetMaterial()->Shade(reflectionRay, reflHInfo, lights, bounceCount);
					}
				}
			}
			// Refraction color
			{
				if (!refraction.IsBlack()) {
					float sinPhi1 = sqrt(1 - cosPhi1 * cosPhi1);
					float sinPhi2 = sinPhi1 * 1.f / ior;
					float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);
					
					Vec3f vTn = -cosPhi2 * vN;
					Vec3f vTp = -sinPhi2 * (vV - cosPhi1 * vN);
					Vec3f vT = vTn + vTp;

					Ray refractionRay;
					refractionRay.p = hInfo.p;
					refractionRay.dir = vT;
					HitInfo refraHInfo;
					hInfo.node->GetNodeObj()->IntersectRay(refractionRay, refraHInfo, 1); // we need back face

				}
			}

			if (cosTheta < 0) {
				// from back side
				continue;
			}
			// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
			Vec3f vH = (vL + vV).GetNormalized();
			Color bdrf = diffuse  + specular * pow(vH.Dot(vN), glossiness) * 1/cosTheta;
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
