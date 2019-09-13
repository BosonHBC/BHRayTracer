#include "materials.h"
#include "scene.h"
#include <math.h>
#include "lights.h"
#define Bias 0.001f
extern Node rootNode;
extern LightList lights;
void recursive(Node* root, Ray ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
void RefractionInternalRecursive(Node* root, const Node* myNode, Ray ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
cy::Color GetRefractionOutHitColor(const Ray& refractionRay_in, const float& sinPhi2_out, const Vec3f& vN_out, const Vec3f& vV_out, const float& RPhi, const Color& refraction);


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
			float R0 = pow((1 - ior) / (1 + ior), 2);
			float RPhi = R0 + (1 - R0)* pow((1 - cosPhi1), 5);
			{
				Color FresnelReflectionFactor = refraction * RPhi;
				Color reflectionFactor = reflection + FresnelReflectionFactor; //add Fresnel Reflection later
				if (!reflectionFactor.IsBlack()) {
					Ray reflectionRay;
					reflectionRay.p = hInfo.p;
					reflectionRay.dir = 2 * cosPhi1 * vN - vV;
					HitInfo reflHInfo;
					bool bReflectionHit = false;
					recursive(&rootNode, reflectionRay, reflHInfo, bReflectionHit, 0);
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
					vT.Normalize();

					Ray refractionRay_in;
					refractionRay_in.dir = vT;
					refractionRay_in.p = hInfo.p + vT* Bias;
					HitInfo refraHInfo_in;
					bool bRefraction_in_hit;
					//recursive(&rootNode, refractionRay_in, refraHInfo_in, bRefraction_in_hit, 0);
					RefractionInternalRecursive(&rootNode,  hInfo.node, refractionRay_in, refraHInfo_in, bRefraction_in_hit, 0);
					if (bRefraction_in_hit) {
						Vec3f vV_out = (hInfo.p - refraHInfo_in.p).GetNormalized();
						Vec3f vN_out = refraHInfo_in.N;
						float cosPhi1_out = vV_out.Dot(-vN_out);
						float sinPhi1_out = sqrt(1 - cosPhi1_out * cosPhi1_out);
						float sinPhi2_out = sinPhi1_out * ior;
						Color refractionColor = Color::Black();
						if (sinPhi2_out <= 1) {
							refractionColor = GetRefractionOutHitColor(refractionRay_in, sinPhi2_out, vN_out, vV_out, RPhi, refraction);
						}
						else {
/*
							// Total internal reflection
							int bounceCount = 3;
							while (bounceCount > 0)
							{
								Vec3f vR_internal = 2 * vV_out.Dot(-vN_out)*-vN_out - vV_out;
								Ray internalRay;
								internalRay.dir = vR_internal.GetNormalized();
								internalRay.p = refraHInfo_in.p + internalRay.dir * Bias;
								HitInfo internalReflectionHitInfo;
								bool b_InternalHit;
								RefractionInternalRecursive(&rootNode, hInfo.node, refractionRay_in, internalReflectionHitInfo, b_InternalHit, 0);
								if (b_InternalHit && internalReflectionHitInfo.node->GetNodeObj() != nullptr) {
									Vec3f vV_internal = (refraHInfo_in.p - internalReflectionHitInfo.p).GetNormalized();
									Vec3f vN_internal = internalReflectionHitInfo.N;
									float cosPhi1_internal = vV_internal.Dot(-vN_internal);
									float sinPhi1_internal = sqrt(1 - cosPhi1_internal * cosPhi1_internal);
									float sinPhi2_internal = ior * sinPhi1_internal;
									if (sinPhi2_internal <= 1) {
										Vec3f vV_internal_out = (internalRay.p - internalReflectionHitInfo.p).GetNormalized();
										refractionColor = GetRefractionOutHitColor(internalRay, sinPhi2_internal, internalReflectionHitInfo.N, vV_internal_out, RPhi, refraction, lights);
										break;
									}
									else {
										bounceCount--;
										vV_out = (internalRay.p - internalReflectionHitInfo.p).GetNormalized();
										vN_out = internalReflectionHitInfo.N;
									}
								}

							}*/
						}
						outColor += refractionColor;
					}
				}
			}

			if (cosTheta < 0) {
				// from back side
				continue;
			}
			// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
			Vec3f vH = (vL + vV).GetNormalized();
			Color bdrf = diffuse + specular * pow(vH.Dot(vN), glossiness) * 1 / cosTheta;
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

void RefractionInternalRecursive(Node* root, const Node* myNode,Ray ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/) {
	if (root->GetNumChild() <= 0) return;
	for (int i = 0; i < root->GetNumChild(); i++)
	{
		Ray transformedRay = root->GetChild(i)->ToNodeCoords(ray);
		if (root->GetChild(i)->GetNodeObj() == myNode->GetNodeObj()) {

			// transform ray to child coordinate
			if (root->GetChild(i)->GetNodeObj()->IntersectRay(transformedRay, outHit, hitSide))
			{
				outHit.node = root->GetChild(i);
				_bHit = true;
				root->GetChild(i)->FromNodeCoords(outHit);
			}
			break;
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


cy::Color GetRefractionOutHitColor(const Ray& refractionRay_in, const float& sinPhi2_out,const Vec3f& vN_out, const Vec3f& vV_out,const float& RPhi, const Color& refraction)
{
	float cosPhi2_out = sqrt(1 - sinPhi2_out * sinPhi2_out);
	Vec3f vTn_out = cosPhi2_out * vN_out;
	Vec3f xNV = vN_out.Cross(vV_out);
	Vec3f vTp_out = sinPhi2_out * vN_out.Cross(xNV).GetNormalized();
	Vec3f vT_out = vTp_out + vTn_out;

	Ray refractionRay_out;
	refractionRay_out.p = refractionRay_in.p;
	refractionRay_out.dir = vT_out.GetNormalized();
	HitInfo refraHinfo_out;
	bool bRefraction_out_Hit = false;
	recursive(&rootNode, refractionRay_out, refraHinfo_out, bRefraction_out_Hit, 0);
	if (bRefraction_out_Hit && refraHinfo_out.node != nullptr) {
		return (1 - RPhi)*refraction * refraHinfo_out.node->GetMaterial()->Shade(refractionRay_out, refraHinfo_out, lights, 0);
	}
}