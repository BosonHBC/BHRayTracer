#include "materials.h"
#include "scene.h"
#include <math.h>
#include <omp.h>

#define ENABLE_Reflection_And_Refraction

#ifdef ENABLE_Reflection_And_Refraction
#define Bias 0.001f
#define EulerN 2.7182818f
#define PI 3.14159265
#define ENABLE_Reflection
#define ENABLE_Refraction

#ifdef ENABLE_Refraction
//#define ENABLE_InternalReflection
#endif
#endif // ENABLE_Reflection_And_Refraction

#define GlossyReflectionSampleCount 16
#define GlossyRefractionSampleCount 16

extern Node rootNode;
extern LightList lights;
extern TexturedColor environment;
void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
void RefractionInternalRecursive(Node* root, const Node* myNode, Ray ray, HitInfo & outHit, bool &_bHit);
// get the ray from sphere inside to outside or doing internal reflection
Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness);
// Diffuse and Specular component
Color DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV);
// Reflection Component
Color Reflection(const Color& reflection, const float& cosPhi1, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int& o_bounceCount, const float& reflectionGlossiness);
// Refraction component
Color Refraction(const Color& refraction, const Color& absorption, const float& cosPhi1, const float& ior, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int& o_bounceCount, const float& refractionGlossiness);

Vec3f GetRandomCrossingVector(const Vec3f& V);
Vec3f GetSampleAlongNormal(const Vec3f& N, float radius);


Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const {
	Color outColor = Color::Black();
	Vec3f vN = hInfo.N.GetNormalized();
	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	// ----------
	// Diffuse and Specular
	outColor += DiffuseNSpecular(diffuse, specular, glossiness, hInfo, vN, vV);
	// ----------
	// Reflection and Refraction
#ifdef ENABLE_Reflection_And_Refraction
	{
		if (bounceCount > 0) {
			bounceCount--;
			// Phi is dot product of View and Normal
			float cosPhi1 = vN.Dot(vV);
			// handle floating point precision issues
			{
				if (cosPhi1 > 1) cosPhi1 = 1;
				if (cosPhi1 <= 0) return outColor;
			}
			float R0 = pow((1 - ior) / (1 + ior), 2);
			// Fresnel Reflection factor
			float fresnelReflectionFactor = R0 + (1 - R0)* pow((1 - cosPhi1), 5);
#ifdef ENABLE_Reflection
			// Reflection color
			outColor += Reflection(fresnelReflectionFactor * refraction.GetColor() + reflection.GetColor(), cosPhi1, hInfo, vN, vV, bounceCount, reflectionGlossiness);
#endif // ENABLE_REFLECTION

#ifdef ENABLE_Refraction
			// Refraction color
			outColor += Refraction((1 - fresnelReflectionFactor)*refraction.GetColor(), absorption, cosPhi1, ior, hInfo, vN, vV, bounceCount, refractionGlossiness);
#endif // ENABLE_RFRACTION
		}
	}
#endif // ENABLE_REFLEC_REFRAC
	// ----------

	if (isnan(outColor.r)) printf("OutColor Has NaN! \n");
	return outColor;
}

Color DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV) {
	Color outColor = Color::Black();
	for (auto it = lights.begin(); it != lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			// transform light
			Vec3f vL = (-(*it)->Direction(hInfo.p)).GetNormalized();
			// Theta is  dot product of Normal and light
			float cosTheta = vL.Dot(vN);

			if (cosTheta < 0) {
				// from back side
				continue;
			}
			// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
			Vec3f vH = (vL + vV).GetNormalized();
			Color bdrf = diffuse.Sample(hInfo.uvw, hInfo.duvw)  * cosTheta + specular.Sample(hInfo.uvw, hInfo.duvw) * pow(vH.Dot(vN), glossiness);
			outColor += bdrf * (*it)->Illuminate(hInfo.p, vN);

		}
		else {
			// it is ambient
			outColor += diffuse.Sample(hInfo.uvw, hInfo.duvw)  * (*it)->Illuminate(hInfo.p, vN);
		}
	}
	return outColor;
}

#ifdef ENABLE_Reflection
Color Reflection(const Color& reflection, const float& cosPhi1, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int& o_bounceCount, const float& reflectionGlossiness) {
	Color reflectionColor = Color::Black();
	if (!reflection.IsBlack()) {

		Color reflectionColorSum = Color::Black();
		int sampleCount = reflectionGlossiness > 0 ? GlossyReflectionSampleCount : 1;
		for (int i = 0; i < sampleCount; ++i)
		{
			Vec3f rndN = reflectionGlossiness > 0 ? GetSampleAlongNormal(vN, reflectionGlossiness * vN.Length()) : vN;
			float newCosPhi1 = rndN.Dot(vV);
			// handle floating point precision issues
			{
				if (newCosPhi1 > 1) newCosPhi1 = 1;
				if (newCosPhi1 <= 0) continue;
			}
			Ray reflectionRay;
			reflectionRay.dir = (2 * newCosPhi1 * rndN - vV).GetNormalized();
			reflectionRay.p = hInfo.p + rndN * Bias;

			HitInfo reflHInfo = HitInfo();
			bool bReflectionHit = false;
			recursive(&rootNode, reflectionRay, reflHInfo, bReflectionHit, 0);
			if (bReflectionHit && reflHInfo.node != nullptr) {
				reflectionColorSum += reflection * reflHInfo.node->GetMaterial()->Shade(reflectionRay, reflHInfo, lights, o_bounceCount);
			}
			else {
				// doesn't bounce to anything
				Vec3f dir_norm = reflectionRay.dir.GetNormalized();
				reflectionColorSum += reflection * environment.SampleEnvironment(dir_norm);
			}
		}
		reflectionColor = reflectionColorSum / sampleCount;
	}
	return reflectionColor;
}

#endif // ENABLE_REFLECTION

#ifdef ENABLE_Refraction
Color Refraction(const Color& refraction, const Color& absorption, const float& cosPhi1, const float& ior, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int& o_bounceCount, const float& refractionGlossiness)
{
	Color refractionColor = Color::Black();
	if (!refraction.IsBlack()) {
		float sinPhi1 = sqrt(1 - cosPhi1 * cosPhi1);
		float sinPhi2 = sinPhi1 / ior;
		float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);

		Color refractionColorSum_in = Color::Black();
		int sampleCount = refractionGlossiness > 0 ? GlossyRefractionSampleCount : 1;
		for (int i = 0; i < sampleCount; ++i)
		{
			Vec3f vN_new = GetSampleAlongNormal(vN, refractionGlossiness *vN.Length());
			Vec3f vTn = -cosPhi2 * vN_new;
			Vec3f vNxV = vN_new.Cross(vV);
			Vec3f vTp = vN_new.Cross(vNxV).GetNormalized()*sinPhi2;
			Vec3f vT = vTn + vTp;

			Ray refractionRay_in;
			refractionRay_in.dir = vT;
			refractionRay_in.p = hInfo.p - vN_new * Bias;
			HitInfo refraHInfo_in = HitInfo();
			bool bRefractionInHit;
			recursive(&rootNode, refractionRay_in, refraHInfo_in, bRefractionInHit, 1);
			if (bRefractionInHit) {
				bool bGoingOut;
				Color refractionColorSum_out = Color::Black();
				int sampleCount = refractionGlossiness > 0 ? GlossyRefractionSampleCount : 1;
				for (int i = 0; i < sampleCount; ++i)
				{
					Ray nextRay = HandleRayWhenRefractionRayOut(refractionRay_in, refraHInfo_in, ior, bGoingOut, refractionGlossiness);
					if (bGoingOut) {
						HitInfo refraHinfo_out = HitInfo();
						bool bRefraction_out_Hit = false;
						recursive(&rootNode, nextRay, refraHinfo_out, bRefraction_out_Hit, 0);
						if (bRefraction_out_Hit && refraHinfo_out.node != nullptr) {
							float absorptionFactorR = pow(EulerN, -absorption.r*refraHinfo_out.z);
							float absorptionFactorG = pow(EulerN, -absorption.g*refraHinfo_out.z);
							float absorptionFactorB = pow(EulerN, -absorption.b*refraHinfo_out.z);
							Color absorptionFactor(absorptionFactorR, absorptionFactorG, absorptionFactorB);
							refractionColorSum_out += refraction * absorptionFactor* refraHinfo_out.node->GetMaterial()->Shade(nextRay, refraHinfo_out, lights, o_bounceCount);
						}
						else {
							// refraction out hit doesn't hit anything
							refractionColorSum_out += refraction * environment.SampleEnvironment(nextRay.dir);
						}
					}
					else {
						// internal reflection
#ifdef ENABLE_InternalReflection
						int bounceCount = 3;
						Ray internalRay = nextRay;
						while (bounceCount > 0)
						{
							HitInfo internalHitInfo;
							bool bInternalHit;
							recursive(&rootNode, internalRay, internalHitInfo, bInternalHit, 1);
							if (bInternalHit) {
								bool bGoOut = false;
								Ray nextRay_internal = HandleRayWhenRefractionRayOut(internalRay, internalHitInfo, ior, bGoOut);
								if (bGoOut) {
									HitInfo refraHinfo_out;
									bool bRefraction_out_Hit = false;
									recursive(&rootNode, nextRay_internal, refraHinfo_out, bRefraction_out_Hit, 0);
									if (bRefraction_out_Hit && refraHinfo_out.node != nullptr) {

										refractionColor = refraction * refraHinfo_out.node->GetMaterial()->Shade(nextRay_internal, refraHinfo_out, lights, 0);
									}
									else {
										// refraction out hit doesn't hit anything
										refractionColor = refraction * environment.SampleEnvironment(nextRay.dir);
									}
									break;
								}
								else {
									internalRay = nextRay_internal;
								}
							}
							bounceCount--;
						}
#else 
						refractionColorSum_in += Color(0, 0, 0);
#endif // ENABLE_INTERNAL_REFLECTION

					}
				}
				refractionColorSum_in += refractionColorSum_out / sampleCount;
			}
		}
		refractionColor = refractionColorSum_in / sampleCount;
	}
	return refractionColor;
}

Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness) {
	Vec3f vN = refractionGlossiness > 0 ? GetSampleAlongNormal(inRayHitInfo.N, refractionGlossiness * inRayHitInfo.N.Length()) : inRayHitInfo.N; // to up

	Vec3f vV = -inRay.dir; // opposite to up

	float cosPhi1 = vV.Dot(-vN);
	float sinPhi1 = sqrt(1 - cosPhi1 * cosPhi1);
	float sinPhi2 = ior * sinPhi1;
	if (sinPhi2 <= 1) {
		// going out
		float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);
		Vec3f vTn = vN * cosPhi2;
		Vec3f vNxV = vN.Cross(vV);
		Vec3f vTp = vN.Cross(vNxV).GetNormalized() * sinPhi2;
		Vec3f vT = vTn + vTp;

		Ray outsideRay;
		outsideRay.dir = vT;
		outsideRay.p = inRayHitInfo.p + vN * Bias;
		toOut = true;
		return outsideRay;
	}
	else {
		// internal reflection
		Vec3f vR = (-2 * cosPhi1 * vN - vV);
		Ray internalRay;
		internalRay.dir = vR;
		internalRay.p = inRayHitInfo.p - vN * Bias;
		toOut = false;
		return internalRay;
	}
}
#ifdef ENABLE_InternalReflection
void RefractionInternalRecursive(Node* root, const Node* myNode, Ray ray, HitInfo & outHit, bool &_bHit) {
	if (root->GetNumChild() <= 0) return;
	for (int i = 0; i < root->GetNumChild(); i++)
	{
		Ray transformedRay = root->GetChild(i)->ToNodeCoords(ray);
		if (root->GetChild(i)->GetNodeObj() == myNode->GetNodeObj()) {

			// transform ray to child coordinate
			if (root->GetChild(i)->GetNodeObj()->IntersectRay(transformedRay, outHit, 0))
			{
				outHit.node = root->GetChild(i);
				_bHit = true;
				root->GetChild(i)->FromNodeCoords(outHit);
			}
		}
		RefractionInternalRecursive(root->GetChild(i), myNode, transformedRay, outHit, _bHit);
	}
	for (int i = 0; i < root->GetNumChild(); i++) {
		if (root->GetChild(i) == outHit.node) {
			root->FromNodeCoords(outHit);
			break;
		}
	}
}

#endif // ENABLE_INTERNAL_REFLECTION
#endif // ENABLE_RFRACTION
cy::Vec3f GetRandomCrossingVector(const Vec3f& V)
{
	Vec3f rndVec = Vec3f(0, 0, 1);
	// if the rndVec is not crossing with V, generate a new rand one
	while (V.Cross(rndVec).IsZero())
	{
		rndVec = Vec3f(((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)), ((double)rand() / (RAND_MAX)));
	}
	return rndVec;
}
cy::Vec3f GetSampleAlongNormal(const Vec3f& N, float R)
{
	float r = ((double)rand() / (RAND_MAX)) * R;
	// Uniform distribution
	r = sqrt(r) * R;
	float theta = ((double)rand() / (RAND_MAX)) * 2 * PI;
	// Random point in a circle
	float x = r * cos(theta);
	float y = r * sin(theta);

	Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
	Vec3f axis2 = axis1.Cross(N);

	return N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
}