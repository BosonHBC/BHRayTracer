#include "materials.h"
#include "lights.h"
#include "scene.h"
#include <math.h>
#include <omp.h>
#include "cyPhotonMap.h"


#define ENABLE_Reflection_And_Refraction
#define Bias 0.0001f
#define EulerN 2.7182818f
#define PI 3.14159265
#define OneOverPI 0.31830988f
#define ENABLE_Refraction
#ifdef ENABLE_Refraction
#define ENABLE_InternalReflection
#endif

#define ENABLE_GI

#define ParallelVecDeterminance 0.001f

#define UsePhotonMapping
#ifdef UsePhotonMapping
//#define USE_PhotonMap
/** Photon Map*/
#define Photon_AbsorbChance 0.3f
#define MAX_PhotonCountInArea 1000
#define MAX_Area 0.5
#define GIBounceCount 3
// -----------------------------------
#endif // UsePhotonMapping



extern Node rootNode;
extern LightList lights;
extern TexturedColor environment;
extern PhotonMap* photonMap;
extern PhotonMap* causticPhotonMap;
extern float allLightIntensity;

float Rnd01() {
	float rnd = ((double)rand() / (RAND_MAX));
	while (rnd == 0.0f || rnd == 1.0f)
	{
		rnd = ((double)rand() / (RAND_MAX));
	}
	return rnd;
}

void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
// get the ray from sphere inside to outside or doing internal reflection


// Path tracing functions
cy::Color PathTracing_DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int i_GIbounceCount);
cy::Color PathTracing_GlobalIllumination(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIbounceCount);
// Refraction
// -----------------
Color PathTracing_Refraction(const Color& refraction, const Color& absorption, const float& ior, const HitInfo& hInfo, const float cosPhi1, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIBounceCount, const float& refractionGlossiness);
Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness);
Color RefractionRecusive(const Color& refraction, const float& ior, const Vec3f& vT, const HitInfo& hInfo, const Vec3f& vN_new, const float& refractionGlossiness, const Color& absorption, int o_bounceCount, int i_GIbounceCount);
cy::Color RefractionOut(const Ray& outRay, const Color& absorption, const Color& refraction, int o_bounceCount, int i_GIbounceCount);
// -----------------

Vec3f GIUseSpecularDirOrDiffuseDir(bool& o_bUseSpecular, const Vec3f& vN, const Vec3f& vV, const float kd, const float ks, const float& glossiness);

const float GetKD(const TexturedColor& diffuse) { return Max(Max(diffuse.GetColor().r, diffuse.GetColor().g), diffuse.GetColor().b); }
const float GetKS(const TexturedColor& specular) { return Max(Max(specular.GetColor().r, specular.GetColor().g), specular.GetColor().b); }
Vec3f GetRandomCrossingVector(const Vec3f& V);
Vec3f GetSampleAlongNormal(const Vec3f& N, float radius);
Vec3f GetSampleInSemiSphere(const Vec3f& N, float& o_theta);
Vec3f GetSampleAlongLightDirection(const Vec3f& N, float alhpa, float& o_theta);
cy::Vec3f GetSampleInLight(const TexturedColor& diffuse, const TexturedColor& specular, Light* light, const HitInfo& hInfo, const float& glossiness);
bool IsVec3Parallel(const Vec3f& i_lVec, const Vec3f& i_rVec);

bool s_debugTrace;

void ClampColorToWhite(Color& o_color) {
	if (o_color.r > 1) o_color.r = 1.f;
	if (o_color.g > 1) o_color.g = 1.f;
	if (o_color.b > 1) o_color.b = 1.f;
}

void PrintDebugColor(const char* i_colorName, const Color& i_color) {
	printf("%s(%f, %f, %f)\n", i_colorName, i_color.r, i_color.g, i_color.b);
}

Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, int i_GIBounceCount) const {
	Color outColor = Color::Black();
	Vec3f vN = hInfo.N.GetNormalized();
	Vec3f vV = (ray.p - hInfo.p).GetNormalized();
	// Phi is dot product of View and Normal
	float cosPhi1 = vN.Dot(vV);
	// handle floating point precision issues
	{
		if (cosPhi1 > 1) cosPhi1 = 1;
		if (cosPhi1 <= 0) cosPhi1 = 0;
	}
	// ----------------
	TexturedColor newSpecular = specular;



#ifdef ENABLE_Refraction
	// Refraction
	float R0 = pow((1 - ior) / (1 + ior), 2);
	// Fresnel Reflection factor
	float fresnelReflectionFactor = R0 + (1 - R0)* pow((1 - cosPhi1), 5);
	Color fresnelSpecularColor = specular.GetColor() + fresnelReflectionFactor * refraction.GetColor();
	ClampColorToWhite(fresnelSpecularColor);
	// Set up fresnelReflection for specular
	newSpecular.SetColor(fresnelSpecularColor);
	newSpecular.SetTexture(const_cast<TextureMap*>(specular.GetTexture()));
	float refraction_glossiness = 0;
	if (glossiness > 50) refraction_glossiness = glossiness;
	outColor += PathTracing_Refraction((1 - fresnelReflectionFactor) *refraction.GetColor(), absorption, ior, hInfo, cosPhi1, vN, vV, bounceCount, i_GIBounceCount, refraction_glossiness);
	if (outColor.r >= 1 && outColor.g >= 1 && outColor.b >= 1) return outColor;
#endif // ENABLE_Refraction

#ifdef ENABLE_GI
	// ----------
	// Global Illumination
	outColor += PathTracing_GlobalIllumination(diffuse, newSpecular, glossiness, hInfo, vN, vV, bounceCount, i_GIBounceCount);
	if (outColor.r >= 1 && outColor.g >= 1 && outColor.b >= 1) return outColor;
#endif // ENABLE_GI

		// ----------
		// Diffuse and Specular
		outColor += PathTracing_DiffuseNSpecular(diffuse, newSpecular, glossiness, hInfo, vN, vV, i_GIBounceCount);
		if (outColor.r >= 1 && outColor.g >= 1 && outColor.b >= 1) return outColor;

	if (isnan(outColor.r)) {
		printf("OutColor Has NaN! \n");
		outColor = Color::NANPurple();
	}
	return outColor;
}

bool MtlBlinn::RandomPhotonBounce(Ray &r, Color &c, const HitInfo &hInfo) const
{
	// rnd for determine absorb, diffuse, specular, refraction 
	float rnd = Rnd01();

	Vec3f vN = hInfo.N.GetNormalized();
	Vec3f vV = -r.dir.GetNormalized();

	// if it is a transmissive material
	if (refraction.GetColor().Gray() > 0) {
		return false;
	}
	else {
		if (rnd < Photon_AbsorbChance) {
			// Not going to generate next bounce
			return false;
		}

		// use kd & ks to determine use diffuse or specular
		bool useSpecular;
		float p_Diff = 0;
		float p_Spec = 0;
		float p_Absorb = Photon_AbsorbChance;
		Vec3f dir;
		{
			// diffuse GI
			float diffuseTheta = 0;
			Vec3f diffuseRayDir = GetSampleInSemiSphere(vN, diffuseTheta).GetNormalized();
			float p_diffuseTheta = /*2* PI * */sin(2 * diffuseTheta);

			// specular GI
			float specularTheta = 0;
			float cosvVvN = vN.Dot(vV);
			Vec3f vR = 2 * cosvVvN * vN - vV;
			Vec3f specualrRayDir = GetSampleAlongLightDirection(vR, glossiness, specularTheta);
			float p_specularTheta = /*2 * OneOverPI *(glossiness + 2)**/pow(cos(specularTheta), glossiness);

			// Probability 
			float P_Diffuse = GetKD() * p_diffuseTheta;
			float P_sum = P_Diffuse + GetKS() * p_specularTheta;

			p_Diff = (P_Diffuse / P_sum) * (1 - Photon_AbsorbChance) + Photon_AbsorbChance;
			p_Spec = (1 - p_Diff)  * (1 - Photon_AbsorbChance) + Photon_AbsorbChance;

			useSpecular = rnd >= p_Diff;

			dir = useSpecular ? specualrRayDir : diffuseRayDir;
		}
		Color kdf = diffuse.GetColor() / p_Diff;
		Color ksf = specular.GetColor() / p_Spec;

		Color newC = c * (useSpecular ? ksf : kdf);
		if (isnan(newC.r)) {
			printf("useSpecular: %s, c: (%f, %f, %f) ,kd: %f, pd: %f, ks: %f, ps: %f\n", useSpecular ? "true" : "false", c.r, c.g, c.b, GetKD(), p_Diff, GetKS(), p_Spec);
		}
		c = newC;

		r.dir = dir;
		r.p = hInfo.p + hInfo.N * Bias;
	}

	return true;
}
bool MtlBlinn::RandomPhotonBounceForCaustic(Ray &r, Color &c, const HitInfo &hInfo) const
{
	// rnd for determine absorb, diffuse, specular, refraction 
	float rnd = Rnd01();

	Vec3f vN = hInfo.N.GetNormalized();
	Vec3f vV = -r.dir.GetNormalized();

	// if it is a transmissive material
	if (refraction.GetColor().Gray() > 0) {
		// let this ray walk through this material

		float cosPhi1 = vN.Dot(vV);
		float sinPhi1 = sqrt(1 - cosPhi1 * cosPhi1);
		float sinPhi2 = sinPhi1 / ior;
		float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);

		Vec3f vTn = -cosPhi2 * vN;
		Vec3f vNxV = vN.Cross(vV);
		Vec3f vTp = vN.Cross(vNxV).GetNormalized()*sinPhi2;
		// perfect refraction direction
		Vec3f vT = vTn + vTp;

		Ray refractionRay_in;
		refractionRay_in.dir = vT;
		refractionRay_in.p = hInfo.p - vN * Bias;
		HitInfo refraHInfo_in = HitInfo();
		bool bRefractionInHit;
		recursive(&rootNode, refractionRay_in, refraHInfo_in, bRefractionInHit, HIT_BACK);

		if (bRefractionInHit && refraHInfo_in.node != nullptr) {
			bool bGoingOut;
			Ray nextRay = HandleRayWhenRefractionRayOut(refractionRay_in, refraHInfo_in, ior, bGoingOut, refractionGlossiness);
			if (bGoingOut) {
				r = nextRay;
				return true;
			}
			else {
				// This photon is not going to bounce
				return false;
			}
		}
		else {
			// This photon is not going to bounce
			return false;
		}
	}
	else {
		if (rnd < Photon_AbsorbChance) {
			// Not going to generate next bounce
			return false;
		}

		// use kd & ks to determine use diffuse or specular
		bool useSpecular;
		float p_Diff = 0;
		float p_Spec = 0;
		float p_Absorb = Photon_AbsorbChance;
		Vec3f dir;
		{
			// diffuse GI
			float diffuseTheta = 0;
			Vec3f diffuseRayDir = GetSampleInSemiSphere(vN, diffuseTheta).GetNormalized();
			float p_diffuseTheta = /*2* PI * */sin(2 * diffuseTheta);

			// specular GI
			float specularTheta = 0;
			float cosvVvN = vN.Dot(vV);
			Vec3f vR = 2 * cosvVvN * vN - vV;
			Vec3f specualrRayDir = GetSampleAlongLightDirection(vR, glossiness, specularTheta);
			float p_specularTheta = /*2 * OneOverPI *(glossiness + 2)**/pow(cos(specularTheta), glossiness);

			// Probability 
			float P_Diffuse = GetKD() * p_diffuseTheta;
			float P_sum = P_Diffuse + GetKS() * p_specularTheta;

			p_Diff = (P_Diffuse / P_sum) * (1 - Photon_AbsorbChance) + Photon_AbsorbChance;
			p_Spec = (1 - p_Diff)  * (1 - Photon_AbsorbChance) + Photon_AbsorbChance;

			useSpecular = rnd >= p_Diff;

			// if it is diffuse bounce stop
			if (!useSpecular) return false;

			dir = specualrRayDir;
		}

		Color ksf = specular.GetColor() / p_Spec;

		Color newC = c * ksf;
		if (isnan(newC.r)) {
			printf("useSpecular: %s, c: (%f, %f, %f) ,kd: %f, pd: %f, ks: %f, ps: %f\n", useSpecular ? "true" : "false", c.r, c.g, c.b, GetKD(), p_Diff, GetKS(), p_Spec);
		}
		c = newC;

		r.dir = dir;
		r.p = hInfo.p + hInfo.N * Bias;
	}

	return true;
}
cy::Color PathTracing_DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int i_GIbounceCount)
{
	Color outColor = Color::Black();
		// Direct shading
		{
			float rnd = Rnd01();
			int i = 0;
			while (rnd > lights[i]->GetIntensity() / allLightIntensity && i < lights.size() - 1)
			{
				i++;
			}

			Light* light = lights[i];
			Vec3f vL = GetSampleInLight(diffuse, specular, light, hInfo, glossiness);
			// Theta is  dot product of Normal and light
			float cosTheta = vL.Dot(vN);

			if (cosTheta > 0) {
				// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
				Vec3f vH = (vL + vV).GetNormalized();
				Color irrad = (light)->Illuminate(hInfo.p, vN);
				Color brdfXCosTheta = diffuse.Sample(hInfo.uvw, hInfo.duvw)  * cosTheta + specular.Sample(hInfo.uvw, hInfo.duvw) * pow(vH.Dot(vN), glossiness);
				outColor += irrad * brdfXCosTheta;
			}
		}
#ifdef USE_PhotonMap
		//Adding caustics for direct color
		{
			Vec3f vL;
			Color causticIndirectIrrad;
			causticPhotonMap->EstimateIrradiance<MAX_PhotonCountInArea>(causticIndirectIrrad, vL, MAX_Area, hInfo.p, &hInfo.N);
			float cosTheta = -vL.Dot(vN);
			if (cosTheta > 0) {
				Vec3f vH = (vL + vV).GetNormalized();
				Color brdf = diffuse.Sample(hInfo.uvw, hInfo.duvw) + specular.Sample(hInfo.uvw, hInfo.duvw) * pow(vH.Dot(vN), glossiness) / cosTheta;
				outColor += brdf * causticIndirectIrrad;
			}
		}
#endif
	ClampColorToWhite(outColor);
	if (isnan(outColor.r)) {
		printf("Diffuse/Specular color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("DiffuseSpecular", outColor);
	return outColor;
}


cy::Vec3f GIUseSpecularDirOrDiffuseDir(bool& o_bUseSpecular, const Vec3f& vN, const Vec3f& vV, const float kd, const float ks, const float& glossiness)
{
	// diffuse GI
	float diffuseTheta = 0;
	Vec3f diffuseRayDir = GetSampleInSemiSphere(vN, diffuseTheta).GetNormalized();
	float p_diffuseTheta = /*2* PI * */sin(2 * diffuseTheta);

	// specular GI
	float specularTheta = 0;
	float cosvVvN = vN.Dot(vV);
	Vec3f vR = 2 * cosvVvN * vN - vV;
	Vec3f specualrRayDir = GetSampleAlongLightDirection(vR, glossiness, specularTheta);
	float p_specularTheta = /*2 * OneOverPI *(glossiness + 2)**/pow(cos(specularTheta), glossiness);

	// Probability 
	float P_Diffuse = kd * p_diffuseTheta;
	float P_sum = P_Diffuse + ks * p_specularTheta;

	float  P_Diffuse_Norm = P_Diffuse / P_sum;

	float rnd = Rnd01();
	o_bUseSpecular = rnd >= P_Diffuse_Norm;

	return o_bUseSpecular ? specualrRayDir : diffuseRayDir;
}


#ifdef ENABLE_GI

cy::Color PathTracing_GlobalIllumination(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIbounceCount /*= 1*/)
{
	// Bound too many times
	if (i_GIbounceCount < 0)return Color::Black();
	Color outColor = Color::Black();
	{
		bool useSpecular;
		Ray GIRay;
		GIRay.dir = GIUseSpecularDirOrDiffuseDir(useSpecular, vN, vV, GetKD(diffuse), GetKS(specular), glossiness);
		GIRay.p = hInfo.p + vN * Bias;

		HitInfo reflHInfo = HitInfo();
		bool bReflectionHit = false;
		recursive(&rootNode, GIRay, reflHInfo, bReflectionHit, HIT_FRONT);
		if (bReflectionHit && reflHInfo.node != nullptr) {
			Color indirectColor = Color::Black();

			/** Mesh intersect situation */
			if (abs(reflHInfo.z) > Bias) {

				indirectColor = reflHInfo.node->GetMaterial()->Shade(GIRay, reflHInfo, lights, o_bounceCount, i_GIbounceCount - 1);
			}

			outColor += indirectColor * (useSpecular ? specular : diffuse).Sample(hInfo.uvw, hInfo.duvw);
		}
		else {
			// doesn't bounce to anything
			Vec3f dir_norm = GIRay.dir;
			if (dir_norm.x == dir_norm.y && dir_norm.x == 0) {
				// nan situation
				outColor += Color::NANPurple();
				printf("Environment nan\n");
			}
			else {
				 Color envColor = environment.SampleEnvironment(dir_norm) * (useSpecular ? specular : diffuse).Sample(hInfo.uvw, hInfo.duvw);
				 if (isnan(envColor.r) || isnan(envColor.g) || isnan(envColor.b)) {
				 }
				 else {
					 outColor += envColor;
				 }
			}
		}
	}

	if (isnan(outColor.r)) {
		printf("GI color has nan! \n");
		return Color::NANPurple();
	}
	ClampColorToWhite(outColor);
	return outColor;
}
#endif // ENABLE_GI


cy::Color PathTracing_Refraction(const Color& refraction, const Color& absorption, const float& ior, const HitInfo& hInfo, const float cosPhi1, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIBounceCount, const float& refractionGlossiness)
{
	Color refractionColor = Color::Black();
	if (o_bounceCount <= 0) return refractionColor;

	if (!refraction.IsBlack()) {

		float sinPhi1 = sqrt(1 - cosPhi1 * cosPhi1);
		float sinPhi2 = sinPhi1 / ior;
		float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);

		Vec3f vTn = -cosPhi2 * vN;
		Vec3f vNxV = vN.Cross(vV);
		Vec3f vTp = vN.Cross(vNxV).GetNormalized()*sinPhi2;
		// perfect refraction direction
		Vec3f vT = vTn + vTp;
		Vec3f vT_sampled = vT.GetNormalized();
		if (refractionGlossiness > 0)
		{
			// Sample a new refraction direction
			float dotSign = 0;
			// if the new vT is going out(sign of dot(vT, vN) > 0), generate a new one
			while (dotSign >= 0)
			{
				float theta = 0;
				vT_sampled = GetSampleAlongLightDirection(vT, refractionGlossiness, theta);
				dotSign = vT_sampled.Dot(vN);
			}
		}

		refractionColor = RefractionRecusive(refraction, ior, vT_sampled.GetNormalized(), hInfo, vN, refractionGlossiness, absorption, o_bounceCount, i_GIBounceCount);
	}

	ClampColorToWhite(refractionColor);

	return refractionColor;
}


cy::Color RefractionRecusive(const Color& refraction, const float& ior, const Vec3f& vT, const HitInfo& hInfo, const Vec3f& vN, const float& refractionGlossiness, const Color& absorption, int o_bounceCount, int i_GIbounceCount)
{
	Ray refractionRay_in;
	refractionRay_in.dir = vT;
	refractionRay_in.p = hInfo.p - vN * Bias;
	HitInfo refraHInfo_in = HitInfo();
	bool bRefractionInHit;
	recursive(&rootNode, refractionRay_in, refraHInfo_in, bRefractionInHit, HIT_FRONT_AND_BACK);
	if (bRefractionInHit && refraHInfo_in.node != nullptr) {
		bool bGoingOut;
		Color refractionColor = Color::Black();

		Ray nextRay = HandleRayWhenRefractionRayOut(refractionRay_in, refraHInfo_in, ior, bGoingOut, refractionGlossiness);
		if (bGoingOut) {
			refractionColor = RefractionOut(nextRay, absorption, refraction, o_bounceCount, i_GIbounceCount);
		}
		else {
			// Total internal reflection
			if (o_bounceCount <= 0)
			{
				// When the bounce count is not enough
				refractionColor = Color::Black();
			}
			else {
				o_bounceCount--;
				refractionColor = RefractionRecusive(refraction, ior, nextRay.dir, refraHInfo_in, refraHInfo_in.N, refractionGlossiness, absorption, o_bounceCount, i_GIbounceCount);
			}
		}

		ClampColorToWhite(refractionColor);
		return refractionColor;
	}
	else {
		// Trace to outside
		//printf("Ray did not hit anything");
		return Color::NANPurple();
	}
}

cy::Color RefractionOut(const Ray& outRay, const Color& absorption, const Color& refraction, int o_bounceCount, int i_GIbounceCount)
{
	Color outColor = Color::Black();
	HitInfo refraHinfo_out = HitInfo();
	bool bRefraction_out_Hit = false;
	recursive(&rootNode, outRay, refraHinfo_out, bRefraction_out_Hit, HIT_FRONT);
	if (bRefraction_out_Hit && refraHinfo_out.node != nullptr) {

		float absorptionFactorR = pow(EulerN, -absorption.r*refraHinfo_out.z);
		float absorptionFactorG = pow(EulerN, -absorption.g*refraHinfo_out.z);
		float absorptionFactorB = pow(EulerN, -absorption.b*refraHinfo_out.z);
		Color absorptionFactor(absorptionFactorR, absorptionFactorG, absorptionFactorB);
		outColor = refraction * absorptionFactor* refraHinfo_out.node->GetMaterial()->Shade(outRay, refraHinfo_out, lights, o_bounceCount, i_GIbounceCount - 1);
	}
	else {
		// refraction out hit doesn't hit anything
		outColor = refraction * environment.SampleEnvironment(outRay.dir);
	}
	ClampColorToWhite(outColor);
	return outColor;
}

Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness) {
	Vec3f vN = inRayHitInfo.N; // to up

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
		// perfect out refraction ray
		Vec3f vT = vTn + vTp;

		Vec3f vT_sampled = vT.GetNormalized();
		if (refractionGlossiness > 0)
		{
			// Sample a new refraction direction
			float dotSign = 0;
			// if the new vT is going out(sign of dot(vT, vN) < 0), generate a new one
			while (dotSign <= 0)
			{
				float theta = 0;
				vT_sampled = GetSampleAlongLightDirection(vT, refractionGlossiness, theta);
				dotSign = vT_sampled.Dot(vN);
			}
		}

		Ray outsideRay;
		outsideRay.dir = vT_sampled.GetNormalized();
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

cy::Vec3f GetRandomCrossingVector(const Vec3f& V)
{
	Vec3f rndVec = Vec3f(0, 0, 1);
	// if the rndVec is not crossing with V, generate a new rand one
	while (V.Cross(rndVec).IsZero())
	{
		rndVec = Vec3f(Rnd01(), Rnd01(), Rnd01());
	}
	return rndVec;
}

cy::Vec3f GetSampleAlongNormal(const Vec3f& N, float R)
{
	float r = Rnd01();
	// Uniform distribution
	r = sqrt(r) * R;
	float theta = Rnd01() * 2 * PI;
	// Random point in a circle
	float x = r * cos(theta);
	float y = r * sin(theta);

	Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
	Vec3f axis2 = axis1.Cross(N);

	Vec3f sampledN = N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
	return sampledN;
}

cy::Vec3f GetSampleAlongLightDirection(const Vec3f& N, float glossiness, float& o_theta)
{
	float weightTheta = ACosSafe(pow(Rnd01(), 1.f / (glossiness + 1.f)));
	o_theta = weightTheta;
	// The radius of the circle
	float R = tan(weightTheta);
	float phi = Rnd01() * 2 * PI;
	// Random point in a circle
	float x = R * cos(phi);
	float y = R * sin(phi);

	Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
	Vec3f axis2 = axis1.Cross(N);

	Vec3f sampledN = N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
	return sampledN;
}

cy::Vec3f GetSampleInLight(const TexturedColor& diffuse, const TexturedColor& specular, Light* light, const HitInfo& hInfo, const float& glossiness)
{
	if (PointLight* pLight = dynamic_cast<PointLight*>(light)) {
		// it is point light, get a random point in the sphere sample
		float kd = GetKD(diffuse);
		float ks = GetKS(specular);
		float p_diffuse = 0;
		Vec3f diffuse_vL;
		float p_specular = 0;
		Vec3f specular_vL;
		// direction from hit point to light center
		Vec3f vL = pLight->GetPosition() - hInfo.p;

		// Sample along the light direction
		{
			float diffuseTheta = 0;
			diffuse_vL = GetSampleAlongLightDirection(vL.GetNormalized(), glossiness, diffuseTheta);

			p_diffuse = /*2 * OneOverPI *(glossiness + 2) **/ pow(cos(diffuseTheta), glossiness);
		}

		if (ks == 0 && kd != 0) return diffuse_vL.GetNormalized();

		// Sample in the circle of light 
		{
			float r = Rnd01();
			float R = sqrt(r) * pLight->GetSize();
			float specularTheta = Rnd01() * 2 * PI;
			// Random point in a circle
			float x = R * cos(specularTheta);
			float y = R * sin(specularTheta);

			Vec3f axis1 = GetRandomCrossingVector(vL).Cross(vL);
			Vec3f axis2 = axis1.Cross(vL);

			specular_vL = vL + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
			p_specular = 2 * r / (R*R);
		}

		if (ks != 0 && kd == 0) return specular_vL.GetNormalized();

		// Probability 

		float P_Diffuse = kd * p_diffuse;
		float P_Specular = ks * p_specular;
		float P_sum = P_Diffuse + P_Specular;

		float  P_Diffuse_Norm = P_Diffuse / P_sum;
		float rnd = Rnd01();
		bool useSpecular = rnd >= P_Diffuse_Norm;

		return useSpecular ? specular_vL.GetNormalized() : diffuse_vL.GetNormalized();
	}
	else
	{
		// Direction Light or other lights
		return -light->Direction(hInfo.p).GetNormalized();
	}
}

cy::Vec3f GetSampleInSemiSphere(const Vec3f& N, float& o_theta)
{
	Vec3f axisY = (N.Cross(GetRandomCrossingVector(N))).GetNormalized();
	Vec3f axisX = N.Cross(axisY);

	// Uniform distribution, phi -> [0 , 2*PI)
	float phi = Rnd01() * 2 * PI;

	float rnd = Rnd01();
	// Uniform distribution, theta -> [0 , PI/2)
	float theta = 0.5f * ACosSafe(1 - 2 * rnd);
	o_theta = theta;
	float sinTheta = sin(theta);

	Vec3f retVec = sinTheta * cos(phi) * axisX + sinTheta * sin(phi) * axisY + cos(theta) * N;
	if (N.Dot(retVec) <= 0) {
		return GetSampleInSemiSphere(N, o_theta);
	}
	return retVec;
}

bool IsVec3Parallel(const Vec3f& i_lVec, const Vec3f& i_rVec)
{
	float l_len = i_lVec.Length();
	float r_len = i_rVec.Length();
	if (l_len == 0 && r_len == 0) return true;

	float cosTheta = i_lVec.Dot(i_rVec) / (l_len * r_len);
	if (abs(cosTheta) > 1 - ParallelVecDeterminance)
		return true;
	else
		return false;
}