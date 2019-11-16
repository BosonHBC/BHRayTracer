#include "materials.h"
#include "lights.h"
#include "scene.h"
#include <math.h>
#include <omp.h>

#define ENABLE_Reflection_And_Refraction
#define Bias 0.0001f
#define EulerN 2.7182818f
#define PI 3.14159265
#define OneOverPI 0.31830988f
#ifdef ENABLE_Reflection_And_Refraction
//#define ENABLE_Reflection
#define ENABLE_Refraction


#ifdef ENABLE_Refraction
#define ENABLE_InternalReflection
#endif
#endif // ENABLE_Reflection_And_Refraction

#define ENABLE_GI

#ifdef ENABLE_GI
#define GISampleCount 16
#endif // ENABLE_GI

#define ParallelVecDeterminance 0.001f

#define GlossyReflectionSampleCount 8
#define GlossyRefractionSampleCount 8

extern Node rootNode;
extern LightList lights;
extern TexturedColor environment;
extern float allLightIntensity;
float Rnd01() { return ((double)rand() / (RAND_MAX)); }

void recursive(Node* root, const Ray& ray, HitInfo & outHit, bool &_bHit, int hitSide /*= HIT_FRONT*/);
// get the ray from sphere inside to outside or doing internal reflection
Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness);
Color RefractionRecusive(const Color& refraction, const float& ior, const Vec3f& vT, const HitInfo& hInfo, const Vec3f& vN_new, const float& refractionGlossiness, const Color& absorption, int o_bounceCount);
cy::Color RefractionOut(const Ray& outRay, const Color& absorption, const Color& refraction, int o_bounceCount);
// Diffuse and Specular component
Color DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV);

Color GlobalIllumination(const TexturedColor& diffuse, const HitInfo& hInfo, const Vec3f& vN, int o_bounceCount, int i_bounceCount = 1);

// Reflection Component
Color Reflection(const Color& reflection, const float& cosPhi1, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int GIBounceCount, const float& reflectionGlossiness);
// Refraction component
Color Refraction(const Color& refraction, const Color& absorption, const float& cosPhi1, const float& ior, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int GIBounceCount, const float& refractionGlossiness);

// Path tracing functions
Color PathTracing_DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV);
cy::Color PathTracing_GlobalIllumination(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIbounceCount);

Vec3f GetRandomCrossingVector(const Vec3f& V);
Vec3f GetSampleAlongNormal(const Vec3f& N, float radius);
Vec3f GetSampleInSemiSphere(const Vec3f& N, float& o_theta);
Vec3f GetSampleAlongLightDirection(const Vec3f& N, float alhpa, float& o_theta);
cy::Vec3f GetSampleInLight(Light* light, const HitInfo& hInfo);
bool IsVec3Parallel(const Vec3f& i_lVec, const Vec3f& i_rVec);

bool s_debugTrace;

void PrintDebugColor(const char* i_colorName, const Color& i_color) {
	printf("%s(%f, %f, %f)\n", i_colorName, i_color.r, i_color.g, i_color.b);
}

Color MtlBlinn::Shade(Ray const &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount, int GIBounceCount) const {
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
	// ----------
	// Diffuse and Specular
	outColor += PathTracing_DiffuseNSpecular(diffuse, specular, glossiness, hInfo, vN, vV);

	// ----------
	// Global Illumination


#ifdef ENABLE_GI
	outColor += PathTracing_GlobalIllumination(diffuse, specular, glossiness, hInfo, vN, vV, bounceCount, GIBounceCount);
#endif // ENABLE_GI
	// ----------
	// Reflection and Refraction
#ifdef ENABLE_Reflection_And_Refraction
	{
		if (bounceCount > 0) {
			bounceCount--;
			float R0 = pow((1 - ior) / (1 + ior), 2);
			// Fresnel Reflection factor
			float fresnelReflectionFactor = R0 + (1 - R0)* pow((1 - cosPhi1), 5);
#ifdef ENABLE_Reflection
			// Reflection color
			outColor += Reflection(fresnelReflectionFactor * refraction.GetColor() + reflection.GetColor(), cosPhi1, hInfo, vN, vV, bounceCount, GIBounceCount, reflectionGlossiness);
#endif // ENABLE_REFLECTION

#ifdef ENABLE_Refraction
			// Refraction color
			outColor += Refraction((1 - fresnelReflectionFactor)*refraction.GetColor(), absorption, cosPhi1, ior, hInfo, vN, vV, bounceCount, GIBounceCount, refractionGlossiness);
#endif // ENABLE_RFRACTION
		}
	}
#endif // ENABLE_REFLEC_REFRAC
	// ----------

	if (isnan(outColor.r)) {

		printf("OutColor Has NaN! \n");
		outColor = Color::NANPurple();
	}
	return outColor;
}

Color DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV) {
	Color outColor = Color::Black();
	// if it is back face, render as black
	for (auto it = lights.begin(); it != lights.end(); ++it)
	{
		if (!(*it)->IsAmbient()) {
			// transform light
			Vec3f vL = (-(*it)->Direction(hInfo.p)).GetNormalized();
			// Theta is  dot product of Normal and light
			float cosTheta = vL.Dot(vN);

			if (cosTheta < -0.0001) {
				// from back side
				continue;
			}
			// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
			Vec3f vH = (vL + vV).GetNormalized();

			Color brdf = diffuse.Sample(hInfo.uvw, hInfo.duvw)  * cosTheta + specular.Sample(hInfo.uvw, hInfo.duvw) * pow(vH.Dot(vN), glossiness);
			outColor += brdf * (*it)->Illuminate(hInfo.p, vN);
		}
		else {
			// it is ambient
			//outColor += diffuse.Sample(hInfo.uvw, hInfo.duvw)  * (*it)->Illuminate(hInfo.p, vN);
		}
	}

	if (isnan(outColor.r)) {
		printf("Diffuse/Specular color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("DiffuseSpecular", outColor);
	return outColor;
}

cy::Vec3f GetSampleAlongLightDirection(const Vec3f& N, float glossiness, float& o_theta)
{
	float weightTheta = ACosSafe(pow(Rnd01(), 1.f / (glossiness + 1.f)));
	o_theta = weightTheta;
	// The radius of the circle
	float R = tan(weightTheta);
	float phi = ((double)rand() / (RAND_MAX)) * 2 * PI;
	// Random point in a circle
	float x = R * cos(phi);
	float y = R * sin(phi);

	Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
	Vec3f axis2 = axis1.Cross(N);

	Vec3f sampledN = N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
	return sampledN;
}


cy::Vec3f GetSampleInLight(Light* light, const HitInfo& hInfo)
{
	if (PointLight* pLight = dynamic_cast<PointLight*>(light)) {
		// it is point light, get a random point in the sphere sample

		// Uniform distribution of the radius 
		float R = sqrt(Rnd01()) * pLight->GetSize();
		float theta = ((double)rand() / (RAND_MAX)) * 2 * PI;
		// Random point in a circle
		float x = R * cos(theta);
		float y = R * sin(theta);
		// direction from hit point to light center
		Vec3f N = pLight->GetPosition() - hInfo.p;

		Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
		Vec3f axis2 = axis1.Cross(N);

		Vec3f sampledN = N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
		return sampledN.GetNormalized();
	}
	else
	{
		// Direction Light or other lights
		return -light->Direction(hInfo.p).GetNormalized();
	}
}


cy::Color PathTracing_DiffuseNSpecular(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV)
{
	Color outColor = Color::Black();
	float rnd = Rnd01();
	int i = 0;
	while (rnd > lights[i]->GetIntensity() / allLightIntensity && i < lights.size() - 1)
	{
		i++;
	}

	Light* light = lights[i];
	//Vec3f vL = orignalLightDir;
	// transform light
	float theta1 = 0;
	//Vec3f vL = (GetSampleAlongLightDirection(-light->Direction(hInfo.p), glossiness, theta1)).GetNormalized();
	//Vec3f vL = -light->Direction(hInfo.p);
	Vec3f vL = GetSampleInLight(light, hInfo);

	// Theta is  dot product of Normal and light
	float cosTheta = vL.Dot(vN);

	if (cosTheta <= 0) {
		// from back side
		return Color::Black();
	}
	// Diffuse & Specular  //  fs = kd + ks * vH.dot(vN) * 1/ Cos(theta)
	Vec3f vH = (vL + vV).GetNormalized();

	Color brdfXCosTheta =OneOverPI * (diffuse.Sample(hInfo.uvw, hInfo.duvw)  * cosTheta +(glossiness *0.5f + 1)*specular.Sample(hInfo.uvw, hInfo.duvw) * pow(vH.Dot(vN), glossiness));
	outColor += brdfXCosTheta * (light)->Illuminate(hInfo.p, vN);

	if (isnan(outColor.r)) {
		printf("Diffuse/Specular color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("DiffuseSpecular", outColor);
	return outColor;
}


#ifdef ENABLE_GI
Color GlobalIllumination(const TexturedColor& diffuse, const HitInfo& hInfo, const Vec3f& vN, int o_bounceCount, int i_bounceCount)
{
	if (diffuse.GetColor().IsBlack()) return Color::Black();
	i_bounceCount--;
	// Bound too many times
	if (i_bounceCount < 0)return Color::Black();

	Color outColor = Color::Black();
	Color GIColorSum = Color::Black();
	for (int i = 0; i < GISampleCount; ++i)
	{
		float theta1 = 0;
		Ray GIRay;
		GIRay.dir = GetSampleInSemiSphere(vN, theta1).GetNormalized();

		GIRay.p = hInfo.p + vN * Bias;
		// float cosTheta = vN.Dot(GIRay.dir);

		HitInfo reflHInfo = HitInfo();
		bool bReflectionHit = false;
		recursive(&rootNode, GIRay, reflHInfo, bReflectionHit, HIT_FRONT);
		if (bReflectionHit && reflHInfo.node != nullptr) {
			Color diffuseColor = diffuse.Sample(hInfo.uvw, hInfo.duvw);
			Color indirectColor = Color::Black();

			/** Mesh intersect situation */
			if (abs(reflHInfo.z) > Bias)
				indirectColor = reflHInfo.node->GetMaterial()->Shade(GIRay, reflHInfo, lights, o_bounceCount, i_bounceCount);
			GIColorSum += indirectColor * diffuseColor;
		}
		else {
			// doesn't bounce to anything
			Vec3f dir_norm = GIRay.dir;
			if (dir_norm.x == dir_norm.y && dir_norm.x == 0) {
				// nan situation
				GIColorSum += Color::NANPurple();
				printf("Environment nan\n");
			}
			else {
				GIColorSum += environment.SampleEnvironment(dir_norm)/* * cosTheta*/ * diffuse.Sample(hInfo.uvw, hInfo.duvw);
			}
		}
	}

	outColor += GIColorSum / GISampleCount;
	if (isnan(outColor.r)) {
		printf("GI color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("GI", outColor);
	return outColor;
}

cy::Color PathTracing_GlobalIllumination(const TexturedColor& diffuse, const TexturedColor& specular, const float& glossiness, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int i_GIbounceCount /*= 1*/)
{
	i_GIbounceCount--;
	// Bound too many times
	if (i_GIbounceCount < 0)return Color::Black();
	Color outColor = Color::Black();

	// diffuse GI
	float diffuseTheta = 0;
	Vec3f diffuseRayDir = GetSampleInSemiSphere(vN, diffuseTheta).GetNormalized();
	float p_diffuseTheta = sin(2 * diffuseTheta);


	// specular GI
	float specularTheta = 0;
	float cosvVvN = vN.Dot(vV);
	Vec3f vR = 2 * cosvVvN * vN - vV;
	Vec3f specualrRayDir = GetSampleAlongLightDirection(vR, glossiness, specularTheta);
	float p_specularTheta = pow(cos(specularTheta), glossiness);

	// Probability 
	float P_Diffuse = diffuse.GetColor().Gray() * p_diffuseTheta;
	float P_sum = P_Diffuse + specular.GetColor().Gray() * p_specularTheta;

	float  P_Diffuse_Norm = P_Diffuse / P_sum;

	float rnd = Rnd01();
	bool useSpecular =rnd >= P_Diffuse_Norm;
	Ray GIRay;
	GIRay.dir = useSpecular ? specualrRayDir : diffuseRayDir;

	GIRay.p = hInfo.p + vN * Bias;
	// float cosTheta = vN.Dot(GIRay.dir);

	HitInfo reflHInfo = HitInfo();
	bool bReflectionHit = false;
	recursive(&rootNode, GIRay, reflHInfo, bReflectionHit, HIT_FRONT);
	if (bReflectionHit && reflHInfo.node != nullptr) {
		Color indirectColor = Color::Black();

		/** Mesh intersect situation */
		if (abs(reflHInfo.z) > Bias)
			indirectColor = reflHInfo.node->GetMaterial()->Shade(GIRay, reflHInfo, lights, o_bounceCount, i_GIbounceCount);
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
			outColor += environment.SampleEnvironment(dir_norm) * (useSpecular ? specular : diffuse).Sample(hInfo.uvw, hInfo.duvw);
		}
	}

	return outColor;
}


#endif // ENABLE_GI


#ifdef ENABLE_Reflection
Color Reflection(const Color& reflection, const float& cosPhi1, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int GIBounceCount, const float& reflectionGlossiness) {
	Color reflectionColor = Color::Black();

	if (!reflection.IsBlack()) {

		Color reflectionColorSum = Color::Black();
		int sampleCount = reflectionGlossiness > 0 ? GlossyReflectionSampleCount : 1;
		for (int i = 0; i < sampleCount; ++i)
		{
			Vec3f rndN = reflectionGlossiness > 0 ? GetSampleAlongNormal(vN, reflectionGlossiness) : vN;
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
			recursive(&rootNode, reflectionRay, reflHInfo, bReflectionHit, HIT_FRONT);
			if (bReflectionHit && reflHInfo.node != nullptr) {
				reflectionColorSum += reflection * reflHInfo.node->GetMaterial()->Shade(reflectionRay, reflHInfo, lights, o_bounceCount, 0);
			}
			else {
				// doesn't bounce to anything
				Vec3f dir_norm = reflectionRay.dir.GetNormalized();
				reflectionColorSum += reflection * environment.SampleEnvironment(dir_norm);
			}
		}
		reflectionColor = reflectionColorSum / sampleCount;
	}
	if (isnan(reflectionColor.r)) {
		printf("Reflection color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("Reflection", reflectionColor);
	return reflectionColor;
}

#endif // ENABLE_REFLECTION

#ifdef ENABLE_Refraction

cy::Color RefractionRecusive(const Color& refraction, const float& ior, const Vec3f& vT, const HitInfo& hInfo, const Vec3f& vN_new, const float& refractionGlossiness, const Color& absorption, int o_bounceCount)
{
	Ray refractionRay_in;
	refractionRay_in.dir = vT;
	refractionRay_in.p = hInfo.p - vN_new * Bias;
	HitInfo refraHInfo_in = HitInfo();
	bool bRefractionInHit;
	recursive(&rootNode, refractionRay_in, refraHInfo_in, bRefractionInHit, HIT_BACK);
	if (bRefractionInHit && refraHInfo_in.node != nullptr) {
		bool bGoingOut;

		int sampleCount = refractionGlossiness > 0 ? GlossyRefractionSampleCount : 1;
		Color refractionColorSum_out = Color::Black();
		for (int i = 0; i < sampleCount; ++i)
		{
			Ray nextRay = HandleRayWhenRefractionRayOut(refractionRay_in, refraHInfo_in, ior, bGoingOut, refractionGlossiness);
			if (bGoingOut) {
				refractionColorSum_out += RefractionOut(nextRay, absorption, refraction, o_bounceCount);
			}
			else {
				if (o_bounceCount <= 0)
				{
					//refractionColorSum_out += Color(0,1,0);
					Ray bouncingCountZeroRay = Ray(refraHInfo_in.p + refraHInfo_in.N*Bias, refractionRay_in.dir);
					refractionColorSum_out += RefractionOut(bouncingCountZeroRay, absorption, refraction, o_bounceCount);
				}
				else {
					o_bounceCount--;
					refractionColorSum_out += RefractionRecusive(refraction, ior, nextRay.dir, refraHInfo_in, refraHInfo_in.N, refractionGlossiness, absorption, o_bounceCount);
				}
			}
		}

		return refractionColorSum_out / sampleCount;
	}
	else {
		printf("Internal Refraction ray doesn't hit anything\n");
		// Trace to outside
		return Color::NANPurple();
	}
}


cy::Color RefractionOut(const Ray& outRay, const Color& absorption, const Color& refraction, int o_bounceCount)
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
		outColor = refraction * absorptionFactor* refraHinfo_out.node->GetMaterial()->Shade(outRay, refraHinfo_out, lights, o_bounceCount, 0);
	}
	else {
		// refraction out hit doesn't hit anything
		outColor = refraction * environment.SampleEnvironment(outRay.dir);
	}
	return outColor;
}


Color Refraction(const Color& refraction, const Color& absorption, const float& cosPhi1, const float& ior, const HitInfo& hInfo, const Vec3f& vN, const Vec3f& vV, int o_bounceCount, int GIBounceCount, const float& refractionGlossiness)
{
	Color refractionColor = Color::Black();
	if (!refraction.IsBlack()) {

		Color refractionColorSum_in = Color::Black();
		int sampleCount = refractionGlossiness > 0 ? GlossyRefractionSampleCount : 1;
		for (int i = 0; i < sampleCount; ++i)
		{
			Vec3f vN_new = refractionGlossiness > 0 ? GetSampleAlongNormal(vN, refractionGlossiness).GetNormalized() : vN;
			float cosPhi1New = vN_new.Dot(vV);
			// handle floating point precision issues
			{
				if (cosPhi1New > 1) cosPhi1New = 1;
				if (cosPhi1New <= 0) continue;
			}
			float sinPhi1 = sqrt(1 - cosPhi1New * cosPhi1New);
			float sinPhi2 = sinPhi1 / ior;
			float cosPhi2 = sqrt(1 - sinPhi2 * sinPhi2);

			Vec3f vTn = -cosPhi2 * vN_new;
			Vec3f vNxV = vN_new.Cross(vV);
			Vec3f vTp = vN_new.Cross(vNxV).GetNormalized()*sinPhi2;
			Vec3f vT = vTn + vTp;

			refractionColorSum_in += RefractionRecusive(refraction, ior, vT, hInfo, vN_new, refractionGlossiness, absorption, o_bounceCount) / sampleCount;
		}
		refractionColor = refractionColorSum_in / sampleCount;
	}

	if (isnan(refractionColor.r)) {
		printf("Refraction color has nan! \n");
		return Color::NANPurple();
	}
	if (s_debugTrace)
		PrintDebugColor("Refraction", refractionColor);
	return refractionColor;
}

Ray HandleRayWhenRefractionRayOut(const Ray& inRay, const HitInfo& inRayHitInfo, const float& ior, bool& toOut, const float& refractionGlossiness) {
	Vec3f vN = refractionGlossiness > 0 ? GetSampleAlongNormal(inRayHitInfo.N, refractionGlossiness).GetNormalized() : inRayHitInfo.N; // to up

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
	float r = ((double)rand() / (RAND_MAX));
	// Uniform distribution
	r = sqrt(r) * R;
	float theta = ((double)rand() / (RAND_MAX)) * 2 * PI;
	// Random point in a circle
	float x = r * cos(theta);
	float y = r * sin(theta);

	Vec3f axis1 = GetRandomCrossingVector(N).Cross(N);
	Vec3f axis2 = axis1.Cross(N);

	Vec3f sampledN = N + axis1.GetNormalized() * x + axis2.GetNormalized() * y;
	return sampledN;
}

cy::Vec3f GetSampleInSemiSphere(const Vec3f& N, float& o_theta)
{
	Vec3f otherVec = Vec3f(1, 0, 0);
	if (IsVec3Parallel(N, otherVec)) {
		otherVec = Vec3f(0, 1, 0);
	}

	Vec3f axisY = (N.Cross(otherVec)).GetNormalized();
	Vec3f axisX = N.Cross(axisY);
	if (axisX.Length() > 1.0001f || axisX.Length() < 0.9999) {
		printf("not normalize\n");
	}

	// Uniform distribution, phi -> [0 , 2*PI)
	float phi = ((double)rand() / (RAND_MAX)) * 2 * PI;

	float rnd = ((double)rand() / (RAND_MAX));
	// Uniform distribution, theta -> [0 , PI/2)
	float theta = 0.5f * ACosSafe(1 - 2 * rnd);
	o_theta = theta;
	float sinTheta = sin(theta);

	Vec3f retVec = sinTheta * cos(phi) * axisX + sinTheta * sin(phi) * axisY + cos(theta) * N;
	if (N.Dot(retVec) < 0) {
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