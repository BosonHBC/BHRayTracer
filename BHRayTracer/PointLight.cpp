#include "lights.h"
#define GlossyShadow

#ifdef GlossyShadow
#define InitialSampleCount 1
#define MaxSampleCount 1
#endif // GlossyShadow
Vec3f GetSampleAlongNormal(const Vec3f& N, float R);

cy::Color PointLight::Illuminate(Vec3f const &p, Vec3f const &N) const
{
	if (size > 0) {
		// Adaptive Sample
		int shadowRayHitCount = 0;
		Vec3f centerDir = position - p;
		for (int i = 0; i < InitialSampleCount; ++i)
		{
			shadowRayHitCount += Shadow(Ray(p, GetSampleAlongNormal(centerDir, size)), 1);
		}
		if (shadowRayHitCount == 0) return Color(0, 0, 0);
		else if (shadowRayHitCount == InitialSampleCount) return intensity;
		else {
			// need more samples
			for (int i = 0; i < MaxSampleCount - InitialSampleCount; ++i)
			{
				shadowRayHitCount += Shadow(Ray(p, GetSampleAlongNormal(centerDir, size)), 1);
			}
			float percentage = (float)shadowRayHitCount / MaxSampleCount;
			return percentage * intensity;
		}
	}
	else
	return Shadow(Ray(p, position - p), 1) * intensity;
}
