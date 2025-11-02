#ifndef FORCES_H
#define FORCES_H

#include <vector>
#include "particle.h"

class Force
{
public:
    Force(void) {}
    virtual ~Force(void) {}

    virtual void apply() = 0;

    void addInfluencedParticle(Particle* p) {
        particles.push_back(p);
    }

    void setInfluencedParticles(const std::vector<Particle*>& vparticles) {
        particles = vparticles;
    }

    void clearInfluencedParticles() {
        particles.clear();
    }

    std::vector<Particle*> getInfluencedParticles() const {
        return particles;
    }

protected:
    std::vector<Particle*>	particles;
};


class ForceConstAcceleration : public Force
{
public:
    ForceConstAcceleration() { acceleration = Vec3(0,0,0); }
    ForceConstAcceleration(const Vec3& a) { acceleration = a; }
    virtual ~ForceConstAcceleration() {}

    virtual void apply();

    void setAcceleration(const Vec3& a) { acceleration = a; }
    Vec3 getAcceleration() const { return acceleration; }

protected:
    Vec3 acceleration;
};


class ForceDrag : public Force
{
public:
    ForceDrag() { klinear = kquadratic = 0; }
    ForceDrag(double k1, double k2) { klinear = k1; kquadratic = k2; }
    virtual ~ForceDrag() {}

    virtual void apply();

    void setDragCoefficients(double k1, double k2) { klinear = k1, kquadratic = k2; }
    double getLinearCoefficient() const { return klinear; }
    double getQuadraticCoefficient() const { return kquadratic; }

protected:
    double klinear, kquadratic;
};


class ForceSpring : public Force
{
public:
    ForceSpring() { ks = kd = 0; }
    ForceSpring(Particle* p1, Particle* p2, double L, double ks, double kd) {
        this->L = L; this->ks = ks; this->kd = kd;
        particles.push_back(p1);
        particles.push_back(p2);
    }
    virtual ~ForceSpring() {}

    virtual void apply();

    void setParticlePair(Particle* p1, Particle* p2) {
        particles.clear();
        particles.push_back(p1);
        particles.push_back(p2);
    }
    Particle* getParticle1() const { return particles[0]; }
    Particle* getParticle2() const { return particles[1]; }

    void setRestLength(double l) { L = l; }
    double getRestLength() const { return L; }

    void setSpringConstant(double k) { ks = k; }
    double getSpringConstant() const { return ks; }

    void setDampingCoeff(double k) { kd = k; }
    double getDampingCoeff() const { return kd; }

protected:
    double L = 0;   // resting length
    double ks = 0;  // spring coeff
    double kd = 0;  // damping coeff
};

class ForceGravitation : public Force
{
public:
    ForceGravitation() { attractor=nullptr; }
    ForceGravitation(const Particle* p, double k) { attractor = p, G = k; }
    virtual ~ForceGravitation() {}

    virtual void apply();

    void setAttractor(const Particle* p) { attractor = p; }
    const Particle* getAttractor() const { return attractor; }
    void setConstant(double k) { G = k; }
    void setSmoothingFactors(double sa, double sb) { a = sa; b = sb; }
    double getConstant() const { return G; }

protected:
    const Particle* attractor;
    double G = 6.6743e-11; // gravitational constant
    double a = 1, b = 1;
};

// Boids algorithm

class ParticleSystem;

class ForceSPH : public Force
{
public:
    ForceSPH(ParticleSystem* system) {
        this->system = system;
        h = 0;
    }
    ForceSPH(ParticleSystem* system, double h) {
        this->system = system;
        setH(h);
    }

    virtual void apply();

    void setH(double h);

    void setRestDensity(double p0) { this->p0 = p0; }
    void setViscosityConstant(double mu) { this->mu = mu; }

protected:
    double kernelPoly6(double rSqr) const;
    Vec3 kernelGradientSpiky(Vec3 r) const;
    double kernelLaplacianViscosity(double r) const;

    double computePressure(double density);

protected:
    ParticleSystem* system;
    double h = 0;
    double poly6Factor = 0;
    double spikyGradFactor = 0;
    double viscosityLaplacianFactor = 0;

    std::vector<double> densities, pressures;

    double p0;
    double mu;
};

inline void ForceSPH::setH(double h) {
    this->h = h;
    this->poly6Factor = 315.0 / (64.0 * M_PI * std::pow(h, 9));
    this->spikyGradFactor = 45.0 / (M_PI * std::pow(h, 6));
    this->viscosityLaplacianFactor = 45.0 / (M_PI * std::pow(h, 5));
}

// Wind force (localized fan/ventilator)
class ForceWind : public Force
{
public:
    ForceWind();
    
    // Create wind with center, cross-section size (width x height), flow direction, length, and strength
    // center: center of the outlet plane where air comes out
    // width, height: cross-section dimensions of the air flow (perpendicular to direction)
    // dir: flow direction (will be snapped to nearest axis)
    // length: how far the air reaches from the outlet
    // strength: acceleration magnitude
    ForceWind(const Vec3& center, double width, double height, const Vec3& dir, double length, double strength);
    
    virtual ~ForceWind() {}

    virtual void apply();

    // Set center (the face where air comes out)
    void setCenter(const Vec3& c) { center = c; updateInternals(); }
    
    // Set cross-section dimensions (width x height perpendicular to flow)
    void setDimensions(double w, double h) { width = w; height = h; updateInternals(); }
    
    // Set flow direction (will be snapped to nearest axis: ±X, ±Y, ±Z)
    void setDirection(const Vec3& d);
    
    // Set how far the air reaches
    void setLength(double l) { length = l; updateInternals(); }
    
    void setStrength(double s) { strength = s; }
    void setFrequency(double f) { frequency = f; }
    void setTimePointer(const double* tp) { timePtr = tp; }

    Vec3 getCenter() const { return center; }
    Vec3 getAABBCenter() const { return aabbCenter; }  // AABB center (for rendering)
    Vec3 getHalfSize() const { return boxHalfSize; }   // AABB half-size (for rendering)
    Vec3 getDirection() const { return direction; }
    double getWidth() const { return width; }
    double getHeight() const { return height; }
    double getLength() const { return length; }
    double getStrength() const { return strength; }
    double getFrequency() const { return frequency; }

protected:
    // Snap direction to nearest cardinal axis (±X, ±Y, or ±Z)
    static Vec3 snapToAxis(const Vec3& d);
    
    // Update internal AABB representation and cached values
    void updateInternals();

protected:
    // User-facing parameters
    Vec3 center;            // center of the outlet face (where air comes out)
    double width, height;   // cross-section dimensions perpendicular to flow
    double length;          // how far the air reaches
    Vec3 direction;         // flow direction (axis-aligned: ±X, ±Y, or ±Z)
    double strength;        // acceleration magnitude
    double frequency;       // pulsation frequency in Hz (0 = steady)
    const double* timePtr;  // pointer to simulation time
    
    // Cached/computed values (for efficiency)
    Vec3 aabbCenter;        // AABB center (center + direction * length/2)
    Vec3 boxHalfSize;       // AABB half-extents
    Vec3 invHalfSize;       // precomputed 1/boxHalfSize (optimization)
};

#endif // FORCES_H
