#include "forces.h"
#include "math.h"
#include "particlesystem.h"
#include <iostream>

void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        p->force += p->mass * this->acceleration;
    }
}

void ForceDrag::apply() {
    for (Particle* p : particles) {
        // Strokes drag:
        p->force += -this->klinear * p->vel;

        // Newton drag:
        //p->force += -this->kquadratic * p->vel.cwiseAbs() * (p->vel);
    }
}

void ForceSpring::apply() {
    if (particles.size() < 2) return;
    Particle* p1 = getParticle1();
    Particle* p2 = getParticle2();

    double norm = (p2->pos - p1->pos).norm();
    Vec3 udir = (p2->pos - p1->pos)/norm;

    Vec3 F = (this->ks * (norm - this->L) + this->kd * (p2->vel - p1->vel).dot(udir)) * (udir);

    p1->force += F;
    p2->force -= F;
}

void ForceGravitation::apply() {
    const Particle* p1 = attractor;
    for (unsigned int j = 0; j < particles.size(); ++j) {
        Particle* p2 = particles[j];

        Vec3 diff = p2->pos - p1->pos;
        double distSquared = diff.squaredNorm();

        Vec3 F = (G * p1->mass * p2->mass) / distSquared * diff / diff.norm();

        F *= 2 / (1 + exp(-a * distSquared / (b * b))) - 1;

        p2->force -= F;
    }
}

void ForceSPH::apply()
{
    system->buildSpatialHashing(h);

    int N = particles.size();

    densities.resize(N);
    pressures.resize(N);

    // Compute density
    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        double density = p->mass * kernelPoly6(0);

        std::vector<int> nbs = system->getNeighbors(p);
        for (int j = 0; j < nbs.size(); ++j) {
            Particle* pj = particles[nbs[j]];
            density += pj->mass * kernelPoly6((pj->pos - p->pos).squaredNorm());
        }
        densities[i] = density;
    }

    // Compute pressure
    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        pressures[i] = computePressure(densities[i]);
    }

    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];

        Vec3 aPressure(0,0,0);
        Vec3 aViscosity(0,0,0);

        // TODO: optimize
        std::vector<int> nbs = system->getNeighbors(p);
        for (int j = 0; j < nbs.size(); ++j) {
            int neighborIndx = nbs[j];
            Particle* pj = particles[neighborIndx];

            Vec3 r = pj->pos - p->pos;

            double rho1 = densities[i];
            double rho2 = densities[neighborIndx];
            
            double Pij = -pj->mass*(pressures[i]/(rho1*rho1) + pressures[neighborIndx]/(rho2*rho2));

            aPressure += Pij * kernelGradientSpiky(r);

            Vec3 Vij = mu * pj->mass * (pj->vel - p->vel)/(rho1*rho2);

            aViscosity += Vij * kernelLaplacianViscosity(r.norm());
        }

        p->force += p->mass * (aPressure + aViscosity);

        // double pressureNorm = std::clamp(pressures[i] / 100.0, 0.0, 100.0);
        // pressureNorm = (pressureNorm + 1.0) * 0.5; // map from [-1,1] to [0,1]
        
        // // Soft colors: light blue (0.5, 0.7, 1.0) -> light red (1.0, 0.5, 0.5)
        // p->color = Vec3(
        //     0.5 + 0.5 * pressureNorm,      // R: 0.5 -> 1.0
        //     0.7 - 0.2 * pressureNorm,      // G: 0.7 -> 0.5
        //     1.0 - 0.5 * pressureNorm       // B: 1.0 -> 0.5
        // );
    }

    // Print densities and pressure stats (min, max, avg, quaredAverage)
    // if (N > 0) {
    //     double minDensity = densities[0];
    //     double maxDensity = densities[0];
    //     double avgDensity = 0.0;
    //     double sqAvgDensity = 0.0;
        
    //     double minPressure = pressures[0];
    //     double maxPressure = pressures[0];
    //     double avgPressure = 0.0;
    //     double sqAvgPressure = 0.0;
        
    //     for (int i = 0; i < N; ++i) {
    //         // Density stats
    //         minDensity = std::min(minDensity, densities[i]);
    //         maxDensity = std::max(maxDensity, densities[i]);
    //         avgDensity += densities[i];
    //         sqAvgDensity += densities[i] * densities[i];
            
    //         // Pressure stats
    //         minPressure = std::min(minPressure, pressures[i]);
    //         maxPressure = std::max(maxPressure, pressures[i]);
    //         avgPressure += pressures[i];
    //         sqAvgPressure += pressures[i] * pressures[i];
    //     }
        
    //     avgDensity /= N;
    //     sqAvgDensity = std::sqrt(sqAvgDensity / N);
    //     avgPressure /= N;
    //     sqAvgPressure = std::sqrt(sqAvgPressure / N);
        
    //     std::cout << "Density  - Min: " << minDensity << ", Max: " << maxDensity 
    //               << ", Avg: " << avgDensity << ", RMS: " << sqAvgDensity << std::endl;
    //     std::cout << "Pressure - Min: " << minPressure << ", Max: " << maxPressure 
    //               << ", Avg: " << avgPressure << ", RMS: " << sqAvgPressure << std::endl;
    // }
}

double ForceSPH::kernelPoly6(double rSqr) const // TODO: If I already have |r|, use r instead of sqr.
{
    if (rSqr > h * h) return 0;

    double term = h*h - rSqr;

    return poly6Factor * term * term * term;
}

Vec3 ForceSPH::kernelGradientSpiky(Vec3 r) const
{
    double rMag = r.norm();

    if (rMag <= 0.0 || rMag > h) {
        return Vec3(0.0, 0.0, 0.0);
    }

    return -r * spikyGradFactor / rMag * (h - rMag) * (h - rMag);
}

double ForceSPH::kernelLaplacianViscosity(double r) const
{
    if (r <= 0.0 || r > h) {
        return 0.0;
    }
    
    return viscosityLaplacianFactor * (1 - r/h);
}

double ForceSPH::computePressure(double density)
{
    const int k = 50;
    return k * (density - p0);
}

ForceWind::ForceWind() {
    center = Vec3(0, 0, 0);
    width = height = 10.0;
    length = 10.0;
    direction = Vec3(1, 0, 0);
    strength = 10.0;
    frequency = 0.0;
    timePtr = nullptr;
    updateInternals();
}

ForceWind::ForceWind(const Vec3& center, double width, double height, const Vec3& dir, double length, double strength)
    : center(center), width(width), height(height), direction(snapToAxis(dir)),
      length(length), strength(strength), frequency(0.0), timePtr(nullptr) {
    updateInternals();
}

void ForceWind::setDirection(const Vec3& d) {
    direction = snapToAxis(d);
    updateInternals();
}

Vec3 ForceWind::snapToAxis(const Vec3& d) {
    Vec3 dn = d.normalized();
    int maxIdx = 0;
    double maxAbs = std::abs(dn[0]);
    for (int i = 1; i < 3; i++) {
        if (std::abs(dn[i]) > maxAbs) {
            maxAbs = std::abs(dn[i]);
            maxIdx = i;
        }
    }
    Vec3 snapped = Vec3(0, 0, 0);
    snapped[maxIdx] = (dn[maxIdx] > 0) ? 1.0 : -1.0;
    return snapped;
}

void ForceWind::updateInternals() {
    // Since direction is axis-aligned (from snapToAxis), find which axis it is
    // by finding the component with absolute value 1
    int flowAxis = 0;
    if (std::abs(direction[1]) > 0.5) flowAxis = 1;
    else if (std::abs(direction[2]) > 0.5) flowAxis = 2;
    
    // Build AABB: aabbCenter is offset from outlet center by length/2 in flow direction
    aabbCenter = center + direction * (length * 0.5);
    
    // Half-size: length/2 along flow axis, width/2 and height/2 on perpendicular axes
    boxHalfSize = Vec3(0, 0, 0);
    boxHalfSize[flowAxis] = length * 0.5;
    
    // Assign width and height to the two perpendicular axes
    int axis1 = (flowAxis + 1) % 3;
    int axis2 = (flowAxis + 2) % 3;
    boxHalfSize[axis1] = width * 0.5;
    boxHalfSize[axis2] = height * 0.5;
    
    // Precompute inverse half-sizes (optimization: only once when geometry changes)
    invHalfSize = Vec3(0, 0, 0);
    for (int i = 0; i < 3; i++) {
        invHalfSize[i] = (boxHalfSize[i] != 0.0) ? 1.0 / boxHalfSize[i] : 0.0;
    }
}

void ForceWind::apply()
{
    // Get current time from the pointer
    double time = (timePtr != nullptr) ? *timePtr : 0.0;
    
    // Cheap pulsation factor (0 = no pulsation). Keep sine (acceptable cost) for natural pulsation.
    double pulseFactor = 1.0;
    if (frequency > 0.0) {
        pulseFactor = 0.5 * (1.0 + std::sin(2.0 * M_PI * frequency * time));
    }

    // invHalfSize is now precomputed in updateInternals() - no need to compute it here!

    // For each influenced particle, check AABB membership and apply a cheap polynomial falloff
    for (Particle* p : particles) {
        // Skip dead particles (important for particle recycling systems)
        if (p->life <= 0.0) continue;
        
        // Compute local vector from AABB center to particle
        Vec3 d = p->pos - aabbCenter;

        // Quick AABB test (cheap): if outside on any axis, skip
        if (std::abs(d[0]) > boxHalfSize[0]) continue;
        if (std::abs(d[1]) > boxHalfSize[1]) continue;
        if (std::abs(d[2]) > boxHalfSize[2]) continue;

        // Compute normalized squared distance within box coordinates (no sqrt)
        double nx = d[0] * invHalfSize[0];
        double ny = d[1] * invHalfSize[1];
        double nz = d[2] * invHalfSize[2];
        double nsq = nx*nx + ny*ny + nz*nz; // in [0,3] for corners

        // If outside a unit sphere in normalized coords, we still keep them because we already
        // passed the AABB check. Use nsq to compute a smooth falloff that is cheap.
        // Map nsq to falloff in [0,1] where center -> 1, outer -> 0.
        // Normalize nsq by 3 to map corner to 1.0
        double s = nsq / 3.0; // s in [0,1]
        if (s >= 1.0) continue; // effectively outside influence (cheap)

        double falloff = 1.0 - s; // linear
        // Smooth the falloff with a cubic smoothstep: f = f*f*(3-2*f)
        double f2 = falloff * falloff;
        falloff = f2 * (3.0 - 2.0 * falloff);

        // Final acceleration magnitude (treat strength as acceleration magnitude)
        double accelMag = strength * falloff * pulseFactor;

        // Apply as F = m * a (consistent with other forces)
        p->force += direction * (accelMag * p->mass);
    }
}

