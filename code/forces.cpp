#include "forces.h"
#include "math.h"
#include "particlesystem.h"
#include <iostream>
#include <cmath>

void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        p->force += p->mass * this->acceleration;
    }
}

void ForceDrag::apply() {
    for (Particle* p : particles) {
        if (useQuadratic) {
            p->force += -this->kquadratic * p->vel.norm() * (p->vel);
        } else {
            p->force += -this->klinear * p->vel;
        }
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

    // Cache neighbor information to avoid recomputing
    std::vector<std::vector<int>> neighborCache(N);

    // Compute density
    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        double density = p->mass * kernelPoly6(0);

        std::vector<int> nbs = system->getNeighbors(p);
        neighborCache[i] = nbs;
        
        for (int j = 0; j < nbs.size(); ++j) {
            Particle* pj = particles[nbs[j]];
            density += pj->mass * kernelPoly6((pj->pos - p->pos).squaredNorm());
        }
        densities[i] = density;
    }

    // Compute pressure
    for(int i = 0; i < N; ++i) {
        pressures[i] = computePressure(densities[i]);
    }

    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];

        Vec3 aPressure(0,0,0);
        Vec3 aViscosity(0,0,0);
        Vec3 aRepulsion(0,0,0);  // Emergency repulsion for collapsed particles


        std::vector<int> nbs = neighborCache[i];
        
        for (int j = 0; j < nbs.size(); ++j) {
            int neighborIndx = nbs[j];
            Particle* pj = particles[neighborIndx];

            Vec3 r = pj->pos - p->pos;
            double dist = r.norm();

            double rho1 = densities[i];
            double rho2 = densities[neighborIndx];
            
            double Pij = -pj->mass*(pressures[i]/(rho1*rho1) + pressures[neighborIndx]/(rho2*rho2));

            aPressure += Pij * kernelGradientSpiky(r);

            Vec3 Vij = mu * pj->mass * (pj->vel - p->vel)/(rho1*rho2);

            aViscosity += Vij * kernelLaplacianViscosity(dist);
            
            double minDist = 1.2 * p->radius;
            if (dist < minDist) {
                double overlap = (minDist - dist) / minDist; // 0 to 1
                double repulsionMagnitude = 500.0 * overlap * overlap;
                Vec3 repulsionDir = r / dist;
                aRepulsion -= repulsionDir * repulsionMagnitude; 
            }
        }

        p->force += p->mass * (aPressure + aViscosity + aRepulsion);
    }
    
    if (colorMode == SPHColorMode::None) {
        if (not hasInitializedDefaultColors) {
            for(int i = 0; i < N; ++i) {
                particles[i]->color = particles[i]->defaultColor;
            }
            hasInitializedDefaultColors = true;
        }
    } else if (colorMode == SPHColorMode::Neighbors) {
        for(int i = 0; i < N; ++i) {
            int numNeighbors = neighborCache[i].size();
            double t = std::min(1.0, numNeighbors / 50.0);
            particles[i]->color = computeColorGradient(t);
        }
    } else if (colorMode == SPHColorMode::Density) {
        for(int i = 0; i < N; ++i) {
            double densityRatio = densities[i] / p0;
            double t = std::clamp((densityRatio - 0.2) / 0.5, 0.0, 1.0);
            particles[i]->color = computeColorGradient(t);
        }
    } else if (colorMode == SPHColorMode::Pressure) {
        for(int i = 0; i < N; ++i) {
            double t = (pressures[i] + 10 / 10);
            particles[i]->color = computeColorGradient(t);
        }
    }
}

double ForceSPH::kernelPoly6(double rSqr) const
{
    if (rSqr > h * h) return 0;

    double term = h*h - rSqr;

    return poly6Factor * term * term * term;
}

Vec3 ForceSPH::kernelGradientSpiky(Vec3 r) const
{
    double rMag = r.norm();

    if (rMag > h) {
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
    int flowAxis = 0;
    if (std::abs(direction[1]) > 0.5) flowAxis = 1;
    else if (std::abs(direction[2]) > 0.5) flowAxis = 2;
    
    aabbCenter = center + direction * (length * 0.5);
    
    // Half-size: length/2 along flow axis, width/2 and height/2 on perpendicular axes
    boxHalfSize[flowAxis] = length * 0.5;
    

    switch (flowAxis) {
        case 0:
            boxHalfSize[1] = height * 0.5; // Height in Y
            boxHalfSize[2] = width * 0.5;
            break;
        case 1:
            boxHalfSize[0] = width * 0.5;
            boxHalfSize[2] = height * 0.5; // Height in Z
            break;
        case 2:
            boxHalfSize[0] = width * 0.5;
            boxHalfSize[1] = height * 0.5; // Height in Y
            break;
    }

    for (int i = 0; i < 3; i++) {
        invHalfSize[i] = (boxHalfSize[i] != 0.0) ? 1.0 / boxHalfSize[i] : 0.0;
    }
}

void ForceWind::apply()
{
    for (Particle* p : particles) {
        Vec3 d = p->pos - aabbCenter;

        // Quick AABB test. if outside on any axis, skip
        if (std::abs(d[0]) > boxHalfSize[0]) continue;
        if (std::abs(d[1]) > boxHalfSize[1]) continue;
        if (std::abs(d[2]) > boxHalfSize[2]) continue;

        // Compute normalized squared distance within box coordinates (no sqrt)
        double nx = d[0] * invHalfSize[0];
        double ny = d[1] * invHalfSize[1];
        double nz = d[2] * invHalfSize[2];
        double nsq = nx*nx + ny*ny + nz*nz; // [0,3]

        double s = nsq / 3.0; // s in [0,1]

        double falloff = 1.0 - s;

        // Smooth the falloff with a cubic smoothstep: f = f*f*(3-2*f)
        double f2 = falloff * falloff;
        falloff = f2 * (3.0 - 2.0 * falloff);

        double accelMag = strength * falloff;

        p->force += direction * (accelMag * p->mass);
    }
}

// Blue (0.0) -> Cyan -> Green -> Yellow -> Red (1.0)
Vec3 ForceSPH::computeColorGradient(double t) const {
    // Clamp t to [0, 1]
    t = std::clamp(t, 0.0, 1.0);
    
    if (t < 0.25) {
        // Blue to Cyan
        double localT = t / 0.25;
        return Vec3(0.0, localT, 1.0);
    } else if (t < 0.5) {
        // Cyan to Green
        double localT = (t - 0.25) / 0.25;
        return Vec3(0.0, 1.0, 1.0 - localT);
    } else if (t < 0.75) {
        // Green to Yellow
        double localT = (t - 0.5) / 0.25;
        return Vec3(localT, 1.0, 0.0);
    } else {
        // Yellow to Red
        double localT = (t - 0.75) / 0.25;
        return Vec3(1.0, 1.0 - localT, 0.0);
    }
}
