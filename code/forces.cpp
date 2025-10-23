#include "forces.h"
#include "math.h"
#include "particlesystem.h"

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
        //p->force += -this->kquadratic * p->vel.cwiseAbs().cwiseProduct(p->vel);
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
            Particle* pj = particles[nbs[j]];

            Vec3 r = pj->pos - p->pos;

            double rho1 = densities[i];
            double rho2 = densities[j];
            
            double Pij = -pj->mass*(pressures[i]/(rho1*rho1) + pressures[j]/(rho2*rho2));

            aPressure += Pij * kernelGradientSpiky(r);

            Vec3 Vij = mu * pj->mass * (pj->vel - p->vel)/(rho1*rho2);

            aViscosity += Vij * kernelLaplacianViscosity(r.norm());
        }

        p->force += aPressure + aViscosity;
    }
}

double ForceSPH::kernelPoly6(double rSqr) // TODO: If I already have |r|, use r instead of sqr.
{
    if (rSqr > h * h) return 0;

    double term = h*h - rSqr;

    return poly6Factor * term * term * term;
}

Vec3 ForceSPH::kernelGradientSpiky(Vec3 r)
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
    
    return viscosiyLaplacianFactor * (1 - r/h);
}

double ForceSPH::computePressure(double density)
{
    return 100*(density - 1);
}
