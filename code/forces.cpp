#include "forces.h"
#include "math.h"

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

double ForceSPH::kernelLaplacianViscosity(double r)
{
    if (r <= 0.0 || r > h) {
        return 0.0;
    }
    
    return viscosiyLaplacianFactor * (1 - r/h);
}
