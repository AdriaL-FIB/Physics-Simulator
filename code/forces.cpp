#include "forces.h"

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
        p->force += -this->kquadratic * p->vel.cwiseAbs().cwiseProduct(p->vel);
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
    for (unsigned int i = 0; i < particles.size(); ++i) {
        for (unsigned int j = i + 1; j < particles.size(); ++j) {
            Particle* p1 = particles[i];
            Particle* p2 = particles[j];

            Vec3 diff = p2->pos - p1->pos;

            Vec3 F = (G * p1->mass * p2->mass) / diff.squaredNorm() * diff / diff.norm();

            p1->force += F;
            p2->force -= F;
        }
    }
}
