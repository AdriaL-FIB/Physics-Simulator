#include "forces.h"

void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        // TODO
    }
}

void ForceDrag::apply() {
    for (Particle* p : particles) {
        // TODO
    }
}

void ForceSpring::apply() {
    if (particles.size() < 2) return;
    Particle* p1 = getParticle1();
    Particle* p2 = getParticle2();

    // TODO
}
