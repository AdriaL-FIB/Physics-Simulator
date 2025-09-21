#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "particlesystem.h"

class Integrator {
public:
    Integrator() {};
    virtual ~Integrator() {};
    virtual void step(ParticleSystem& system, double dt) = 0;
};


class IntegratorEuler : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
};


class IntegratorSymplecticEuler : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
};

class IntegratorImplicitEuler : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
    const int MAX_NEWTON_ITERATIONS = 100;
    const double CONVERGENCE_THRESHOLD = 1e-3;

    Vecd evaluateResiudal(ParticleSystem& system, const Vecd& XGuess, double dt);
    Matd computeJacobi(ParticleSystem& system, const Vecd& XGuess, double dt);
};


class IntegratorMidpoint : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
};

class IntegratorRK2 : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
};

class IntegratorRK4 : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
};

class IntegratorVerlet : public Integrator {
public:
    virtual void step(ParticleSystem& system, double dt);
    double kd = 1;
};


#endif // INTEGRATORS_H
