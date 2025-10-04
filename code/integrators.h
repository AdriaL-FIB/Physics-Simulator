#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "particlesystem.h"

#include <iostream>

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
    IntegratorImplicitEuler(int method = 0) { this->method = method; }

    virtual void step(ParticleSystem& system, double dt);
    const int MAX_NEWTON_ITERATIONS = 100;
    const double CONVERGENCE_THRESHOLD = 1e-3;

	int method; // 0: Taylor, 1: Newton-Raphson

    Vecd evaluateResiudal(ParticleSystem& system, const Vecd& XGuess, double dt);
    Matd computeJacobi(ParticleSystem& system, const Vecd& XGuess, double dt);

    Matd computeStiffnessMatrix(ParticleSystem& system);
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

class IntegratorEmpty : public Integrator {
public:
    IntegratorEmpty() {
        std::cout << "Empty integrator created" << std::endl;
    }
    virtual void step(ParticleSystem& system, double dt) {};
};


#endif // INTEGRATORS_H
