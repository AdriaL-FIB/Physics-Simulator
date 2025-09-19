#include "integrators.h"


void IntegratorEuler::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd dx = system.getDerivative();
    Vecd x1 = x0 + dt*dx;
    system.setState(x1);
    system.setTime(t0+dt);
    system.updateForces();
}


void IntegratorSymplecticEuler::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd dx = system.getDerivative();

    Vecd a0 = system.getAccelerations();

    Vecd x1 = x0 + dt*dx;

    system.setState(x1);
    Vecd v1 = system.getVelocities();

    system.setState(x0);
    system.setVelocities(v1);

    dx = system.getDerivative(); // v1 , a0

    Vecd x = x0 + dt * dx;

    system.setState(x);
    system.setTime(t0+dt);
    system.updateForces();
}


void IntegratorMidpoint::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd dx = system.getDerivative();
    Vecd x1 = x0 + (dt/2)*dx; // updates position and velocity
    system.setState(x1);

    // Required?
    system.setTime(t0+dt/2);
    system.updateForces();
    // ---------
    
    Vecd dx_mid = system.getDerivative();
    Vecd x = x0 + dt * dx_mid;
    system.setState(x);
    system.setTime(t0+dt);
    system.updateForces();
}


void IntegratorRK2::step(ParticleSystem &system, double dt) {
    double t1 = system.getTime();
    Vecd x1 = system.getState();
    Vecd k1 = system.getDerivative();
    Vecd x2 = x1 + dt * k1;
    system.setState(x2);

    // Required?
    system.setTime(t1+dt);
    system.updateForces();
    // ---------

    Vecd k2 = system.getDerivative();
    Vecd x = x1 + dt/2 * (k1 + k2);
    system.setState(x);
    system.setTime(t1+dt);
    system.updateForces();
}


void IntegratorRK4::step(ParticleSystem &system, double dt) {
    double t0 = system.getTime();
    Vecd x0 = system.getState();
    Vecd k1 = system.getDerivative();
    
    Vecd xk2 = x0 + dt/2 * k1;

    system.setState(xk2);
    system.setTime(t0 + dt/2);
    system.updateForces();

    Vecd k2 = system.getDerivative();

    Vecd xk3 = x0 + dt/2 * k2;

    system.setState(xk3);
    system.setTime(t0 + dt/2);
    system.updateForces();

    Vecd k3 = system.getDerivative();

    Vecd xk4 = x0 + dt * k3;

    system.setState(xk4);
    system.setTime(t0 + dt);
    system.updateForces();

    Vecd k4 = system.getDerivative();

    Vecd x = x0 + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4);

    system.setState(x);
    system.setTime(t0 + dt);
    system.updateForces();
}


void IntegratorVerlet::step(ParticleSystem &system, double dt) {
    // TODO
}
