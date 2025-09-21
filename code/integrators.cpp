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

	Vecd v1 = system.getVelocities() + dt * system.getAccelerations(); // updates velocities
	Vecd p1 = system.getPositions() + dt * v1; // updates positions


	system.setPositions(p1);
	system.setVelocities(v1);
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
    double t0 = system.getTime();
	Vecd p0 = system.getPositions();
    Vecd prev = system.getPreviousPositions();
	Vecd a0 = system.getAccelerations();
    Vecd v0 = system.getVelocities();

    if (t0 < dt) {
        prev = p0 - v0 * dt;
    }


	Vecd p1 = p0 + this->kd * (p0 - prev) + a0 * dt * dt; // updates positions

	Vecd v1 = (p1 - p0) / (dt); // updates velocities

	system.setPositions(p1);
	system.setVelocities(v1);
    system.setPreviousPositions(p0);
	system.setTime(t0 + dt);
	system.updateForces();
}
