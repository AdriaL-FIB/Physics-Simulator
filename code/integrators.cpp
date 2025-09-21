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


void IntegratorImplicitEuler::step(ParticleSystem& system, double dt) {
    double t0 = system.getTime();
    Vecd X0 = system.getState();

    // Newton-Raphson method

    Vecd dx = system.getDerivative();
    Vecd X_guess = X0 + dt * dx;

    for (int i = 0; i < this->MAX_NEWTON_ITERATIONS; ++i) {
        Vecd b = -evaluateResiudal(system, X_guess, dt);
        
        Matd A = computeJacobi(system, X_guess, dt);

        Vecd deltaX = A.colPivHouseholderQr().solve(b);

        X_guess += deltaX;

        if (deltaX.norm() < this->CONVERGENCE_THRESHOLD) {
            break;
        }
    }

    system.setState(X_guess);
    system.setTime(t0 + dt);
    system.updateForces();

}

Vecd IntegratorImplicitEuler::evaluateResiudal(ParticleSystem& system, const Vecd& XGuess, double dt) {
    double t0 = system.getTime();
    Vecd X0 = system.getState();

    system.setState(XGuess);
    system.setTime(t0 + dt);
    system.updateForces();

    Vecd dXGuess = system.getDerivative(); // (V_guess, F_guess/m)

    Vecd G = XGuess - X0 - dt * dXGuess;

    system.setState(X0);
    system.setTime(t0);
    system.updateForces();

    return G;
}

Matd IntegratorImplicitEuler::computeJacobi(ParticleSystem& system, const Vecd& XGuess, double dt) {
    double t0 = system.getTime();
    Vecd X0 = system.getState();

    int N = system.getStateSize();
    const double epsilon = 1e-6;

    system.setState(XGuess);
    system.setTime(t0 + dt);
    system.updateForces();

    Vecd dX_base = system.getDerivative();

    Matd J_f(N, N);

    // Change 1 component at a time
    for (int j = 0; j < N; ++j) {
        Vecd epsilonState = XGuess;
        epsilonState[j] += epsilon;
        
        system.setState(epsilonState);
        system.updateForces();

        Vecd dX_epsilon = system.getDerivative();

        J_f.col(j) = (dX_epsilon - dX_base) / epsilon;
    }

    Matd A = Matd::Identity(N, N) - dt * J_f;

    // Restore state
    system.setState(X0);
    system.setTime(t0);
    system.updateForces();

    return A;
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
