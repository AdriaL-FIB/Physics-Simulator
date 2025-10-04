#include "scene.h"

Integrator *Scene::createIntegrator(int type) {
    switch (type) {
        case 0: return new IntegratorEuler();
        case 1: return new IntegratorSymplecticEuler();
        case 2: return new IntegratorImplicitEuler();
        case 3: return new IntegratorImplicitEuler(1);
        case 4: return new IntegratorMidpoint();
        //case 5:
        case 6: return new IntegratorRK2();
        case 7: return new IntegratorRK4();
        case 8: return new IntegratorVerlet();
        default: return new IntegratorEmpty();
    }
}