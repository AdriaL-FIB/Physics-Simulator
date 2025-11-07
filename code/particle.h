#ifndef PARTICLE_H
#define PARTICLE_H

#include "defines.h"

class Particle
{
public:

    static const int PhaseDimension = 6;

    Vec3 pos, prevPos;
    Vec3 vel;
    Vec3 force;
    double mass;
    double invMass;
    double radius = 1.0;
    double life   = 0.0;
    Vec3 color    = Vec3(153 / 255.0, 217 / 255.0, 234 / 255.0);
    Vec3 defaultColor = Vec3(153 / 255.0, 217 / 255.0, 234 / 255.0);
    unsigned int id = 0;

    Particle() {
        pos	    = Vec3(0.0, 0.0, 0.0);
        vel	    = Vec3(0.0, 0.0, 0.0);
        force   = Vec3(0.0, 0.0, 0.0);
        prevPos = pos;
        mass    = 1.0;
        invMass = 1.0;
    }

    Particle(const Vec3& p, const Vec3& v, float m) {
        pos		= p;
        vel		= v;
        force	= Vec3(0.0, 0.0, 0.0);
        prevPos = pos;
        mass	= m;
        invMass = 1 / m;
    }

    Particle(const Particle& p) {
        id      = p.id;
        pos     = p.pos;
        vel     = p.vel;
        force   = p.force;
        prevPos = p.pos;
        mass    = p.mass;
        invMass = p.invMass;
        color   = p.color;
        radius  = p.radius;
        life    = p.life;
    }

    void setMass(double m) { mass = m; invMass = 1 / m; }
    void fix(bool b) {
        if (b) invMass = 0;
        else invMass = 1 / mass;
    }

    ~Particle() {
    }
};


#endif // PARTICLE_H
