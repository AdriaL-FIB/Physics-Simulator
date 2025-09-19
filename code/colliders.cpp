#include "colliders.h"
#include <cmath>


/*
 * Generic function for collision response from contact plane
 */
void Collider::resolveCollision(Particle* p, const Collision& col, double kElastic, double kFriction) const
{
    // TODO
}



/*
 * Plane
 */
bool ColliderPlane::isInside(const Particle* p) const
{
    // TODO
    return false;
}


bool ColliderPlane::testCollision(const Particle* p, Collision& colInfo) const
{
    // TODO
    return false;
}



/*
 * Sphere
 */
bool ColliderSphere::isInside(const Particle* p) const
{
    // TODO
    return false;
}


bool ColliderSphere::testCollision(const Particle* p, Collision& colInfo) const
{
    // TODO
    return false;
}



/*
 * AABB
 */
bool ColliderAABB::isInside(const Particle* p) const
{
    // TODO
    return false;
}


bool ColliderAABB::testCollision(const Particle* p, Collision& colInfo) const
{
    // TODO
    return false;
}
