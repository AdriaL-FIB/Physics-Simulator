#include "colliders.h"
#include <cmath>


/*
 * Generic function for collision response from contact plane
 */
void Collider::resolveCollision(Particle* p, const Collision& col, double kElastic, double kFriction) const
{
    double d = -col.normal.dot(col.position);
    double h = col.normal.dot(p->pos) + d;

    Vec3 v_N = (col.normal.dot(p->vel)) * col.normal;
    Vec3 v_T = p->vel - v_N;

    Vec3 v_N1 = -kElastic * v_N;
    Vec3 v_T1 = (1 - kFriction) * v_T;

    Vec3 new_vel = v_N1 + v_T1;

    Vec3 new_pos = p->pos - ((1 + kElastic) * h) * col.normal;

    p->pos = new_pos;
    p->vel = new_vel;
}



/*
 * Plane
 */
bool ColliderPlane::isInside(const Particle* p) const
{
    return this->planeN.dot(p->pos) + this->planeD < 0 + p->radius;
}


bool ColliderPlane::testCollision(const Particle* p, Collision& colInfo) const
{
    Vec3 ray = p->pos - p->prevPos;

    double lambda = - (this->planeN.dot(p->prevPos) + this->planeD - p->radius) / this->planeN.dot(ray);

    colInfo.normal = this->planeN;
    colInfo.position = p->prevPos + lambda * ray;

    return 0 <= lambda and lambda <= 1;
}



/*
 * Sphere
 */
bool ColliderSphere::isInside(const Particle* p) const
{
    return (p->pos - this->center).squaredNorm() < (this->radius + p->radius)*(this->radius + p->radius);
}


bool ColliderSphere::testCollision(const Particle* p, Collision& colInfo) const
{
    Vec3 v = p->pos - p->prevPos;

    double a = v.dot(v);
    double b = 2 * v.dot(p->prevPos - this->center);
    double c = this->center.squaredNorm() + p->prevPos.squaredNorm() - 2*p->prevPos.dot(this->center) - (this->radius + p->radius) * (this->radius + p->radius);

    double D = b*b - 4.0*a*c;
    if (D < 0) return false;

    double sqrtD = std::sqrt(D);

    double l1 = (-b + sqrtD) / (2*a);
    double l2 = (-b - sqrtD) / (2*a);

    // double l2 = std::min(l1, l2);
    // double l1 = std::max(l1, l2);
    double l = (0.0 <= l2 && l2 <= 1.0) ? l2 : ((0.0 <= l1 && l1 <= 1.0) ? l1 : -1.0);

    if (l < 0) return false;

    Vec3 point = p->prevPos + l * v;
    colInfo.normal = (point - this->center).normalized();
    colInfo.position = point;

    return true;
}



/*
 * AABB
 */
bool ColliderAABB::isInside(const Particle* p) const
{
    return ( (p->pos.array() > bmin.array()) and (p->pos.array() < bmax.array()) ).all();
}

bool ColliderAABB::isInside(const Vec3& pos) const
{
    return ((pos.array() > bmin.array()) and (pos.array() < bmax.array())).all();

}


bool ColliderAABB::testCollision(const Particle* p, Collision& colInfo) const
{
    if (hollow and isInside(p->prevPos)) {
        Vec3 shrunkenBmin = this->bmin + Vec3(p->radius, p->radius, p->radius);
        Vec3 shrunkenBmax = this->bmax - Vec3(p->radius, p->radius, p->radius);

        Vec3 v = p->pos - p->prevPos;

        Vec3 lmin = (shrunkenBmin - p->prevPos).cwiseQuotient(v);
        Vec3 lmax = (shrunkenBmax - p->prevPos).cwiseQuotient(v);

        Vec3 l_near = lmin.cwiseMin(lmax);
        Vec3 l_far = lmin.cwiseMax(lmax);

        Eigen::Index hitAxis;

        double l_in = l_near.maxCoeff();
        double l_out = l_far.minCoeff(&hitAxis);
        bool collision = (l_in <= l_out) && 0 <= l_out && l_out <= 1;

        if (not collision) return false;

        colInfo.position = p->prevPos + l_out * v;
        colInfo.normal = Vec3::Zero();
        colInfo.normal[hitAxis] = v[hitAxis] > 0 ? 1 : -1;

        return true;
    }
    else {
        Vec3 expandedBmin = this->bmin - Vec3(p->radius, p->radius, p->radius);
        Vec3 expandedBmax = this->bmax + Vec3(p->radius, p->radius, p->radius);

        Vec3 v = p->pos - p->prevPos;

        Vec3 lmin = (expandedBmin - p->prevPos).cwiseQuotient(v);
        Vec3 lmax = (expandedBmax - p->prevPos).cwiseQuotient(v);

        Vec3 l_near = lmin.cwiseMin(lmax);
        Vec3 l_far  = lmin.cwiseMax(lmax);

        Eigen::Index hitAxis;

        double l_in = l_near.maxCoeff(&hitAxis);
        bool collision = (l_in <= l_far.minCoeff()) && 0 <= l_in && l_in <= 1;

        if (not collision) return false;

        colInfo.position = p->prevPos + l_in * v;
        colInfo.normal = Vec3::Zero();
        colInfo.normal[hitAxis] = v[hitAxis] > 0 ? -1 : 1;

        return true;
    }
}