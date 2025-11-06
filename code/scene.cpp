#include "scene.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void Scene::mousePressed(const QMouseEvent* e, const Camera& cam) {
    grabX = e->pos().x();
    grabY = e->pos().y();

    if (e->modifiers() & Qt::ControlModifier) {

        Vec3 rayDir = cam.getRayDir(grabX, grabY);
        Vec3 origin = cam.getPos();

        selectedCollider = -1;

        double closest = INFINITY;

        for (int i = 0; i < sceneColliders.size(); i++) {

            Collider* c = sceneColliders[i];
            double lambda;
            if (c->testRayIntersection(origin, rayDir, lambda)) {
                if (lambda < closest) {
                    closest = lambda;
                    selectedCollider = i;
                    
                }
            }
        }
    }
}
void Scene::mouseMoved(const QMouseEvent* e, const Camera& cam) {
    dx = e->pos().x() - grabX;
    dy = e->pos().y() - grabY;
    grabX = e->pos().x();
    grabY = e->pos().y();

    if (e->modifiers() & Qt::ControlModifier) {
        Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());
        if (selectedCollider >= 0) {
            sceneColliders[selectedCollider]->move(disp);
        }
    }
}
void Scene::mouseReleased(const QMouseEvent* e, const Camera& cam) {
    selectedCollider = -1;
}

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

QMatrix4x4 Scene::getPlaneModelMatrix(const ColliderPlane& plane, float scale) {
    Vec3 normal = plane.getNormal();
    double d = plane.getD();

    // Point on the plane: p = -d * normal (closest point to origin)
    Vec3 planePoint = -d * normal;

    QMatrix4x4 modelMat;
    modelMat.translate(planePoint[0], planePoint[1], planePoint[2]);

    // Rotate to align the quad with the plane normal
    // The quad is in XZ plane (normal = Y), we need to rotate it to match plane normal
    Vec3 up(0, 1, 0);  // quad's default normal

    // If normal is not aligned with Y axis, compute rotation
    if ((normal - up).norm() > 1e-6 && (normal + up).norm() > 1e-6) {
        Vec3 axis = up.cross(normal).normalized();
        double angle = std::acos(std::max(-1.0, std::min(1.0, up.dot(normal))));
        modelMat.rotate(angle * 180.0 / M_PI, axis[0], axis[1], axis[2]);
    }
    else if ((normal + up).norm() < 1e-6) {
        // normal points down, rotate 180 degrees
        modelMat.rotate(180.0, 1, 0, 0);
    }

    modelMat.scale(scale, 1, scale);
    return modelMat;
}

QMatrix4x4 Scene::getSphereModelMatrix(const ColliderSphere& sphere) {
    Vec3 center = sphere.getCenter();
    double radius = sphere.getRadius();

    QMatrix4x4 modelMat;
    modelMat.translate(center[0], center[1], center[2]);
    modelMat.scale(radius);
    return modelMat;
}

QMatrix4x4 Scene::getAABBModelMatrix(const ColliderAABB& box) {
    Vec3 center = box.getCenter();
    Vec3 halfSize = 0.5 * box.getSize();

    QMatrix4x4 modelMat;
    modelMat.translate(center[0], center[1], center[2]);
    modelMat.scale(halfSize[0], halfSize[1], halfSize[2]);
    return modelMat;
}