#ifndef SCENE_H
#define SCENE_H

#include <QWidget>
#include <QMouseEvent>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include "camera.h"
#include "integrators.h"
#include "colliders.h"

class Scene : public QObject
{
    Q_OBJECT
public:
    Scene() {}
    virtual ~Scene() {};

    virtual void initialize() = 0;
    virtual void reset() = 0;
    virtual void update(double dt) = 0;
    virtual void paint(const Camera& cam) = 0;

    virtual void mousePressed (const QMouseEvent*, const Camera&);
    virtual void mouseMoved   (const QMouseEvent*, const Camera&);
    virtual void mouseReleased(const QMouseEvent*, const Camera&);
    virtual void keyPressed   (const QKeyEvent*,   const Camera&) {};

    virtual void getSceneBounds(Vec3& bmin, Vec3& bmax) = 0;
    virtual unsigned int getNumParticles() { return 0; }

    virtual QWidget* sceneUI() = 0;

    void clearBuffers() {
        for (QOpenGLBuffer* buf : buffers) {
            buf->destroy();
            delete buf;
        }
        buffers.clear();
    }
public:
    Integrator* createIntegrator(int type);
    inline const std::vector<std::string> getIntegratorNames() {
        return {
            "Euler", "Symplectic Euler", "Implicit Euler (Taylor)", "Implicit Euler (Newton)",
            "Midpoint", "RK2", "RK4", "Verlet"
        };
    }

    // Helper functions to get model matrices for rendering colliders
    static QMatrix4x4 getPlaneModelMatrix(const ColliderPlane& plane, float scale = 100.0f);
    static QMatrix4x4 getSphereModelMatrix(const ColliderSphere& sphere);
    static QMatrix4x4 getAABBModelMatrix(const ColliderAABB& box);

protected:
    std::vector<QOpenGLBuffer*> buffers;

    std::vector<Collider*> sceneColliders;
    int selectedCollider = -1;

    // mouse interaction
    int grabX, grabY, dx, dy;
};

#endif // SCENE_H
