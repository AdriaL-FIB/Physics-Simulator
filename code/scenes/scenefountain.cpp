#include "scenefountain.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>


SceneFountain::SceneFountain() {
    widget = new WidgetFountain();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
    widget->setSolverTypes(getIntegratorNames());
}


SceneFountain::~SceneFountain() {
    if (widget)     delete widget;
    if (shader)     delete shader;
    if (vaoFloor)   delete vaoFloor;
    if (vaoSphereH) delete vaoSphereH;
    if (vaoSphereL) delete vaoSphereL;
    if (vaoCube)    delete vaoCube;
    if (fGravity)   delete fGravity;
}


void SceneFountain::initialize() {
    // load shader
    shader = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");

    // create floor VAO
    Model quad = Model::createQuad();
    vaoFloor = glutils::createVAO(shader, &quad, buffers);
    glutils::checkGLError();

    // create particle VAOs
    Model sphereLowres = Model::createIcosphere(1);
    vaoSphereL = glutils::createVAO(shader, &sphereLowres, buffers);
    numFacesSphereL = sphereLowres.numFaces();
    glutils::checkGLError();

    // create sphere VAO
    Model sphere = Model::createIcosphere(3);
    vaoSphereH = glutils::createVAO(shader, &sphere, buffers);
    numFacesSphereH = sphere.numFaces();
    glutils::checkGLError();

    // create box VAO
    Model cube = Model::createCube();
    vaoCube = glutils::createVAO(shader, &cube, buffers);
    glutils::checkGLError();

    // create forces
    fGravity = new ForceConstAcceleration();
    system.addForce(fGravity);

    // scene description
    fountainPos = Vec3(0, 80, 0);    
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);    
    colliderRamp.setPlane(Vec3(0, std::sqrt(3.0)/2.0, 0.5), 6);
    colliderSphere.setCenter(Vec3(0,0,0));
    colliderSphere.setRadius(20);
    colliderBox.setFromBounds(Vec3(30,0,20), Vec3(50,10,60));

    colliderBall1.setCenter(Vec3(-40, 15, -30));
    colliderBall1.setRadius(12);
    colliderBall2.setCenter(Vec3(40, 25, 30));
    colliderBall2.setRadius(15);
    colliderBall3.setCenter(Vec3(-30, 35, 40));
    colliderBall3.setRadius(10);

    colliderPillar1.setFromBounds(Vec3(-60, 0, -60), Vec3(-50, 30, -40));
    colliderPillar2.setFromBounds(Vec3(50, 0, -20), Vec3(60, 25, 0));
    colliderPillar3.setFromBounds(Vec3(-10, 10, 50), Vec3(10, 35, 60));

    colliderPlatform1.setFromBounds(Vec3(-25, 40, -25), Vec3(-5, 42, -5));
    colliderPlatform2.setFromBounds(Vec3(10, 50, 10), Vec3(30, 52, 30));

    colliderWallLeft.setFromBounds(Vec3(-70, 0, -70), Vec3(-68, 40, 70));
    colliderWallRight.setFromBounds(Vec3(68, 0, -70), Vec3(70, 40, 70));

    sceneColliders.push_back(&colliderFloor);
    sceneColliders.push_back(&colliderRamp);
    sceneColliders.push_back(&colliderSphere);
    sceneColliders.push_back(&colliderBox);
    sceneColliders.push_back(&colliderBall1);
    sceneColliders.push_back(&colliderBall2);
    sceneColliders.push_back(&colliderBall3);
    sceneColliders.push_back(&colliderPillar1);
    sceneColliders.push_back(&colliderPillar2);
    sceneColliders.push_back(&colliderPillar3);
    sceneColliders.push_back(&colliderPlatform1);
    sceneColliders.push_back(&colliderPlatform2);
    sceneColliders.push_back(&colliderWallLeft);
    sceneColliders.push_back(&colliderWallRight);
}


void SceneFountain::reset()
{
    // update values from UI
    updateSimParams();

    // reset random seed
    Random::seed(1337);

    // erase all particles
    fGravity->clearInfluencedParticles();
    system.deleteParticles();
    deadParticles.clear();
}


void SceneFountain::updateSimParams()
{
    if (integrator) delete integrator;
    integrator = createIntegrator(widget->getIntegratorType());

    // get gravity from UI and update force
    double g = widget->getGravity();
    fGravity->setAcceleration(Vec3(0, -g, 0));

    // get other relevant UI values and update simulation params
    kBounce = widget->getKBounce();
    kFriction = widget->getKFriction();
    maxParticleLife = widget->getMaxParticleLife();
    emitRate = widget->getEmissionRate();
}


void SceneFountain::paint(const Camera& camera) {

    QOpenGLFunctions* glFuncs = nullptr;
    glFuncs = QOpenGLContext::currentContext()->functions();

    shader->bind();

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shader->setUniformValue("ProjMatrix", camProj);
    shader->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 1;
    const QVector3D lightPosWorld[numLights] = {QVector3D(100,500,100)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView.map(lightPosWorld[i]);  // map = matrix * vector
    }
    shader->setUniformValue("numLights", numLights);
    shader->setUniformValueArray("lightPos", lightPosCam, numLights);
    shader->setUniformValueArray("lightColor", lightColor, numLights);

    // draw floor
    vaoFloor->bind();
    QMatrix4x4 modelMat = getPlaneModelMatrix(colliderFloor, 100.0f);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw ramp
    modelMat = getPlaneModelMatrix(colliderRamp, 100.0f);
    shader->setUniformValue("ModelMatrix", modelMat);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw the particles
    vaoSphereL->bind();
    for (const Particle* particle : system.getParticles()) {
        Vec3   p = particle->pos;
        Vec3   c = particle->color;
        double r = particle->radius;

        modelMat = QMatrix4x4();
        modelMat.translate(p[0], p[1], p[2]);
        modelMat.scale(r);
        shader->setUniformValue("ModelMatrix", modelMat);

        shader->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shader->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereL, GL_UNSIGNED_INT, 0);
    }

    // draw sphere
    vaoSphereH->bind();
    modelMat = getSphereModelMatrix(colliderSphere);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.4f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    // draw box
    vaoCube->bind();
    modelMat = getAABBModelMatrix(colliderBox);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.4f, 0.8f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    vaoSphereH->bind();
    modelMat = getSphereModelMatrix(colliderBall1);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.9f, 0.7f, 0.2f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    modelMat = getSphereModelMatrix(colliderBall2);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.2f, 0.7f, 0.9f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    modelMat = getSphereModelMatrix(colliderBall3);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.9f, 0.2f, 0.7f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    vaoCube->bind();
    modelMat = getAABBModelMatrix(colliderPillar1);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.6f, 0.4f, 0.3f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderPillar2);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.3f, 0.6f, 0.4f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderPillar3);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.4f, 0.3f, 0.6f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderPlatform1);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.7f, 0.5f, 0.3f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderPlatform2);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.5f, 0.7f, 0.3f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderWallLeft);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.5f, 0.5f, 0.5f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    modelMat = getAABBModelMatrix(colliderWallRight);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.5f, 0.5f, 0.5f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    vaoCube->release();
    shader->release();
}


void SceneFountain::update(double dt) {

    // emit new particles, reuse dead ones if possible
    int emitParticles = std::max(1, int(std::round(emitRate * dt)));
    
    for (int i = 0; i < emitParticles; i++) {
        Particle* p;
        if (!deadParticles.empty()) {
            // reuse one dead particle
            p = deadParticles.front();
            deadParticles.pop_front();
        }
        else {
            // create new particle
            p = new Particle();
            system.addParticle(p);

            // don't forget to add particle to forces that affect it
            fGravity->addInfluencedParticle(p);
        }

        p->color = Vec3(153/255.0, 217/255.0, 234/255.0);
        p->radius = 1.0;
        p->life = maxParticleLife;

        double x = Random::get(-20.0, 20.0);
        double y = 0;
        double z = Random::get(-20.0, 20.0);
        p->pos = Vec3(x, y, z) + fountainPos;
        p->vel = Vec3(0,0,0);
    }

    // integration step
    Vecd ppos = system.getPositions();
    integrator->step(system, dt);
    system.setPreviousPositions(ppos);

    // collisions
    Collision colInfo;
    for (Particle* p : system.getParticles()) {
        for (Collider* c : sceneColliders) {
            if (c->testCollision(p, colInfo)) {
                c->resolveCollision(p, colInfo, kBounce, kFriction);
            }
        }
    }


    //system.buildSpatialHashing(2);

    //const std::vector<Particle*>& particles = system.getParticles();
    //for (int i = 0; i < particles.size(); ++i) {
    //    Particle* p1 = particles[i];

    //    // Get neighbors using spatial hashing (only check nearby particles)
    //    std::vector<int> neighbors = system.getNeighbors(i);

    //    for (int j : neighbors) {
    //        // Only process each pair once (avoid duplicate checks)
    //        if (j <= i) continue;

    //        Particle* p2 = particles[j];

    //        // Create temporary sphere collider for particle 1
    //        ColliderSphere sphereP1;
    //        sphereP1.setCenter(p1->prevPos);
    //        sphereP1.setRadius(p1->radius);

    //        // Test collision between p1 and p2
    //        Collision colInfo1;
    //        if (sphereP1.testCollision(p2, colInfo1)) {
    //            // Resolve collision for p2 with normal pointing from p1 to p2
    //            sphereP1.resolveCollision(p2, colInfo1, kBounce, kFriction);

    //            // Resolve collision for p1 with negated normal (Newton's 3rd law)
    //            Collision colInfo2 = colInfo1;
    //            colInfo2.normal = -colInfo1.normal;
    //            sphereP1.resolveCollision(p1, colInfo2, kBounce, kFriction);
    //        }
    //    }
    //}

    // check dead particles
    for (Particle* p : system.getParticles()) {
        if (p->life > 0) {
            p->life -= dt;
            if (p->life < 0) {
                deadParticles.push_back(p);
            }
        }
    }
}

//void SceneFountain::mousePressed(const QMouseEvent* e, const Camera&)
//{
//    mouseX = e->pos().x();
//    mouseY = e->pos().y();
//}

void SceneFountain::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    Scene::mouseMoved(e, cam);

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    // example
    if (e->buttons() & Qt::RightButton) {
        if (!e->modifiers()) {
            // move fountain
            fountainPos += disp;
        }
    }
}
