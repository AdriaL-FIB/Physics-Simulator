#include "scenefluidsim.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>
#include <iostream>



SceneFluidSim::SceneFluidSim() {
    widget = new WidgetFluidSim();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
    particleRadius = 1.0;
    widget->setSolverTypes(getIntegratorNames());
    widget->setIntegratorType(1);
#if ENABLE_WIND
    windEnabled = false;  // Wind is created but disabled by default
#endif
}


SceneFluidSim::~SceneFluidSim() {
    if (widget)     delete widget;
    if (shader)     delete shader;
    if (vaoFloor)   delete vaoFloor;
    if (vaoSphereH) delete vaoSphereH;
    if (vaoSphereL) delete vaoSphereL;
    if (vaoCube)    delete vaoCube;
    if (fGravity)   delete fGravity;
#if ENABLE_WIND
    if (fWind)      delete fWind;
#endif
}


void SceneFluidSim::initialize() {
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

    fSPH = new ForceSPH(&system, 3.0 * particleRadius);
    system.addForce(fSPH);

#if ENABLE_WIND
    // create wind/fan force
    fWind = new ForceWind();
    // Note: Force is created but not added to system yet (controlled by windEnabled)
#endif

    // scene description
    fountainPos = Vec3(0, 80, 0);    
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);    
    colliderRamp.setPlane(Vec3(0, std::sqrt(3.0)/2.0, 0.5), 6);
    colliderSphere.setCenter(Vec3(0,0,0));
    colliderSphere.setRadius(20);
    colliderBox.setFromBounds(Vec3(30,0,20), Vec3(50,10,60));

    colliderBounds.setFromBounds(Vec3(-50, 0, -50), Vec3(50, 100, 50));
    colliderBounds.setHollow(true);

#if ENABLE_WIND
    // Setup wind/fan at ground level in positive Z corner
    // Bounds: [-50,0,-50] to [50,100,50]
    // Center at corner (-55, 5, 40), blowing towards center (positive X)
    Vec3 windCenter(-55, 5, 40);     // center of outlet face (where air comes out)
    double windWidth = 10.0;         // cross-section width (Z axis)
    double windHeight = 10.0;        // cross-section height (Y axis)
    double windLength = 70.0;        // how far the air reaches (X axis)
    fWind->setCenter(windCenter);
    fWind->setDimensions(windWidth, windHeight);
    fWind->setDirection(Vec3(1, 0, 0));  // blow towards +X (direction snaps to axis)
    fWind->setLength(windLength);
    fWind->setStrength(100.0);
    fWind->setFrequency(0.0);
    fWind->setTimePointer(system.getTimePointer());
#endif

    //colliders.push_back(&colliderFloor);
    //colliders.push_back(&colliderRamp);
    colliders.push_back(&colliderSphere);
    //colliders.push_back(&colliderBox);
    colliders.push_back(&colliderBounds);
}


void SceneFluidSim::reset()
{
    // update values from UI
    updateSimParams();

    // reset random seed
    Random::seed(1337);

    // erase all particles
    fGravity->clearInfluencedParticles();
    fSPH->clearInfluencedParticles();
#if ENABLE_WIND
    if (fWind) fWind->clearInfluencedParticles();
#endif
    system.deleteParticles();
    deadParticles.clear();
}


void SceneFluidSim::updateSimParams()
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
    particleRadius = widget->getParticleRadius();
    particleMass = widget->getParticleMass();
    emitParticles = widget->emitParticles();

    fSPH->setRestDensity(widget->getRestDensity());
    fSPH->setViscosityConstant(widget->getViscosityConstant());

#if ENABLE_WIND
    // Control wind force activation from UI
    windEnabled = widget->getEnableWind();
    
    bool windInSystem = false;
    for (Force* f : system.getSystemForces()) {
        if (f == fWind) {
            windInSystem = true;
            break;
        }
    }
    
    if (windEnabled && !windInSystem) {
        // Add wind to system
        system.addForce(fWind);
        std::cout << "WIND ADDED - Total forces: " << system.getNumForces() << std::endl;
    } else if (!windEnabled && windInSystem) {
        // Remove wind from system
        system.removeForce(fWind);
        std::cout << "WIND REMOVED - Total forces: " << system.getNumForces() << std::endl;
    }
#endif
}


void SceneFluidSim::paint(const Camera& camera) {

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
    QMatrix4x4 modelMat;
    modelMat.scale(100, 1, 100);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw ramp
    modelMat = QMatrix4x4();
    modelMat.rotate(30.0, QVector3D(1, 0, 0));
    modelMat.translate(0, -6, 0);
    modelMat.scale(100, 1, 100);
    modelMat.translate(0, 0, -1);
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
    Vec3 cc = colliderSphere.getCenter();
    modelMat = QMatrix4x4();
    modelMat.translate(cc[0], cc[1], cc[2]);
    modelMat.scale(colliderSphere.getRadius());
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.4f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereH, GL_UNSIGNED_INT, 0);

    // draw box
    vaoCube->bind();
    cc = colliderBox.getCenter();
    Vec3 hs = 0.5*colliderBox.getSize();
    modelMat = QMatrix4x4();
    modelMat.translate(cc[0], cc[1], cc[2]);
    modelMat.scale(hs[0], hs[1], hs[2]);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.4f, 0.8f, 0.4f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

#if ENABLE_WIND
    // draw wind/fan box (the physical fan, not the air)
    if (fWind && windEnabled) {
        Vec3 outletCenter = fWind->getCenter();      // where air comes out
        Vec3 windDir = fWind->getDirection();
        double windWidth = fWind->getWidth();
        double windHeight = fWind->getHeight();
        
        // The fan box should be behind the outlet, opposite to wind direction
        // Fan dimensions: same height and width as air region, but fixed thickness
        double fanThickness = 2.0; // fixed depth of the fan box
        
        // Position the fan box behind the outlet
        Vec3 fanCenter = outletCenter - windDir * fanThickness;
        
        // Fan half-size: thickness in the direction of flow, same width/height as air cross-section
        Vec3 fanHalfSize;
        // Determine which axis is flow direction and set dimensions accordingly
        if (std::abs(windDir[0]) > 0.5) {
            // Flow along X axis
            fanHalfSize = Vec3(fanThickness, windHeight * 0.5, windWidth * 0.5);
        } else if (std::abs(windDir[1]) > 0.5) {
            // Flow along Y axis
            fanHalfSize = Vec3(windWidth * 0.5, fanThickness, windHeight * 0.5);
        } else {
            // Flow along Z axis
            fanHalfSize = Vec3(windWidth * 0.5, windHeight * 0.5, fanThickness);
        }
        
        modelMat = QMatrix4x4();
        modelMat.translate(fanCenter[0], fanCenter[1], fanCenter[2]);
        modelMat.scale(fanHalfSize[0], fanHalfSize[1], fanHalfSize[2]);
        shader->setUniformValue("ModelMatrix", modelMat);
        
        // Dark gray/metallic color for the fan casing
        shader->setUniformValue("matdiff", 0.3f, 0.3f, 0.35f);
        shader->setUniformValue("matspec", 0.6f, 0.6f, 0.6f);
        shader->setUniformValue("matshin", 50.0f);
        
        glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);
    }
#endif
    
    vaoCube->release();
    shader->release();
}


void SceneFluidSim::update(double dt) {

    // emit new particles, reuse dead ones if possible
    if (emitParticles) {
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
                fSPH->addInfluencedParticle(p);
#if ENABLE_WIND
                if (fWind) fWind->addInfluencedParticle(p);
#endif
            }
    
            p->color = Vec3(153/255.0, 217/255.0, 234/255.0);
            p->radius = particleRadius;
            p->life = maxParticleLife;
            p->mass = particleMass;
    
            double x = Random::get(-20.0, 20.0);
            double y = 0;
            double z = Random::get(-20.0, 20.0);
            p->pos = Vec3(x, y, z) + fountainPos;
            p->vel = Vec3(0,0,0);
        }
    }

    // integration step
    Vecd ppos = system.getPositions();
    integrator->step(system, dt);
    system.setPreviousPositions(ppos);

    // collisions
    Collision colInfo;
    for (Particle* p : system.getParticles()) {
        for (Collider* c : colliders) {
            while (c->testCollision(p, colInfo)) {
                c->resolveCollision(p, colInfo, kBounce, kFriction);
            }
        }
    }

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

void SceneFluidSim::mousePressed(const QMouseEvent* e, const Camera&)
{
    mouseX = e->pos().x();
    mouseY = e->pos().y();
}

void SceneFluidSim::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - mouseX;
    int dy = e->pos().y() - mouseY;
    mouseX = e->pos().x();
    mouseY = e->pos().y();

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    // example
    if (e->buttons() & Qt::RightButton) {
        if (!e->modifiers()) {
            // move fountain
            fountainPos += disp;
        }
        else if (e->modifiers() & Qt::ShiftModifier){
            // move box
            colliderBox.setFromCenterSize(
                        colliderBox.getCenter() + disp,
                        colliderBox.getSize());
        }
    }
}
