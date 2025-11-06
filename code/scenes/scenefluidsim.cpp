#include "scenefluidsim.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>
#include <iostream>
#include <QOpenGLVersionFunctionsFactory>





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
    // Fan positioned on the left wall (X = -50), blowing towards +X
    Vec3 windCenter(-55, 5, 40);     // center of outlet face (where air comes out)
    double windWidth = 10.0;         // cross-section width (Z axis)
    double windHeight = 10.0;        // cross-section height (Y axis)
    double windLength = 70.0;        // how far the air reaches (X axis)
    double fanThickness = 2.0;       // thickness of the fan box itself
    
    // Calculate fan center and size for rendering
    fanCenter = windCenter - Vec3(fanThickness, 0, 0);  // wind blows in +X, so fan is at -X
    fanHalfSize = Vec3(fanThickness * 0.5, windHeight * 0.5, windWidth * 0.5);
    
    // Set up wind force
    fWind->setCenter(windCenter);
    fWind->setDimensions(windWidth, windHeight);
    fWind->setDirection(Vec3(1, 0, 0));  // blow towards +X (direction snaps to axis)
    fWind->setLength(windLength);
    fWind->setStrength(100.0);
    fWind->setFrequency(0.0);
    fWind->setTimePointer(system.getTimePointer());
#endif

    //sceneColliders.push_back(&colliderFloor);
    //sceneColliders.push_back(&colliderRamp);
    sceneColliders.push_back(&colliderSphere);
    sceneColliders.push_back(&colliderBox);
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

    QOpenGLContext* context = QOpenGLContext::currentContext();

    QOpenGLFunctions_3_3_Core* glFuncs = QOpenGLVersionFunctionsFactory::get<QOpenGLFunctions_3_3_Core>(context);
    glFuncs->initializeOpenGLFunctions();

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

    QMatrix4x4 modelMat;
    //// draw floor
    //vaoFloor->bind();
    //modelMat.scale(100, 1, 100);
    //shader->setUniformValue("ModelMatrix", modelMat);
    //shader->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    //shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    //shader->setUniformValue("matshin", 0.0f);
    //glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    //// draw ramp
    //modelMat = QMatrix4x4();
    //modelMat.rotate(30.0, QVector3D(1, 0, 0));
    //modelMat.translate(0, -6, 0);
    //modelMat.scale(100, 1, 100);
    //modelMat.translate(0, 0, -1);
    //shader->setUniformValue("ModelMatrix", modelMat);
    //glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

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

#if ENABLE_WIND
    // draw wind/fan box (the physical fan)
    if (fWind && windEnabled) {
        vaoCube->bind();
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
    
    // Draw colliderBounds with semi-transparent wireframe effect
    // We draw it with a very low alpha so it's barely visible but shows the boundaries
    vaoCube->bind();
    modelMat = getAABBModelMatrix(colliderBounds);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.5f, 0.5f, 0.5f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    
    glFuncs->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glFuncs->glDrawElements(GL_TRIANGLES, 3*2*6, GL_UNSIGNED_INT, 0);

    glFuncs->glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
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
        for (Collider* c : sceneColliders) {
            while (c->testCollision(p, colInfo)) {
                c->resolveCollision(p, colInfo, kBounce, kFriction);
            }
        }

        while (colliderBounds.testCollision(p, colInfo)) {
            colliderBounds.resolveCollision(p, colInfo, kBounce, kFriction);
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

//void SceneFluidSim::mousePressed(const QMouseEvent* e, const Camera&)
//{
//    mouseX = e->pos().x();
//    mouseY = e->pos().y();
//}

void SceneFluidSim::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    Scene::mouseMoved(e, cam);
    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    if (e->buttons() & Qt::RightButton) {
        if (!e->modifiers()) {
            // move fountain
            fountainPos += disp;
        }
        else if (e->modifiers() & Qt::AltModifier) {
            // move bounds
            colliderBounds.move(disp);
            
#if ENABLE_WIND
            // When bounds move, update fan position to keep it on the wall
            // Get the new wall position (min X of bounds)
            Vec3 boundsMin = colliderBounds.getMin();
            Vec3 boundsMax = colliderBounds.getMax();
            
            // Update fan X position to be on the left wall (min X)
            // Keep Y and Z relative positions within bounds
            double wallX = boundsMin[0];
            double windOutletX = wallX - fanHalfSize[0] * 2.0; // outlet is at wallX, fan is behind
            
            // Clamp fan position to be within bounds height and depth
            double minY = boundsMin[1] + fanHalfSize[1];
            double maxY = boundsMax[1] - fanHalfSize[1];
            double minZ = boundsMin[2] + fanHalfSize[2];
            double maxZ = boundsMax[2] - fanHalfSize[2];
            
            fanCenter[0] = windOutletX - fanHalfSize[0];
            fanCenter[1] = std::max(minY, std::min(maxY, fanCenter[1]));
            fanCenter[2] = std::max(minZ, std::min(maxZ, fanCenter[2]));
            
            // Update wind outlet position
            Vec3 windCenter = fanCenter + Vec3(fanHalfSize[0] * 2.0, 0, 0);
            fWind->setCenter(windCenter);
#endif
        }
#if ENABLE_WIND
        else if (e->modifiers() & Qt::ShiftModifier) {
            // move fan on the wall plane (YZ plane, X is fixed to wall)
            // Project displacement to YZ plane
            Vec3 planeDisp = disp;
            planeDisp[0] = 0;  // no movement in X direction
            
            fanCenter += planeDisp;
            
            // Clamp fan position to be within bounds
            Vec3 boundsMin = colliderBounds.getMin();
            Vec3 boundsMax = colliderBounds.getMax();
            
            double wallX = boundsMin[0];
            double minY = boundsMin[1] + fanHalfSize[1];
            double maxY = boundsMax[1] - fanHalfSize[1];
            double minZ = boundsMin[2] + fanHalfSize[2];
            double maxZ = boundsMax[2] - fanHalfSize[2];
            
            // Keep fan on the wall and within bounds
            fanCenter[0] = wallX - fanHalfSize[0];
            fanCenter[1] = std::max(minY, std::min(maxY, fanCenter[1]));
            fanCenter[2] = std::max(minZ, std::min(maxZ, fanCenter[2]));
            
            // Update wind force position
            Vec3 windCenter = fanCenter + Vec3(fanHalfSize[0] * 2.0, 0, 0);
            fWind->setCenter(windCenter);
        }
#endif
    }
}
