#include "forces.h"
#include "math.h"
#include "particlesystem.h"
#include <iostream>
#include <cmath>

void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        p->force += p->mass * this->acceleration;
    }
}

void ForceDrag::apply() {
    for (Particle* p : particles) {
        // Strokes drag:
        p->force += -this->klinear * p->vel;

        // Newton drag:
        //p->force += -this->kquadratic * p->vel.cwiseAbs() * (p->vel);
    }
}

void ForceSpring::apply() {
    if (particles.size() < 2) return;
    Particle* p1 = getParticle1();
    Particle* p2 = getParticle2();

    double norm = (p2->pos - p1->pos).norm();
    Vec3 udir = (p2->pos - p1->pos)/norm;

    Vec3 F = (this->ks * (norm - this->L) + this->kd * (p2->vel - p1->vel).dot(udir)) * (udir);

    p1->force += F;
    p2->force -= F;
}

void ForceGravitation::apply() {
    const Particle* p1 = attractor;
    for (unsigned int j = 0; j < particles.size(); ++j) {
        Particle* p2 = particles[j];

        Vec3 diff = p2->pos - p1->pos;
        double distSquared = diff.squaredNorm();

        Vec3 F = (G * p1->mass * p2->mass) / distSquared * diff / diff.norm();

        F *= 2 / (1 + exp(-a * distSquared / (b * b))) - 1;

        p2->force -= F;
    }
}

void ForceSPH::apply()
{
    system->buildSpatialHashing(h);

    int N = particles.size();

    densities.resize(N);
    pressures.resize(N);

    // DEBUG: Track specific particle across frames
    static int debugCounter = 0;
    static int trackedParticleID = -1;
    static int framesSinceTracking = 0;
    debugCounter++;
    
    bool printDebug = (debugCounter % 120 == 0); // Every 120 calls (~2 seconds)
    bool printTracking = (trackedParticleID >= 0 && trackedParticleID < N);
    
    int particleWithMostNeighbors = -1;
    int maxNeighborCount = 0;

    // Compute density
    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        double density = p->mass * kernelPoly6(0);

        std::vector<int> nbs = system->getNeighbors(p);
        
        // Track particle with most neighbors for detailed analysis
        if (nbs.size() > maxNeighborCount) {
            maxNeighborCount = nbs.size();
            particleWithMostNeighbors = i;
        }
        
        for (int j = 0; j < nbs.size(); ++j) {
            Particle* pj = particles[nbs[j]];
            density += pj->mass * kernelPoly6((pj->pos - p->pos).squaredNorm());
        }
        densities[i] = density;
    }

    // Compute pressure
    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        pressures[i] = computePressure(densities[i]);
    }

    // DETAILED ANALYSIS of the most crowded particle
    if (printDebug && particleWithMostNeighbors >= 0) {
        Particle* pCrowded = particles[particleWithMostNeighbors];
        std::vector<int> nbsCrowded = system->getNeighbors(pCrowded);
        
        std::cout << "\n========================================================" << std::endl;
        std::cout << "      SPH ANALYSIS - Most Crowded Particle" << std::endl;
        std::cout << "========================================================" << std::endl;
        std::cout << "Particle ID: " << particleWithMostNeighbors << std::endl;
        std::cout << "Neighbors: " << nbsCrowded.size() << std::endl;
        std::cout << "Rest Density (p0): " << p0 << std::endl;
        std::cout << "\n--- DENSITY ANALYSIS ---" << std::endl;
        std::cout << "Calculated Density: " << densities[particleWithMostNeighbors] << std::endl;
        std::cout << "Ratio to rest (rho/p0): " << (densities[particleWithMostNeighbors]/p0) << "x" << std::endl;
        
        // Analyze neighbor distances
        int distBins[5] = {0, 0, 0, 0, 0}; // 0-0.2h, 0.2-0.4h, 0.4-0.6h, 0.6-0.8h, 0.8-1.0h
        double minDist = 1e10, maxDist = 0;
        double contribByRange[5] = {0, 0, 0, 0, 0};
        
        for (int j = 0; j < nbsCrowded.size(); ++j) {
            Particle* pj = particles[nbsCrowded[j]];
            double rSqr = (pj->pos - pCrowded->pos).squaredNorm();
            double r = std::sqrt(rSqr);
            double normalizedR = r / h;
            
            if (r < minDist) minDist = r;
            if (r > maxDist) maxDist = r;
            
            double contrib = pj->mass * kernelPoly6(rSqr);
            
            // Bin by distance
            if (normalizedR < 0.2) { distBins[0]++; contribByRange[0] += contrib; }
            else if (normalizedR < 0.4) { distBins[1]++; contribByRange[1] += contrib; }
            else if (normalizedR < 0.6) { distBins[2]++; contribByRange[2] += contrib; }
            else if (normalizedR < 0.8) { distBins[3]++; contribByRange[3] += contrib; }
            else { distBins[4]++; contribByRange[4] += contrib; }
        }
        
        double totalNeighborContrib = 0;
        for (int i = 0; i < 5; i++) totalNeighborContrib += contribByRange[i];
        
        std::cout << "\nMin distance to neighbor: " << minDist << " (" << (minDist/h) << "*h)" << std::endl;
        std::cout << "Max distance to neighbor: " << maxDist << " (" << (maxDist/h) << "*h)" << std::endl;
        
        std::cout << "\n--- NEIGHBOR DISTRIBUTION ---" << std::endl;
        std::cout << "Distance Range | Count | Density Contribution" << std::endl;
        std::cout << "0.0-0.2*h      | " << distBins[0] << " | " << contribByRange[0] << " (" << (100.0*contribByRange[0]/totalNeighborContrib) << "%)" << std::endl;
        std::cout << "0.2-0.4*h      | " << distBins[1] << " | " << contribByRange[1] << " (" << (100.0*contribByRange[1]/totalNeighborContrib) << "%)" << std::endl;
        std::cout << "0.4-0.6*h      | " << distBins[2] << " | " << contribByRange[2] << " (" << (100.0*contribByRange[2]/totalNeighborContrib) << "%)" << std::endl;
        std::cout << "0.6-0.8*h      | " << distBins[3] << " | " << contribByRange[3] << " (" << (100.0*contribByRange[3]/totalNeighborContrib) << "%)" << std::endl;
        std::cout << "0.8-1.0*h      | " << distBins[4] << " | " << contribByRange[4] << " (" << (100.0*contribByRange[4]/totalNeighborContrib) << "%)" << std::endl;
        
        std::cout << "\n--- PRESSURE ANALYSIS ---" << std::endl;
        std::cout << "Pressure: " << pressures[particleWithMostNeighbors] << std::endl;
        std::cout << "Pressure Formula: P = k * (rho - p0)" << std::endl;
        std::cout << "  where k = 50" << std::endl;
        std::cout << "Calculated: 50 * (" << densities[particleWithMostNeighbors] << " - " << p0 << ") = " << pressures[particleWithMostNeighbors] << std::endl;
    }

    for(int i = 0; i < N; ++i) {
        Particle* p = particles[i];

        Vec3 aPressure(0,0,0);
        Vec3 aViscosity(0,0,0);
        double minDistToNeighbor = h; // Track closest neighbor distance

        // TODO: optimize
        std::vector<int> nbs = system->getNeighbors(p);
        for (int j = 0; j < nbs.size(); ++j) {
            int neighborIndx = nbs[j];
            Particle* pj = particles[neighborIndx];

            Vec3 r = pj->pos - p->pos;
            double dist = r.norm();
            if (dist < minDistToNeighbor) {
                minDistToNeighbor = dist;
            }

            double rho1 = densities[i];
            double rho2 = densities[neighborIndx];
            
            double Pij = -pj->mass*(pressures[i]/(rho1*rho1) + pressures[neighborIndx]/(rho2*rho2));

            aPressure += Pij * kernelGradientSpiky(r);

            Vec3 Vij = mu * pj->mass * (pj->vel - p->vel)/(rho1*rho2);

            aViscosity += Vij * kernelLaplacianViscosity(r.norm());
        }

        p->force += p->mass * (aPressure + aViscosity);
        
        // FORCE ANALYSIS for the most crowded particle
        if (printDebug && i == particleWithMostNeighbors) {
            std::cout << "\n--- FORCE ANALYSIS ---" << std::endl;
            std::cout << "Pressure Acceleration: (" << aPressure[0] << ", " << aPressure[1] << ", " << aPressure[2] << ")" << std::endl;
            std::cout << "Pressure Accel Magnitude: " << aPressure.norm() << std::endl;
            std::cout << "Viscosity Acceleration: (" << aViscosity[0] << ", " << aViscosity[1] << ", " << aViscosity[2] << ")" << std::endl;
            std::cout << "Viscosity Accel Magnitude: " << aViscosity.norm() << std::endl;
            
            Vec3 totalAccel = aPressure + aViscosity;
            std::cout << "Total SPH Acceleration: " << totalAccel.norm() << std::endl;
            std::cout << "Total SPH Force: " << (p->mass * totalAccel.norm()) << std::endl;
            
            // ============================================================
            // ANÁLISIS DETALLADO: ¿Por qué no se separan?
            // ============================================================
            
            std::cout << "\n========================================" << std::endl;
            std::cout << "  DIAGNÓSTICO DE SEPARACIÓN" << std::endl;
            std::cout << "========================================" << std::endl;
            
            // 1. Estado de la partícula
            std::cout << "\n[1] ESTADO ACTUAL:" << std::endl;
            std::cout << "  - Vecinos: " << nbs.size() << std::endl;
            std::cout << "  - Densidad: " << densities[i] << " (ratio: " << (densities[i]/p0) << "x rest)" << std::endl;
            std::cout << "  - Presion: " << pressures[i] << std::endl;
            std::cout << "  - Velocidad: " << p->vel.norm() << " m/s" << std::endl;
            std::cout << "  - Vecino mas cercano: " << minDistToNeighbor << " (overlap: " << (2*p->radius - minDistToNeighbor) << ")" << std::endl;
            
            // 2. Fuerzas externas (antes de SPH)
            Vec3 externalForce = p->force - p->mass * (aPressure + aViscosity);
            Vec3 externalAccel = externalForce / p->mass;
            std::cout << "\n[2] FUERZAS EXTERNAS (gravedad, viento, etc):" << std::endl;
            std::cout << "  - Aceleracion externa: (" << externalAccel[0] << ", " << externalAccel[1] << ", " << externalAccel[2] << ")" << std::endl;
            std::cout << "  - Magnitud: " << externalAccel.norm() << std::endl;
            
            // 3. Fuerzas SPH
            std::cout << "\n[3] FUERZAS SPH:" << std::endl;
            std::cout << "  - Aceleracion presion: (" << aPressure[0] << ", " << aPressure[1] << ", " << aPressure[2] << ")" << std::endl;
            std::cout << "  - Magnitud presion: " << aPressure.norm() << std::endl;
            std::cout << "  - Aceleracion viscosidad: " << aViscosity.norm() << std::endl;
            
            totalAccel = externalAccel + aPressure + aViscosity;
            std::cout << "\n[4] ACELERACION TOTAL:" << std::endl;
            std::cout << "  - Total: (" << totalAccel[0] << ", " << totalAccel[1] << ", " << totalAccel[2] << ")" << std::endl;
            std::cout << "  - Magnitud: " << totalAccel.norm() << std::endl;
            
            // 5. Análisis crítico
            std::cout << "\n[5] ANALISIS CRITICO:" << std::endl;
            
            // Check if pressure force is being dominated by external forces
            double ratio = aPressure.norm() / (externalAccel.norm() + 0.001);
            std::cout << "  - Ratio Presion/Externa: " << ratio << std::endl;
            
            if (ratio < 1.0) {
                std::cout << "  >>> PROBLEMA: Fuerza externa DOMINA sobre presion! <<<" << std::endl;
                std::cout << "      Gravedad/viento estan empujando mas fuerte que SPH." << std::endl;
            }
            
            // Check if forces point in opposite directions (cancelation)
            double dotProduct = aPressure.normalized().dot(externalAccel.normalized());
            std::cout << "  - Direccion presion vs externa (dot): " << dotProduct << std::endl;
            if (dotProduct < -0.5) {
                std::cout << "  >>> PROBLEMA: Fuerzas apuntan en DIRECCIONES OPUESTAS! <<<" << std::endl;
                std::cout << "      SPH empuja hacia arriba, gravedad/viento empujan hacia abajo." << std::endl;
            }
            
            // Check velocity - are particles even moving?
            if (p->vel.norm() < 0.1) {
                std::cout << "  >>> PROBLEMA: Particula casi ESTATICA (vel < 0.1)! <<<" << std::endl;
                std::cout << "      Las fuerzas se estan CANCELANDO completamente." << std::endl;
            }
            
            // Check if pressure formula is effective
            double pressureCoeff = pressures[i] / (densities[i] * densities[i]);
            std::cout << "  - Coeficiente P/(rho^2): " << pressureCoeff << std::endl;
            if (pressureCoeff < 0.1) {
                std::cout << "  >>> PROBLEMA: P/(rho^2) es muy PEQUEÑO! <<<" << std::endl;
                std::cout << "      Densidad alta (" << densities[i] << ") hace que rho^2 = " << (densities[i]*densities[i]) << std::endl;
                std::cout << "      Esto DEBILITA enormemente la fuerza de presion!" << std::endl;
            }
            
            std::cout << "\n========================================\n" << std::endl;
            
            // Start tracking this particle if it has many neighbors
            if (nbs.size() > 100) {
                trackedParticleID = i;
                framesSinceTracking = 0;
                std::cout << "\n*** TRACKING STARTED for Particle " << i << " ***" << std::endl;
                std::cout << "    Will print updates every frame for next 100 frames\n" << std::endl;
            }
        }
        
        // ============================================================
        // FRAME-BY-FRAME TRACKING of problematic particle
        // ============================================================
        if (printTracking && i == trackedParticleID && framesSinceTracking < 100) {
            framesSinceTracking++;
            
            // Calculate center of mass of neighbors to check symmetry
            Vec3 neighborCenterOfMass(0, 0, 0);
            for (int j : nbs) {
                neighborCenterOfMass += particles[j]->pos;
            }
            if (nbs.size() > 0) {
                neighborCenterOfMass /= nbs.size();
            }
            
            // Track previous position/velocity for comparison
            static Vec3 prevPos = p->pos;
            static Vec3 prevVel = p->vel;
            static Vec3 prevCOM = neighborCenterOfMass;
            
            Vec3 displacement = p->pos - prevPos;
            Vec3 velocityChange = p->vel - prevVel;
            
            std::cout << "\n[FRAME " << framesSinceTracking << "] Particle " << i << " Status:" << std::endl;
            std::cout << "  Pos: (" << p->pos[0] << ", " << p->pos[1] << ", " << p->pos[2] << ")" << std::endl;
            std::cout << "  Vel: (" << p->vel[0] << ", " << p->vel[1] << ", " << p->vel[2] << ") mag=" << p->vel.norm() << std::endl;
            std::cout << "  Neighbors: " << nbs.size() << " | Density: " << densities[i] << " | Pressure: " << pressures[i] << std::endl;
            
            Vec3 externalForce = p->force - p->mass * (aPressure + aViscosity);
            Vec3 externalAccel = externalForce / p->mass;
            
            std::cout << "  SPH Accel: (" << aPressure[0] << ", " << aPressure[1] << ", " << aPressure[2] << ") mag=" << aPressure.norm() << std::endl;
            std::cout << "  Ext Accel: (" << externalAccel[0] << ", " << externalAccel[1] << ", " << externalAccel[2] << ") mag=" << externalAccel.norm() << std::endl;
            
            Vec3 totalAccel = externalAccel + aPressure + aViscosity;
            std::cout << "  Tot Accel: (" << totalAccel[0] << ", " << totalAccel[1] << ", " << totalAccel[2] << ") mag=" << totalAccel.norm() << std::endl;
            std::cout << "  Closest neighbor: " << minDistToNeighbor << " (overlap: " << (2*p->radius - minDistToNeighbor) << ")" << std::endl;
            
            // Analyze neighbor distribution symmetry
            Vec3 directionAwayFromCOM = (p->pos - neighborCenterOfMass).normalized();
            double dotSPHvsAway = aPressure.normalized().dot(directionAwayFromCOM);
            
            std::cout << "\n  --- SPATIAL ANALYSIS ---" << std::endl;
            std::cout << "  Neighbor COM: (" << neighborCenterOfMass[0] << ", " << neighborCenterOfMass[1] << ", " << neighborCenterOfMass[2] << ")" << std::endl;
            std::cout << "  Distance to COM: " << (p->pos - neighborCenterOfMass).norm() << std::endl;
            std::cout << "  Direction AWAY from COM: (" << directionAwayFromCOM[0] << ", " << directionAwayFromCOM[1] << ", " << directionAwayFromCOM[2] << ")" << std::endl;
            std::cout << "  SPH Accel direction: (" << aPressure.normalized()[0] << ", " << aPressure.normalized()[1] << ", " << aPressure.normalized()[2] << ")" << std::endl;
            std::cout << "  SPH vs AWAY alignment: " << dotSPHvsAway << " (1=perfect repulsion, -1=attraction)" << std::endl;
            
            if (framesSinceTracking > 1) {
                std::cout << "\n  --- MOVEMENT ANALYSIS (since last frame) ---" << std::endl;
                std::cout << "  Displacement: (" << displacement[0] << ", " << displacement[1] << ", " << displacement[2] << ") mag=" << displacement.norm() << std::endl;
                
                // Check if particle is moving away from or towards COM
                if (displacement.norm() > 0.001) {
                    double dotDispVsAway = displacement.normalized().dot(directionAwayFromCOM);
                    std::cout << "  Displacement vs AWAY: " << dotDispVsAway << " (>0 = escaping, <0 = approaching)" << std::endl;
                }
                
                if (p->vel.norm() > 0.001) {
                    double dotVelVsAway = p->vel.normalized().dot(directionAwayFromCOM);
                    std::cout << "  Velocity vs AWAY: " << dotVelVsAway << " (>0 = escaping, <0 = approaching)" << std::endl;
                    
                    if (dotVelVsAway < -0.3) {
                        std::cout << "  >>> WARNING! Particle moving TOWARDS COM despite high pressure! <<<" << std::endl;
                    }
                }
                
                // Distance change to COM
                double prevDistToCOM = (prevPos - prevCOM).norm();
                double currDistToCOM = (p->pos - neighborCenterOfMass).norm();
                double distChange = currDistToCOM - prevDistToCOM;
                std::cout << "  Distance change to COM: " << distChange << " (>0 = separating, <0 = collapsing)" << std::endl;
                
                // Critical analysis
                if (dotSPHvsAway < -0.1) {
                    std::cout << "  >>> BUG! SPH force points TOWARDS COM (attraction instead of repulsion)! <<<" << std::endl;
                }
                if (distChange < -0.01) {
                    std::cout << "  >>> COLLAPSING! Distance to COM decreasing! <<<" << std::endl;
                } else if (distChange > 0.01) {
                    std::cout << "  >>> ESCAPING! Distance to COM increasing! <<<" << std::endl;
                }
            }
            
            // Update previous values for next frame
            prevPos = p->pos;
            prevVel = p->vel;
            prevCOM = neighborCenterOfMass;
            
            // Calculate expected SPH magnitude vs actual
            double avgPressurePerNeighbor = aPressure.norm() / (nbs.size() + 1);
            std::cout << "\n  Avg SPH per neighbor: " << avgPressurePerNeighbor << std::endl;
            if (avgPressurePerNeighbor < 0.1 && nbs.size() > 100) {
                std::cout << "  >>> FUERZAS SE CANCELAN! " << nbs.size() << " vecinos pero contribucion promedio tiny <<<" << std::endl;
            }
            
            // Check if particle is escaping or staying collapsed
            if (nbs.size() < 50) {
                std::cout << "  >>> ESCAPING! Neighbors dropped below 50 <<<" << std::endl;
                trackedParticleID = -1; // Stop tracking
            } else if (p->vel.norm() < 0.5) {
                std::cout << "  >>> STUCK! Velocity very low despite forces <<<" << std::endl;
            } else if (aPressure.norm() < 10.0) {
                std::cout << "  >>> WEAK FORCES! SPH pressure < 10 <<<" << std::endl;
            }
            
            if (framesSinceTracking >= 100) {
                std::cout << "\n*** TRACKING ENDED after 100 frames ***\n" << std::endl;
                trackedParticleID = -1;
            }
        }

        // double pressureNorm = std::clamp(pressures[i] / 100.0, 0.0, 100.0);
        // pressureNorm = (pressureNorm + 1.0) * 0.5; // map from [-1,1] to [0,1]
        
        // // Soft colors: light blue (0.5, 0.7, 1.0) -> light red (1.0, 0.5, 0.5)
        // p->color = Vec3(
        //     0.5 + 0.5 * pressureNorm,      // R: 0.5 -> 1.0
        //     0.7 - 0.2 * pressureNorm,      // G: 0.7 -> 0.5
        //     1.0 - 0.5 * pressureNorm       // B: 1.0 -> 0.5
        // );
    }

    // Print densities and pressure stats (min, max, avg, quaredAverage)
    // if (N > 0) {
    //     double minDensity = densities[0];
    //     double maxDensity = densities[0];
    //     double avgDensity = 0.0;
    //     double sqAvgDensity = 0.0;
        
    //     double minPressure = pressures[0];
    //     double maxPressure = pressures[0];
    //     double avgPressure = 0.0;
    //     double sqAvgPressure = 0.0;
        
    //     for (int i = 0; i < N; ++i) {
    //         // Density stats
    //         minDensity = std::min(minDensity, densities[i]);
    //         maxDensity = std::max(maxDensity, densities[i]);
    //         avgDensity += densities[i];
    //         sqAvgDensity += densities[i] * densities[i];
            
    //         // Pressure stats
    //         minPressure = std::min(minPressure, pressures[i]);
    //         maxPressure = std::max(maxPressure, pressures[i]);
    //         avgPressure += pressures[i];
    //         sqAvgPressure += pressures[i] * pressures[i];
    //     }
        
    //     avgDensity /= N;
    //     sqAvgDensity = std::sqrt(sqAvgDensity / N);
    //     avgPressure /= N;
    //     sqAvgPressure = std::sqrt(sqAvgPressure / N);
        
    //     std::cout << "Density  - Min: " << minDensity << ", Max: " << maxDensity 
    //               << ", Avg: " << avgDensity << ", RMS: " << sqAvgDensity << std::endl;
    //     std::cout << "Pressure - Min: " << minPressure << ", Max: " << maxPressure 
    //               << ", Avg: " << avgPressure << ", RMS: " << sqAvgPressure << std::endl;
    // }
}

double ForceSPH::kernelPoly6(double rSqr) const // TODO: If I already have |r|, use r instead of sqr.
{
    if (rSqr > h * h) return 0;

    double term = h*h - rSqr;

    return poly6Factor * term * term * term;
}

Vec3 ForceSPH::kernelGradientSpiky(Vec3 r) const
{
    double rMag = r.norm();

    if (rMag <= 0.0 || rMag > h) {
        return Vec3(0.0, 0.0, 0.0);
    }

    return -r * spikyGradFactor / rMag * (h - rMag) * (h - rMag);
}

double ForceSPH::kernelLaplacianViscosity(double r) const
{
    if (r <= 0.0 || r > h) {
        return 0.0;
    }
    
    return viscosityLaplacianFactor * (1 - r/h);
}

double ForceSPH::computePressure(double density)
{
    const int k = 50;
    return k * (density - p0);
}

ForceWind::ForceWind() {
    center = Vec3(0, 0, 0);
    width = height = 10.0;
    length = 10.0;
    direction = Vec3(1, 0, 0);
    strength = 10.0;
    frequency = 0.0;
    timePtr = nullptr;
    updateInternals();
}

ForceWind::ForceWind(const Vec3& center, double width, double height, const Vec3& dir, double length, double strength)
    : center(center), width(width), height(height), direction(snapToAxis(dir)),
      length(length), strength(strength), frequency(0.0), timePtr(nullptr) {
    updateInternals();
}

void ForceWind::setDirection(const Vec3& d) {
    direction = snapToAxis(d);
    updateInternals();
}

Vec3 ForceWind::snapToAxis(const Vec3& d) {
    Vec3 dn = d.normalized();
    int maxIdx = 0;
    double maxAbs = std::abs(dn[0]);
    for (int i = 1; i < 3; i++) {
        if (std::abs(dn[i]) > maxAbs) {
            maxAbs = std::abs(dn[i]);
            maxIdx = i;
        }
    }
    Vec3 snapped = Vec3(0, 0, 0);
    snapped[maxIdx] = (dn[maxIdx] > 0) ? 1.0 : -1.0;
    return snapped;
}

void ForceWind::updateInternals() {
    // Since direction is axis-aligned (from snapToAxis), find which axis it is
    // by finding the component with absolute value 1
    int flowAxis = 0;
    if (std::abs(direction[1]) > 0.5) flowAxis = 1;
    else if (std::abs(direction[2]) > 0.5) flowAxis = 2;
    
    // Build AABB: aabbCenter is offset from outlet center by length/2 in flow direction
    aabbCenter = center + direction * (length * 0.5);
    
    // Half-size: length/2 along flow axis, width/2 and height/2 on perpendicular axes
    boxHalfSize = Vec3(0, 0, 0);
    boxHalfSize[flowAxis] = length * 0.5;
    
    // Assign width and height to the two perpendicular axes
    int axis1 = (flowAxis + 1) % 3;
    int axis2 = (flowAxis + 2) % 3;
    boxHalfSize[axis1] = width * 0.5;
    boxHalfSize[axis2] = height * 0.5;
    
    // Precompute inverse half-sizes (optimization: only once when geometry changes)
    invHalfSize = Vec3(0, 0, 0);
    for (int i = 0; i < 3; i++) {
        invHalfSize[i] = (boxHalfSize[i] != 0.0) ? 1.0 / boxHalfSize[i] : 0.0;
    }
}

void ForceWind::apply()
{
    // Get current time from the pointer
    double time = (timePtr != nullptr) ? *timePtr : 0.0;
    
    // Cheap pulsation factor (0 = no pulsation). Keep sine (acceptable cost) for natural pulsation.
    double pulseFactor = 1.0;
    if (frequency > 0.0) {
        pulseFactor = 0.5 * (1.0 + std::sin(2.0 * M_PI * frequency * time));
    }

    // invHalfSize is now precomputed in updateInternals() - no need to compute it here!

    // For each influenced particle, check AABB membership and apply a cheap polynomial falloff
    for (Particle* p : particles) {
        
        // Compute local vector from AABB center to particle
        Vec3 d = p->pos - aabbCenter;

        // Quick AABB test (cheap): if outside on any axis, skip
        if (std::abs(d[0]) > boxHalfSize[0]) continue;
        if (std::abs(d[1]) > boxHalfSize[1]) continue;
        if (std::abs(d[2]) > boxHalfSize[2]) continue;

        // Compute normalized squared distance within box coordinates (no sqrt)
        double nx = d[0] * invHalfSize[0];
        double ny = d[1] * invHalfSize[1];
        double nz = d[2] * invHalfSize[2];
        double nsq = nx*nx + ny*ny + nz*nz; // in [0,3] for corners

        // If outside a unit sphere in normalized coords, we still keep them because we already
        // passed the AABB check. Use nsq to compute a smooth falloff that is cheap.
        // Map nsq to falloff in [0,1] where center -> 1, outer -> 0.
        // Normalize nsq by 3 to map corner to 1.0
        double s = nsq / 3.0; // s in [0,1]
        if (s >= 1.0) continue; // effectively outside influence (cheap)

        double falloff = 1.0 - s; // linear
        // Smooth the falloff with a cubic smoothstep: f = f*f*(3-2*f)
        double f2 = falloff * falloff;
        falloff = f2 * (3.0 - 2.0 * falloff);

        // Final acceleration magnitude (treat strength as acceleration magnitude)
        double accelMag = strength * falloff * pulseFactor;

        // Apply as F = m * a (consistent with other forces)
        p->force += direction * (accelMag * p->mass);
    }
}

