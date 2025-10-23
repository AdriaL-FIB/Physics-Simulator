#include "particlesystem.h"
#include <iostream>

Vecd ParticleSystem::getState() const {
    Vecd state(this->getStateSize());
    for (unsigned int i = 0; i < particles.size(); i++) {
        state[Particle::PhaseDimension*i    ] = particles[i]->pos[0];
        state[Particle::PhaseDimension*i + 1] = particles[i]->pos[1];
        state[Particle::PhaseDimension*i + 2] = particles[i]->pos[2];
        state[Particle::PhaseDimension*i + 3] = particles[i]->vel[0];
        state[Particle::PhaseDimension*i + 4] = particles[i]->vel[1];
        state[Particle::PhaseDimension*i + 5] = particles[i]->vel[2];
    }
    return state;
}

Vecd ParticleSystem::getDerivative() const {
    Vecd deriv(this->getStateSize());
    for (unsigned int i = 0; i < particles.size(); i++) {
        deriv[Particle::PhaseDimension*i    ] = particles[i]->vel[0];
        deriv[Particle::PhaseDimension*i + 1] = particles[i]->vel[1];
        deriv[Particle::PhaseDimension*i + 2] = particles[i]->vel[2];
        deriv[Particle::PhaseDimension*i + 3] = particles[i]->force[0]/particles[i]->mass;
        deriv[Particle::PhaseDimension*i + 4] = particles[i]->force[1]/particles[i]->mass;
        deriv[Particle::PhaseDimension*i + 5] = particles[i]->force[2]/particles[i]->mass;
    }
    return deriv;
}

Vecd ParticleSystem::getSecondDerivative() const {
    Vecd deriv(this->getStateSize());
    for (unsigned int i = 0; i < particles.size(); i++) {
        deriv[Particle::PhaseDimension*i + 0] = particles[i]->force[0]/particles[i]->mass;
        deriv[Particle::PhaseDimension*i + 1] = particles[i]->force[1]/particles[i]->mass;
        deriv[Particle::PhaseDimension*i + 2] = particles[i]->force[2]/particles[i]->mass;
        deriv[Particle::PhaseDimension*i + 3] = 0;
        deriv[Particle::PhaseDimension*i + 4] = 0;
        deriv[Particle::PhaseDimension*i + 5] = 0;
    }
    return deriv;
}

void ParticleSystem::setState(const Vecd& state) {
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->pos[0]  = state[Particle::PhaseDimension*i    ];
        particles[i]->pos[1]  = state[Particle::PhaseDimension*i + 1];
        particles[i]->pos[2]  = state[Particle::PhaseDimension*i + 2];
        particles[i]->vel[0]  = state[Particle::PhaseDimension*i + 3];
        particles[i]->vel[1]  = state[Particle::PhaseDimension*i + 4];
        particles[i]->vel[2]  = state[Particle::PhaseDimension*i + 5];
    }
}

void ParticleSystem::updateForces() {
    // clear force accumulators
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->force = Vec3(0.0, 0.0, 0.0);
    }
    // apply forces
    for (unsigned int i = 0; i < forces.size(); i++) {
        forces[i]->apply();
    }
}

Vecd ParticleSystem::getPositions() const {
    Vecd res(3*this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[3*i  ] = particles[i]->pos[0];
        res[3*i+1] = particles[i]->pos[1];
        res[3*i+2] = particles[i]->pos[2];
    }
    return res;
}

Vecd ParticleSystem::getVelocities() const {
    Vecd res(3*this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[3*i  ] = particles[i]->vel[0];
        res[3*i+1] = particles[i]->vel[1];
        res[3*i+2] = particles[i]->vel[2];
    }
    return res;
}

Vecd ParticleSystem::getAccelerations() const {
    Vecd res(3*this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[3*i  ] = particles[i]->force[0]/particles[i]->mass;
        res[3*i+1] = particles[i]->force[1]/particles[i]->mass;
        res[3*i+2] = particles[i]->force[2]/particles[i]->mass;
    }
    return res;
}

Vecd ParticleSystem::getPreviousPositions() const {
    Vecd res(3*this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[3*i  ] = particles[i]->prevPos[0];
        res[3*i+1] = particles[i]->prevPos[1];
        res[3*i+2] = particles[i]->prevPos[2];
    }
    return res;
}

Vecd ParticleSystem::getMasses() const
{
    Vecd res(this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[i] = particles[i]->mass;
    }
    return res;
}

Vecd ParticleSystem::getForces() const
{
    Vecd res(3*this->getNumParticles());
    for (unsigned int i = 0; i < particles.size(); i++) {
        res[3*i  ] = particles[i]->force[0];
        res[3*i+1] = particles[i]->force[1];
        res[3*i+2] = particles[i]->force[2];
    }
    return res;
}

void ParticleSystem::setPositions(const Vecd& pos) {
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->pos[0] = pos[3*i    ];
        particles[i]->pos[1] = pos[3*i + 1];
        particles[i]->pos[2] = pos[3*i + 2];
    }
}

void ParticleSystem::setVelocities(const Vecd& vel) {
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->vel[0] = vel[3*i    ];
        particles[i]->vel[1] = vel[3*i + 1];
        particles[i]->vel[2] = vel[3*i + 2];
    }
}

void ParticleSystem::setPreviousPositions(const Vecd& ppos) {
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->prevPos[0] = ppos[3*i    ];
        particles[i]->prevPos[1] = ppos[3*i + 1];
        particles[i]->prevPos[2] = ppos[3*i + 2];
    }
}

Vec3i ParticleSystem::intCoord(Vec3 pos) const {
    return (pos / h).array().floor().cast<int>();
}

int ParticleSystem::hashPos(Vec3 pos) const
{
    Vec3i gridPos = intCoord(pos);
    return hashCoords(gridPos(0), gridPos(1), gridPos(2));
}

void ParticleSystem::buildSpatialHashing(double h) {
    this->h = h;
    int N = getNumParticles();
    tableSize = 2 * N;
    hashTable.assign(tableSize + 1, 0); // Guard
    particleIndices.assign(N, 0);

    for (int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        int hash = hashPos(p->pos);
        hashTable[hash]++; // Increase count
    }

    // Partial sums
    for (int i = 1; i < hashTable.size(); ++i) {
        hashTable[i] = hashTable[i-1] + hashTable[i];
    }

    // Fill object ids
    for (int i = 0; i < N; ++i) {
        Particle* p = particles[i];
        int hash = hashPos(p->pos);
        hashTable[hash]--;
        particleIndices[hashTable[hash]] = i;
    }
}

std::vector<int> ParticleSystem::getNeighbors(int pi) const
{
    return getNeighbors(particles[pi]);
}

std::vector<int> ParticleSystem::getNeighbors(Particle *p) const
{
    std::vector<int> nb;

    Vec3i gridPos = intCoord(p->pos);

    std::vector<bool> visited(tableSize, false);

    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int hash = hashCoords( gridPos(0) + i, gridPos(1) + j, gridPos(2) + k );
                
                if (visited[hash]) continue;

                visited[hash] = true;
                // Cell
                int start = hashTable[hash];
                int count = hashTable[hash + 1] - start;

                for (int pidx = start; pidx < start + count; ++pidx) {
                    int pIndex = particleIndices[pidx];
                    Particle* p2 = particles[pIndex];

                    if (p != p2 and (p2->pos - p->pos).squaredNorm() <= h * h ) {
                        nb.push_back(pIndex);
                    }
                }

            }
        }
    }

    return nb;
}
