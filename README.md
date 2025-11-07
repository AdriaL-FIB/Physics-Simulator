# Physics Simulator

Interactive particle-based physics simulations written in C++ with Qt and OpenGL.

It includes several numerical integrators, force models (gravity, springs, drag, wind, SPH) and a set of scenes to compare them visually.

## Scenes

- **Projectiles**: projectile motion compared against the analytic solution.
- **N-Body**: gravitational interaction between particles.
- **Fountain**: particle fountain.
- **Rope**: chain of particles connected by springs.
- **Cloth**: cloth simulation with interactive particle pinning.
- **Fluid Simulation**: SPH fluid simulation.
- **Test Integrators**: comparison of the integrators.
- **Test Colliders**: collision tests.

## Integrators

Euler, Symplectic Euler, Implicit Euler (Taylor and Newton variants), Midpoint, RK2, RK4 and Verlet.

## Controls

### General Controls

- **Move Colliders**: `Ctrl + Right Click` while hovering over the collider

### Cloth Scene

- **Select Particle**: `Right Click` on a particle
- **Fix/Unfix Particle**: Press `F` while a particle is selected

### Fluid Simulation Scene

- **Move Fan**: `Shift + Right Click`
- **Move Bounds**: `Alt + Right Click`

## Important Notes

- **Implicit Integrators**: Implicit integration methods (Euler-Taylor, Euler-Newton) work best with a small number of particles. Performance degrades significantly with larger particle systems.

## Building

Requirements: Qt 6 and a C++ compiler with OpenGL support.

With qmake:

```bash
qmake6 Simulations.pro
make -j
./Simulations
```

On Windows you can also open `CA-Simulations.sln` in Visual Studio.
