# CA-Simulations

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

```bash
qmake6 Simulations.pro & make -j & ./Simulations
```

Or open `CA-Simulations.sln` in Visual Studio.
