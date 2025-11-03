# VTK-LPT-Framework
A solver-agnostic Lagrangian particle tracking framework for ParaView

## Overview
VTK-LPT-Framework extends ParaView/VTK with a solver-agnostic Lagrangian particle tracking workflow. It provides a drop-in replacement for the stock `vtkLagrangianParticleTracker` and adds the `vtkLagrangianFinallyModel` integration model to capture gravity, buoyancy, drag (including Cunningham slip), shear and rotation lift, pressure-gradient forces, added mass, Brownian motion, and turbulent dispersion within a single plugin-friendly package.

## Highlights
- Rich force model: toggle gravity, buoyancy, drag, Cunningham slip, shear lift, rotation lift, pressure gradient, added mass, Brownian diffusion, turbulent dispersion, and radius-aware wall interactions individually.
- Multi-file and multi-loop support: sequentially load multiple flow datasets (`FileListNumber` + `FileName`) and reuse seeds for repeated integration loops (`IfUseMultiLoop`, `MultiLoopNumber`).
- Enhanced particle state: records instantaneous flow velocity, relative velocity, and angular velocity components along the path for advanced diagnostics.
- Configurable physics: override viscosity, temperature, normal/tangential restitution, and integration step time via ParaView properties.
- ParaView-native integration: retains the familiar `LagrangianParticleTracker` UI (seed and surface helpers, property panels) while exposing additional properties and outputs (`ParticlePaths`, `ParticleInteractions`).
- BSD 3-Clause license for permissive reuse.

## Repository Layout
```
VTK-LPT-Framework/
  src/                          # Custom C++ implementation of tracker and integration model
    vtkLagrangianFinallyModel.{h,cxx}
    vtkLagrangianParticleTracker.{h,cxx}
  plugin/                       # ParaView plugin XML and module declaration
    CMakeLists.txt               # vtk_module definition consumed by paraview_add_plugin
    LagrangianParticleTracker   # ServerManager XML (GUI and property definitions)
  LICENSE                       # BSD 3-Clause license
  README.md
```

## Requirements
- ParaView 5.10 or newer (ParaView 5.11+/VTK 9.2+ recommended for `vtk_module_add_module` support).
- CMake 3.21 or newer.
- C++14 capable compiler (the standard ParaView toolchains are sufficient).
- Flow, seed, and surface datasets that expose the required scalar/vector fields.

## Build and Installation

The repository does not ship a top-level `CMakeLists.txt`. Choose one of the integration options below.

### Option 1: build as an external ParaView plugin (recommended)
1. Clone the repository:
   ```bash
   git clone https://github.com/yimengduhhh/VTK-LPT-Framework.git
   ```
2. Create a driver `CMakeLists.txt` next to the clone (or in a parent directory):
   ```cmake
   cmake_minimum_required(VERSION 3.21)
   project(VTKLPTFrameworkPlugin LANGUAGES CXX)

   find_package(ParaView REQUIRED COMPONENTS vtkRemotingServerManager)
   include(${ParaView_USE_FILE})

   paraview_plugin_scan(
     PLUGIN_FILES    "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/plugin/LagrangianParticleTracker"
     PLUGINS_FOUND   lpt_plugins
     AUTOLOAD        lpt_autoload
   )

   paraview_add_plugin(VTKLPTFramework
     VERSION                "0.1"
     MODULES                LagrangianParticleTracker::vtkLagrangianParticleTracker
     MODULE_FILES           "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/plugin/CMakeList.txt"
     SERVER_MANAGER_XML     "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/plugin/LagrangianParticleTracker"
     SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/src/vtkLagrangianFinallyModel.cxx"
       "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/src/vtkLagrangianParticleTracker.cxx"
     SERVER_MANAGER_SOURCES
       "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/src/vtkLagrangianFinallyModel.h"
       "${CMAKE_CURRENT_SOURCE_DIR}/VTK-LPT-Framework/src/vtkLagrangianParticleTracker.h"
   )
   ```
   ParaView 5.12 introduced the `paraview_plugin_add` macro. Replace the final block if you target the newer API.
3. Configure and build:
   ```bash
   cmake -S . -B build -DParaView_DIR=/path/to/paraview-build
   cmake --build build --config Release
   ```
4. Deploy the resulting plugin library (under `build/plugins/VTKLPTFramework`) by adding it to `PARAVIEW_PLUGIN_PATH` or loading it through *Tools -> Manage Plugins* in ParaView.

### Option 2: patch ParaView/VTK sources directly
1. Copy the contents of `src/` into the ParaView or VTK source tree under `VTK/Filters/FlowPaths/`, replacing the stock implementations.
2. Replace the stock plugin files in `Plugins/LagrangianParticleTracker/` with the versions in `plugin/`.
3. Reconfigure and rebuild ParaView.
4. Verify downstream filters that rely on `vtkLagrangianParticleTracker`, because this approach updates the core library.

## Required Data Fields

### Flow input (connected to `Flow Input`)
| Default name         | Components | Location     | Required | Purpose                                              |
|----------------------|------------|--------------|----------|------------------------------------------------------|
| FlowVelocity         | 3          | cell or point| yes      | Local fluid velocity                                 |
| FlowDensity          | 1          | cell data    | yes      | Fluid density                                        |
| FlowDynamicViscosity | 1          | cell data    | yes*     | Dynamic viscosity (*may be overridden manually)      |
| dxVelocitydy         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| dxVelocitydz         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| dyVelocitydx         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| dyVelocitydz         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| dzVelocitydx         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| dzVelocitydy         | 1          | cell data    | yes      | Velocity gradient tensor component                   |
| Temperature          | 1          | cell data    | yes*     | Fluid temperature (*may be overridden manually)      |
| dpressuredx          | 1          | cell data    | yes**    | Pressure gradient x-component (**if pressure forces) |
| dpressuredy          | 1          | cell data    | yes**    | Pressure gradient y-component                        |
| dpressuredz          | 1          | cell data    | yes**    | Pressure gradient z-component                        |
| turbulentKE          | 1          | cell data    | yes***   | Turbulent kinetic energy (***if turbulent dispersion)|
| turbulentDR          | 1          | cell data    | yes***   | Turbulent dissipation rate                           |
| fillId               | 1          | cell data    | optional | Region or phase identifier                           |
| RMSE_X               | 1          | cell data    | yes***   | Velocity fluctuation RMS (x)                         |
| RMSE_Y               | 1          | cell data    | yes***   | Velocity fluctuation RMS (y)                         |
| RMSE_Z               | 1          | cell data    | yes***   | Velocity fluctuation RMS (z)                         |
| Vorticity            | 3          | cell data    | yes      | Vorticity vector (recomputed internally but required)|

### Seed input (connected to `Particle Seeds`)
| Default name     | Components | Location   | Required | Purpose                             |
|------------------|------------|------------|----------|-------------------------------------|
| ParticleDiameter | 1          | point data | yes      | Particle diameter (mass and radius) |
| ParticleDensity  | 1          | point data | yes      | Particle density                     |
| TrackTime        | 1          | point data | yes      | Accumulated tracking time budget     |
| Initial velocity | 3          | point data | optional | Overrides the default initial speed  |

Any 3-component point-data array can be mapped as the initial velocity through the property panel.

### Surface input (optional, connected to `Surface`)
- Supplies wall or obstacle geometry for particle interactions.
- Use together with `IfConsiderRadius`, `NormalRemaining`, and `TangentialRemaining` to control restitution during impacts.

## Usage Workflow in ParaView
1. Load the compiled plugin via *Tools -> Manage Plugins* and enable auto-load if desired.
2. Add the flow dataset to the pipeline and connect it to `Flow Input`.
3. Provide a seed dataset (point cloud, volumetric generator, etc.) and connect it to `Particle Seeds`, ensuring the required point-data arrays exist.
4. Optionally connect a surface dataset to `Surface`.
5. Choose **Finally Model** in the *Integration Model* drop-down and map each required array. The default Matida model remains available if you prefer stock behaviour.
6. Configure integration parameters:
   - Set `OneStepTime`, `NumberOfSteps`, and `MaximumIntegrationTime`.
   - Define multi-file playback with `FileListNumber` and `FileName` (the file name should point to the first file; subsequent files must share the prefix and use incrementing suffixes).
   - Enable or disable individual force contributions (`IfOpenGravity`, `IfOpenDrag`, `IfOpenTurbulentDispersion`, etc.).
   - Turn on `IfManualViscosity` or `IfManualTemperature` when the flow dataset does not provide the corresponding fields.
7. Click *Apply* to execute the integration. Inspect the outputs:
   - `ParticlePaths`: `vtkPolyData` with polylines and per-point attributes (position, flow velocity, relative velocity, angular velocity, integration time).
   - `ParticleInteractions`: composite data capturing events such as surface hits.
8. Post-process with standard ParaView filters (for example `Plot Over Line`, `Temporal Statistics`, `Integrate Variables`, or `Append Attributes` on the interaction output).

## Tips and Troubleshooting
- Array mismatches: run-time errors usually indicate a missing or mis-mapped array. Revisit the property panel and confirm attribute locations.
- Multi-file ingestion: file names must follow a consistent prefix plus incrementing numeric suffix (for example `case0.vtu`, `case1.vtu`). Set `FileListNumber` to the total number of files.
- Premature termination: ensure `TrackTime` values are consistent with `OneStepTime`, and that `MaximumNumberOfSteps` is large enough for the intended path length.
- Turbulent dispersion stability: when stochastic fluctuations dominate, reduce `OneStepTime`, verify RMS magnitudes, or disable `IfOpenTurbulentDispersion`.
- Wall penetration: enable `IfConsiderRadius`, check surface normals, and adjust `NormalRemaining` and `TangentialRemaining`.

## Contributing
Issues and pull requests are welcome. When reporting a problem, include:
- A short description of the dataset (mesh type, units, relevant value ranges).
- Which force toggles were enabled.
- ParaView/VTK version, operating system, and compiler details.

## License
Released under the [BSD 3-Clause License](LICENSE). Suggested citation: Yimeng Du, Yan Cui, and VTK-LPT-Framework contributors (2025).

## Contact
- Issues: https://github.com/yimengduhhh/VTK-LPT-Framework/issues
- For academic collaboration or feature requests, contact the maintainers via the repository information.
