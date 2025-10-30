# AI Assistant Instructions for AE543 Project

## Project Overview
This project focuses on structural dynamics analysis of aircraft components, specifically the T-33 wing vibrations. The project combines theoretical analysis with computational modeling using both Jupyter notebooks and finite element analysis (FEM).

## Repository Structure
- `Project-1.ipynb`: Main Jupyter notebook containing the analysis workflow
- `FEM/`: Directory containing finite element analysis files
  - `CantileverBeam.bdf`: Nastran input file for structural analysis
  - `cantileverbeam.f04`, `cantileverbeam.f06`: Nastran output files
  - `cantileverbeam.h5`: HDF5 results file
  - `patran.ses.01`: Patran session file
- `Images/`: Supporting images for documentation

## Key Parameters
When working with this codebase, be aware of these critical parameters:
- Wing length: 5.0 m
- Tip mass (mp): 100 kg
- Chord: 0.70 m
- Height: 0.25 m
- Wall thickness: 0.010 m
- Material: Aluminum 6061-T6
  - E = 69 GPa
  - ρ = 2700 kg/m²
  - ν = 0.33

## Analysis Workflows
1. **Jupyter Notebook Analysis**
   - Mathematical modeling of wing structure
   - Both SDoF and MDoF system analysis
   - Integration of Lagrangian and Newtonian mechanics
   - Visualization requirements: All plots must include labeled axes, grids, titles, and legends

2. **FEM Analysis**
   - Uses MSC.Nastran/Patran for structural analysis
   - Primary analysis types:
     - Linear static analysis (SOL 101)
     - Multiple load cases including 1G down force

## Conventions and Standards
- All figures must include:
  - Properly labeled axes
  - Grid lines
  - Title
  - Legend (where applicable)
- Extra credit available for additional engineering analysis figures (up to 5 points)
- All additional visualizations must include engineering interpretation

## Integration Points
- Jupyter notebook serves as both computation and reporting platform
- FEM results should be referenced and compared with analytical solutions
- Results should demonstrate connection between solid mechanics and vibration fundamentals

Remember to maintain academic rigor in analysis and presentation, as this is an educational project in structural dynamics.