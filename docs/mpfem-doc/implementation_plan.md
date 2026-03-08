# MPFEM Implementation Checklist

## Naming Registry

- `LogLevel` (enum class): source `src/logger.hpp`
- `Logger` (class): source `src/logger.hpp`
- `IdRangeParser` (class): source `src/id_range_parser.hpp`
- `CaseXmlSummary` (struct): source `src/case_xml_summary_reader.hpp`
- `CaseXmlSummaryReader` (class): source `src/case_xml_summary_reader.hpp`
- `CaseXmlReader` (class): source `src/io/case_xml_reader.hpp`
- `MaterialXmlReader` (class): source `src/io/material_xml_reader.hpp`
- `ComsolResultReader` (class): source `src/io/comsol_result_reader.hpp`
- `ElectrostaticsSolver` (class): source `src/physics/electrostatics_solver.hpp`
- `HeatTransferSolver` (class): source `src/physics/heat_transfer_solver.hpp`
- `SolidMechanicsSolver` (class): source `src/physics/solid_mechanics_solver.hpp`
- `SegregatedCouplingEngine` (class): source `src/coupling/segregated_coupling_engine.hpp`
- `ComsolTextExporter` (class): source `src/output/comsol_text_exporter.hpp`
- `VtuExporter` (class): source `src/output/vtu_exporter.hpp`

## Phase 0: Build Skeleton

- [x] Create `src/` and add static library target `mpfem`.
- [x] Wire `example/` to link against `mpfem`.
- [x] Enable CTest and add `tests/` subdirectory.

## Phase 1: Parsing Baseline (TDD)

- [x] Add id-range parser for domain and boundary expressions.
- [x] Add lightweight busbar case summary reader.
- [x] Add parser tests for range logic and busbar counts.
- [x] Extend reader to full case schema (`physics`, `source`, `coupledPhysics`).
- [x] Add material XML reader with SI-value extraction.

## Phase 2: Physics Modules

- [x] Electrostatics module baseline.
- [x] Heat transfer module baseline.
- [x] Solid mechanics module baseline.

## Phase 3: Coupling

- [x] Segregated loop: Electrostatics -> Heat -> Mechanics.
- [x] Joule heating source update from electrostatics field.
- [x] Thermal expansion load update from temperature field.

## Phase 4: Output and Validation

- [x] VTU exporter for field outputs.
- [x] COMSOL-compatible text exporter (`x y z V T disp`).
- [x] Comparison script with `L2`, `Linf`, max relative error metrics.

## Notes

- Baseline Phase 2/3 solvers currently use simplified point-sampled approximations and are explicitly marked with `TODO` in code where FEM assembly must replace approximations.
