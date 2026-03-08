# MPFEM Implementation Checklist

## Naming Registry

- `LogLevel` (enum class): source `src/logger.hpp`
- `Logger` (class): source `src/logger.hpp`
- `IdRangeParser` (class): source `src/id_range_parser.hpp`
- `CaseXmlSummary` (struct): source `src/case_xml_summary_reader.hpp`
- `CaseXmlSummaryReader` (class): source `src/case_xml_summary_reader.hpp`

## Phase 0: Build Skeleton

- [x] Create `src/` and add static library target `mpfem`.
- [x] Wire `example/` to link against `mpfem`.
- [x] Enable CTest and add `tests/` subdirectory.

## Phase 1: Parsing Baseline (TDD)

- [x] Add id-range parser for domain and boundary expressions.
- [x] Add lightweight busbar case summary reader.
- [x] Add parser tests for range logic and busbar counts.
- [ ] Extend reader to full case schema (`physics`, `source`, `coupledPhysics`).
- [ ] Add material XML reader with SI-value extraction.

## Phase 2: Physics Modules

- [ ] Electrostatics module baseline.
- [ ] Heat transfer module baseline.
- [ ] Solid mechanics module baseline.

## Phase 3: Coupling

- [ ] Segregated loop: Electrostatics -> Heat -> Mechanics.
- [ ] Joule heating source update from electrostatics field.
- [ ] Thermal expansion load update from temperature field.

## Phase 4: Output and Validation

- [ ] VTU exporter for field outputs.
- [ ] COMSOL-compatible text exporter (`x y z V T disp`).
- [ ] Comparison script with `L2`, `Linf`, max relative error metrics.
