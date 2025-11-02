# MFEM Architecture Analysis: Lessons and Improvements

**Document Purpose:** Analysis of MFEM's architecture to identify best practices to adopt and potential improvements for our HPC-FEM solver development.

**Date:** November 2, 2025  
**Author:** Code Analysis  
**Target Project:** hpc-fem-playground

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [MFEM Architectural Strengths](#mfem-architectural-strengths)
3. [Key Design Patterns to Learn](#key-design-patterns-to-learn)
4. [Areas for Potential Improvement](#areas-for-potential-improvement)
5. [Recommendations for Our Implementation](#recommendations-for-our-implementation)

---

## Executive Summary

MFEM is a mature, production-grade finite element library with over a decade of development by LLNL. After analyzing its source code, documentation, and architecture, the following key insights emerge:

**Strengths:**
- Highly modular and extensible architecture
- Comprehensive assembly level abstraction (Legacy, Full, Partial, Element, Matrix-Free)
- Excellent device/GPU support through abstraction layers
- Rich integrator system with pluggable physics
- Strong support for high-order elements and various element types
- Well-designed memory management with host/device awareness

**Opportunities for Improvement:**
- Some over-engineering in areas that could be simplified
- Template metaprogramming could reduce runtime overhead in critical paths
- Modern C++ features (C++17/20) could improve type safety and performance
- Physics-solver coupling could be more decoupled
- Domain-specific optimizations for specific problem classes

---

## MFEM Architectural Strengths

### 1. Assembly Level Abstraction

**What MFEM Does:**
```cpp
enum class AssemblyLevel {
   LEGACY = 0,    // Fully assembled sparse matrix (CPU)
   FULL,          // Fully assembled (device-compatible)
   ELEMENT,       // Element matrices stored
   PARTIAL,       // Data at quadrature points only
   NONE,          // Matrix-free, compute on-the-fly
};
```

**Why It's Brilliant:**
- Allows runtime selection of assembly strategy based on hardware
- Progressive optimization path: start with LEGACY, optimize to PA/EA/MF
- Memory-performance tradeoffs explicit and controllable
- Each level has appropriate Extension classes (PA, EA, FA, MF)

**Learning:**
✅ **We should adopt this pattern** but start simpler:
- Begin with FULL and PARTIAL only
- Add others as needed
- Our current architecture already has physics/solver separation, so we can inject this at the `assemble()` stage

### 2. Integrator Architecture

**What MFEM Does:**
```cpp
class BilinearFormIntegrator : public NonlinearFormIntegrator {
   // Different assembly methods for different levels
   virtual void AssemblePA(const FiniteElementSpace &fes);
   virtual void AssembleElementMatrix(const FiniteElement &el, ...);
   virtual void AssembleFaceMatrix(...);
   virtual void AddMultPA(const Vector &x, Vector &y) const;
};
```

**Key Features:**
- Single integrator class handles multiple assembly strategies
- Polymorphic: DiffusionIntegrator, MassIntegrator, etc.
- Coefficient objects inject physics parameters
- Face integrators for DG methods

**Learning:**
✅ **Excellent pattern to adopt:**
```cpp
// For our implementation:
class IntegratorInterface {
    virtual void assemble_full(...) = 0;
    virtual void assemble_partial(...) = 0;
    virtual void add_mult_PA(...) = 0;
};

class DiffusionIntegrator : public IntegratorInterface {
    Coefficient* diffusion_coeff;
    // Implementation...
};
```

### 3. Coefficient System

**What MFEM Does:**
```cpp
class Coefficient {
   virtual real_t Eval(ElementTransformation &T, 
                      const IntegrationPoint &ip) = 0;
};

// Derived types:
class ConstantCoefficient : public Coefficient { ... }
class PWConstCoefficient : public Coefficient { ... }
class FunctionCoefficient : public Coefficient { ... }
class GridFunctionCoefficient : public Coefficient { ... }
```

**Why It's Powerful:**
- Complete separation of material properties from discretization
- Can couple solutions (GridFunctionCoefficient for nonlinear/coupled problems)
- Time-dependent coefficients through `SetTime()`
- Vectorial and matrix coefficients for anisotropic materials

**Learning:**
✅ **Must have for flexible physics:**
- Our current implementation lacks this flexibility
- Joule heating coefficients (electrical conductivity depending on temperature) would benefit enormously
- Would eliminate current tight coupling in `physics_joule_heating.cpp`

### 4. Device Abstraction

**What MFEM Does:**
```cpp
struct Backend {
   enum Id {
      CPU, OMP, CUDA, HIP, RAJA_CUDA, RAJA_OMP, 
      OCCA_CUDA, CEED_CUDA, DEBUG_DEVICE, ...
   };
};

class Device {
   static Device& Get() { return device_singleton; }
   void Configure(const std::string &device);
   // Memory management integration
};
```

**Architecture:**
- Single global device singleton
- Compile-time backend selection
- Memory manager aware of device/host distinction
- Vector operations automatically use appropriate backend

**Learning:**
⚠️ **Partially relevant:**
- For HPC clusters with GPUs: essential
- For our current CPU-focused development: overkill initially
- **Recommendation:** Design memory layout for future GPU support, but don't implement backends yet
- Use `MemoryClass` concept (HOST vs DEVICE) in API even if both map to HOST initially

### 5. Restriction/Prolongation Operators

**What MFEM Does:**
```cpp
class ElementRestriction : public Operator {
   void Mult(const Vector &x, Vector &y) const override;  // L-vector → E-vector
   void MultTranspose(const Vector &x, Vector &y) const;  // E-vector → L-vector
};

class L2ElementRestriction : public ElementRestrictionOperator { ... }
class FaceRestriction : public Operator { ... }
```

**Purpose:**
- **L-vectors:** Degrees of freedom (conforming, global numbering)
- **E-vectors:** Element-wise discontinuous data (local, for assembly)
- Efficient transformation between representations
- Critical for PA/EA assembly where element-wise data is preferred

**Learning:**
✅ **Important for performance:**
- Enables efficient partial assembly
- Reduces scatter/gather overhead in GPU kernels
- Our current implementation works with global DOFs directly - inefficient for PA
- Should add when implementing PA assembly level

### 6. Finite Element Space Management

**What MFEM Does:**
```cpp
class FiniteElementSpace {
   Mesh *mesh;
   FiniteElementCollection *fec;
   int vdim;  // vector dimension
   Ordering::Type ordering;  // byNODES or byVDIM
   
   // DOF management
   virtual void GetElementDofs(int i, Array<int> &dofs) const;
   virtual void GetBdrElementDofs(int i, Array<int> &dofs) const;
   
   // Operators
   virtual const Operator *GetProlongationMatrix() const;
   virtual const Operator *GetRestrictionMatrix() const;
   const ElementRestriction *GetElementRestriction(...);
};
```

**Key Insights:**
- Separation of concerns: Mesh topology vs. DOF management vs. Element types
- Support for both scalar and vector fields (vdim)
- Two orderings: XXX...YYY...ZZZ (byNODES) vs XYZ,XYZ,XYZ (byVDIM)
- Lazy initialization of restriction operators
- Non-conforming mesh support through prolongation/restriction

**Learning:**
✅ **Our current `FiniteElementSpace` wrapper is too thin:**
- We should add ordering support (important for performance)
- Element DOF queries should be cached/optimized
- Add restriction operator support for PA assembly

### 7. Form Assembly Architecture

**What MFEM Does:**
```cpp
class BilinearForm : public Matrix {
   FiniteElementSpace *fes;
   Array<BilinearFormIntegrator*> domain_integs;
   Array<BilinearFormIntegrator*> boundary_integs;
   Array<BilinearFormIntegrator*> interior_face_integs;
   
   void Assemble();  // Based on assembly level
   void FormSystemMatrix(const Array<int> &ess_tdof_list, ...);
   void FormLinearSystem(...);
};

class LinearForm : public Vector {
   Array<LinearFormIntegrator*> domain_integs;
   Array<DeltaLFIntegrator*> domain_delta_integs;
   // ...
};
```

**Pattern:**
- Forms aggregate integrators
- Integrators don't know about global system
- Form handles:
  - Integration over all elements/faces
  - Assembly into global structures
  - Boundary condition application
  - Static condensation (optional)

**Learning:**
✅ **Better than our current approach:**
```cpp
// Current (in physics classes):
void assemble(...) {
   // Everything mixed together
}

// Better (MFEM-inspired):
class BilinearForm {
   void add_domain_integrator(Integrator* integ);
   void add_boundary_integrator(Integrator* integ);
   void assemble(AssemblyLevel level);
};
```

### 8. QuadratureInterpolator

**What MFEM Does:**
```cpp
class QuadratureInterpolator {
   enum EvalFlags {
      VALUES = 1 << 0,
      DERIVATIVES = 1 << 1,
      DETERMINANTS = 1 << 2,
      PHYSICAL_DERIVATIVES = 1 << 3,
      PHYSICAL_VALUES = 1 << 4,
   };
   
   void Mult(const Vector &e_vec, unsigned eval_flags,
            Vector &q_val, Vector &q_der, Vector &q_det) const;
};
```

**Purpose:**
- Evaluate FE functions at quadrature points efficiently
- Batch evaluation across all elements
- Optimized tensor-product evaluation for hex/quad elements
- Critical component for PA assembly

**Learning:**
✅ **Essential for PA implementation:**
- We don't have this yet
- Required for efficient PA assembly
- Tensor-product optimization gives 10x+ speedup for high-order elements
- Should be added before PA assembly

### 9. Mesh and Element Transformation

**What MFEM Does:**
```cpp
class Mesh {
   Array<Element*> elements;
   Array<Vertex> vertices;
   Array<Element*> boundary;
   Array<Element*> faces;
   
   // Connectivity tables
   Table *el_to_edge;
   Table *el_to_face;
   Table *face_edge;
   
   // Geometric data
   mutable GeometricFactors *geom_factors;
};

class ElementTransformation {
   void SetIntPoint(const IntegrationPoint *ip);
   void Transform(const IntegrationPoint &ip, Vector &x);
   const DenseMatrix &Jacobian();
   real_t Weight();
};
```

**Key Design:**
- Lazy computation of geometric quantities
- Caching of Jacobians at integration points
- Support for curved elements through `Nodes` GridFunction
- Separate classes for different element types

**Learning:**
⚠️ **We rely on MFEM's mesh - good decision:**
- Don't reinvent mesh management
- Focus on higher-level abstractions
- Can wrap MFEM mesh with domain-specific convenience functions

### 10. Extension Pattern for Assembly Levels

**What MFEM Does:**
```cpp
class BilinearFormExtension : public Operator {
   virtual void Assemble() = 0;
   virtual void FormSystemMatrix(...) = 0;
   virtual void Mult(const Vector &x, Vector &y) const = 0;
};

class PABilinearFormExtension : public BilinearFormExtension { ... }
class EABilinearFormExtension : public PABilinearFormExtension { ... }
class FABilinearFormExtension : public EABilinearFormExtension { ... }
class MFBilinearFormExtension : public BilinearFormExtension { ... }
```

**Pattern:**
- Strategy pattern for assembly
- BilinearForm delegates to appropriate extension
- Each extension optimized for its assembly level
- Extensions created on-demand based on `SetAssemblyLevel()`

**Learning:**
✅ **Clean separation of concerns:**
- We should adopt similar pattern
- Allows assembly-level-specific optimizations
- Better than massive switch statements

---

## Key Design Patterns to Learn

### 1. **Integrator + Coefficient Pattern**

**Pattern:**
```
Physics = Integrator (knows PDE operator) + Coefficient (knows material properties)
```

**Example:**
```cpp
// Bad (our current approach):
class ThermalPhysics {
   double thermal_conductivity;  // Hardcoded scalar!
   void assemble() { /* uses thermal_conductivity everywhere */ }
};

// Good (MFEM approach):
class DiffusionIntegrator {
   Coefficient* diffusion_coeff;  // Can be spatially varying!
   void AssembleElementMatrix(const FiniteElement &el, ...) {
      for (integration points ip) {
         double k = diffusion_coeff->Eval(trans, ip);
         // Use k in assembly
      }
   }
};

// Usage:
PWConstCoefficient k_by_material(num_materials);
k_by_material(1) = 1.0;   // Material 1: k=1.0
k_by_material(2) = 100.0; // Material 2: k=100.0
DiffusionIntegrator integ(&k_by_material);
```

**Benefits:**
- Material properties completely decoupled from discretization
- Easy to handle heterogeneous materials
- Nonlinear problems: use GridFunctionCoefficient
- Coupled problems: one coefficient depends on another solution

### 2. **Assembly Level Strategy Pattern**

**Pattern:**
```
Form + AssemblyLevel → Extension (strategy) → Assembled Operator
```

**Example:**
```cpp
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(k));

// Strategy selection:
a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
a.Assemble();  // Creates PABilinearFormExtension internally

// Usage is uniform:
a.Mult(x, y);  // Regardless of assembly level
```

**Benefits:**
- User code doesn't change based on assembly level
- Can benchmark different assembly strategies easily
- Progressive optimization path

### 3. **Operator Algebra**

**Pattern:**
```cpp
class Operator {
   virtual void Mult(const Vector &x, Vector &y) const = 0;
   virtual void MultTranspose(const Vector &x, Vector &y) const;
   virtual const Operator* GetProlongation() const;
   virtual const Operator* GetRestriction() const;
};
```

**Benefits:**
- Uniform interface for matrices, operators, preconditioners
- Composability: RAP for coarse-grid operators
- Solver algorithms work with any Operator, not just SparseMatrix

**Learning:**
✅ **We should adopt this more thoroughly:**
- Our `SolverInterface` could accept `Operator` instead of specific matrix types
- Would enable matrix-free methods
- Would enable easy preconditioning composition

### 4. **Lazy Initialization**

**Pattern:**
```cpp
class FiniteElementSpace {
   mutable Table *el_to_edge;  // mutable!
   
   const Table& GetElementToEdgeTable() const {
      if (!el_to_edge) { BuildElementToEdgeTable(); }
      return *el_to_edge;
   }
};
```

**Benefits:**
- Don't compute what you don't need
- Reduces memory footprint
- Const methods can still build cache (via mutable)

**Learning:**
✅ **Use for expensive computations:**
- Connectivity tables
- Prolongation/restriction operators
- Geometric factors

### 5. **Marker Arrays for Selective Assembly**

**Pattern:**
```cpp
Array<int> ess_bdr_marker(mesh->bdr_attributes.Max());
ess_bdr_marker = 0;
ess_bdr_marker[0] = 1;  // Mark attribute 1 as essential

a.AddBoundaryIntegrator(new BoundaryIntegrator(...), ess_bdr_marker);
```

**Benefits:**
- Selective assembly on element/boundary attributes
- No need for separate loops
- Clear and declarative

**Learning:**
✅ **Much better than our current approach:**
```cpp
// Current: hardcoded loops over specific boundaries
for (int i = 0; i < mesh->GetNBE(); i++) {
   if (mesh->GetBdrAttribute(i) == 1) { /* ... */ }
}

// Better: declarative marker arrays
```

---

## Areas for Potential Improvement

### 1. **Over-Reliance on Raw Pointers**

**Issue:**
```cpp
class BilinearForm {
   Array<BilinearFormIntegrator*> domain_integs;  // Raw pointers!
   // Who owns these? When are they deleted?
};
```

**Problems:**
- Unclear ownership semantics
- Manual memory management error-prone
- Difficult to reason about lifetimes

**Our Improvement:**
```cpp
class BilinearForm {
   std::vector<std::unique_ptr<BilinearFormIntegrator>> domain_integs;
   // Ownership clear: BilinearForm owns the integrators
};
```

**Benefits:**
- RAII: automatic cleanup
- Move semantics: efficient transfers
- Compile-time ownership guarantees

### 2. **Limited Use of Templates for Performance**

**Issue:**
MFEM uses virtual functions extensively:
```cpp
virtual void AssembleElementMatrix(...);  // Virtual call per element!
```

**Performance Cost:**
- Virtual function call overhead
- No inlining across virtual boundary
- No template specialization optimizations

**Our Improvement:**
```cpp
template<typename Integrator>
class BilinearForm {
   void Assemble() {
      for (int i = 0; i < num_elements; i++) {
         integrator.AssembleElementMatrix(...);  // Can be inlined!
      }
   }
};
```

**Tradeoffs:**
- ✅ Faster: 10-30% speedup in hot loops
- ✅ Better optimization: inlining, loop unrolling
- ❌ Code bloat: instantiation per integrator type
- ❌ Longer compile times

**Recommendation:** Use templates for innermost loops (element assembly), keep virtual for high-level interfaces

### 3. **C++11 Limitations**

**Issue:**
MFEM targets C++11 for compatibility. This limits:
- No `std::optional` for maybe-present values
- No structured bindings
- No `if constexpr` for compile-time branches
- No concepts for template constraints
- Limited constexpr

**Our Improvement:**
We're using C++17, can leverage:
```cpp
// MFEM:
const IntegrationRule* GetRule(...) { 
   if (IntRule) return IntRule;
   else return GetDefaultRule(...);
}

// C++17:
std::optional<IntegrationRule> GetRule(...) {
   return IntRule ? std::make_optional(*IntRule) : std::nullopt;
}

// Or with if constexpr for compile-time:
template<bool HasRule>
const IntegrationRule& GetRule() {
   if constexpr (HasRule) {
      return *IntRule;
   } else {
      return GetDefaultRule();
   }
}
```

### 4. **Global State (Device Singleton)**

**Issue:**
```cpp
static Device device_singleton;  // Global mutable state
```

**Problems:**
- Difficult to test with different configurations
- Not thread-safe without locks
- Hidden dependencies

**Our Improvement:**
```cpp
class ComputeContext {
   DeviceConfig device_config;
   MemoryManager mem_manager;
   // All "global" state here
};

// Pass context explicitly
void Assemble(const ComputeContext& ctx, ...);
```

**Benefits:**
- Explicit dependencies
- Testable: each test can have its own context
- Thread-safe: each thread can have context

**Caveat:** More verbose. Acceptable tradeoff for testability.

### 5. **Monolithic Headers**

**Issue:**
```cpp
#include "fem.hpp"  // Includes EVERYTHING
```

Includes:
- All finite element types
- All integrators
- All forms
- All special features (NURBS, AMR, etc.)

**Problems:**
- Long compile times
- Unnecessary dependencies
- Hard to understand what's actually needed

**Our Improvement:**
```cpp
// Modular headers:
#include "hpcfem/core/finite_element_space.hpp"
#include "hpcfem/integrators/diffusion.hpp"
#include "hpcfem/forms/bilinear_form.hpp"
// Only what's needed
```

**Benefits:**
- Faster compilation
- Clearer dependencies
- Easier to navigate

### 6. **Error Handling**

**Issue:**
```cpp
MFEM_ABORT("Error message");  // Macro that calls exit()
```

**Problems:**
- Cannot recover from errors
- Difficult to test error conditions
- Unidiomatic C++

**Our Improvement:**
```cpp
class FEMException : public std::runtime_error {
   using std::runtime_error::runtime_error;
};

if (invalid_input) {
   throw FEMException("Invalid input: ...");
}
```

**Benefits:**
- Stack unwinding cleans up resources
- Can catch and recover
- Can test exception conditions
- Idiomatic C++

### 7. **Documentation vs. Code Clarity**

**Observation:**
MFEM has excellent Doxygen documentation, but the code itself sometimes requires it:

```cpp
void Mult(const Vector &x, Vector &y) const override;
// What does this do? Need to read docs...
```

**Our Improvement:**
```cpp
void apply_operator(const Vector& input, Vector& output) const;
// Self-documenting name

// Or:
void L_vector_to_E_vector(const Vector& L_vec, Vector& E_vec) const;
// Specific about what transformation is performed
```

### 8. **Physics-Specific Optimizations Missing**

**Issue:**
MFEM is general-purpose, so misses domain-specific optimizations:

**Example: Joule Heating (our domain)**
- Electrical conductivity is temperature-dependent: σ(T)
- Thermal conductivity is temperature-independent (often)
- Electrical solve is linear if T is fixed
- Can do fixed-point iteration: solve E, update T, repeat

**MFEM Approach:**
```cpp
// Treat as fully nonlinear system
NonlinearForm a(...);
NewtonSolver newton;
newton.Solve(...);  // Full Newton, expensive Jacobian
```

**Domain-Optimized Approach:**
```cpp
// Segregated fixed-point iteration
for (iter) {
   solve_electrical(T_fixed);    // Linear solve!
   solve_thermal(joule_heating); // Linear solve!
   check_convergence();
}
```

**Speedup:** 3-10x for weakly coupled problems

**Learning:**
✅ **Keep physics-specific optimizations in physics classes:**
- Don't force everything into MFEM's nonlinear framework
- Use MFEM for building blocks, not for solver strategy

### 9. **Verbosity of Assembly Loops**

**Issue:**
MFEM's assembly code is verbose:

```cpp
for (int i = 0; i < fes->GetNE(); i++) {
   const FiniteElement &fe = *fes->GetFE(i);
   ElementTransformation *eltrans = fes->GetElementTransformation(i);
   fes->GetElementVDofs(i, vdofs);
   
   for (int j = 0; j < integrators.Size(); j++) {
      integrators[j]->AssembleElementMatrix(fe, *eltrans, elmat);
      mat.AddSubMatrix(vdofs, vdofs, elmat, skip_zeros);
   }
}
```

**Our Improvement:**
Use ranges and algorithms:
```cpp
// C++20 ranges:
auto elements = fes.elements() | views::enumerate;
for (auto [i, elem] : elements) {
   for (auto& integ : integrators) {
      auto elmat = integ.assemble_element(elem);
      mat.add_submatrix(elem.vdofs(), elmat);
   }
}
```

**Benefits:**
- More concise
- Intent clearer
- Easier to parallelize (ranges are composable)

### 10. **Limited Compile-Time Optimization Opportunities**

**Issue:**
Many things known at compile time are handled at runtime:

```cpp
// Runtime:
if (element_type == QUAD) { /* ... */ }
else if (element_type == HEX) { /* ... */ }
```

**Our Improvement:**
```cpp
// Compile-time with tag dispatch:
template<ElementType EType>
void assemble_element() {
   if constexpr (EType == QUAD) {
      // Quad-specific code, fully optimized
   } else if constexpr (EType == HEX) {
      // Hex-specific code, fully optimized
   }
}
```

**Benefits:**
- Dead code elimination
- Better vectorization
- Smaller code size (per instantiation)

**Tradeoff:** Code bloat if overused. Use judiciously.

---

## Recommendations for Our Implementation

### Phase 1: Foundation (Current)
**Status:** ✅ Mostly complete

✅ **What we have:**
- Clean physics/solver interface separation
- MFEM mesh integration
- Basic assembly working
- MPI support

❌ **What we should add:**

1. **Coefficient System**
   ```cpp
   class Coefficient {
       virtual double eval(const ElementTransformation& T,
                          const IntegrationPoint& ip) const = 0;
   };
   ```
   - Implement ConstantCoefficient, PWConstCoefficient
   - Refactor physics classes to use coefficients
   - Add GridFunctionCoefficient for coupled problems

2. **Integrator Abstraction**
   ```cpp
   class BilinearIntegrator {
       virtual void assemble_element_matrix(...) = 0;
   };
   class DiffusionIntegrator : public BilinearIntegrator { ... };
   class MassIntegrator : public BilinearIntegrator { ... };
   ```
   - Move assembly logic from physics classes to integrators
   - Physics classes become aggregators of integrators

3. **Form Classes**
   ```cpp
   class BilinearForm {
       void add_domain_integrator(BilinearIntegrator* integ);
       void add_boundary_integrator(BilinearIntegrator* integ);
       void assemble();
   };
   ```

### Phase 2: Performance (Next 3 months)

1. **Implement Partial Assembly**
   - Add AssemblyLevel enum
   - Implement PABilinearFormExtension
   - Add ElementRestriction operators
   - Add QuadratureInterpolator
   - Benchmark: expect 2-5x speedup for high-order elements

2. **Template-Based Optimization**
   - Template-ize innermost assembly loops
   - Use tag dispatch for element types
   - Profile and verify speedup

3. **Better Memory Management**
   - Replace raw pointers with smart pointers
   - Add move semantics throughout
   - Design for future GPU support (even if not implementing yet)

### Phase 3: Advanced Features (6+ months)

1. **Matrix-Free Assembly**
   - Implement MFBilinearFormExtension
   - For large-scale problems (millions of DOFs)
   - Critical for GPU efficiency

2. **GPU Support**
   - Add Device abstraction (simplified version of MFEM's)
   - Port PA kernels to CUDA
   - Memory manager with host/device awareness

3. **Advanced Solvers**
   - Multigrid with restriction/prolongation operators
   - Matrix-free preconditioners
   - Domain decomposition methods

### Code Organization Recommendations

**Directory Structure:**
```
src/hpcfem/
├── core/
│   ├── operator.hpp              # Base operator class
│   ├── finite_element_space.hpp  # Enhanced FE space
│   ├── coefficient.hpp           # Coefficient system
│   └── assembly_level.hpp        # Assembly level enum
├── integrators/
│   ├── integrator_interface.hpp
│   ├── diffusion_integrator.hpp
│   ├── mass_integrator.hpp
│   └── convection_integrator.hpp
├── forms/
│   ├── bilinear_form.hpp
│   ├── linear_form.hpp
│   └── extensions/
│       ├── pa_extension.hpp
│       ├── ea_extension.hpp
│       └── mf_extension.hpp
├── physics/
│   ├── physics_base.hpp          # Use forms + integrators
│   ├── electrostatics.hpp
│   └── thermal.hpp
└── solvers/
    ├── solver_interface.hpp
    ├── hypre_wrapper.hpp
    └── iterative_solvers.hpp
```

### Development Priorities

**Priority 1 (Critical):**
1. Coefficient system
2. Integrator abstraction
3. Refactor physics classes to use integrators

**Priority 2 (Important):**
1. BilinearForm/LinearForm classes
2. Partial assembly implementation
3. Template-based optimization of hot paths

**Priority 3 (Nice to have):**
1. Matrix-free assembly
2. Advanced preconditioners
3. GPU support

### Testing Strategy

1. **Unit Tests:**
   - Test each integrator independently
   - Test coefficient evaluation
   - Test assembly at each level produces same result

2. **Integration Tests:**
   - Compare with MFEM results for standard problems
   - Verify convergence rates
   - Performance benchmarks vs. MFEM

3. **Performance Tests:**
   - Assembly time vs. assembly level
   - Solve time vs. problem size
   - Scalability tests (strong/weak scaling)

### Documentation Strategy

**Learn from MFEM:**
- ✅ Excellent: Example-driven documentation
- ✅ Excellent: Clear description of concepts (L-vector vs E-vector)
- ❌ Weakness: Sometimes code needs documentation to understand

**Our Approach:**
1. Self-documenting code through:
   - Clear names
   - Type-safe interfaces
   - Modern C++ idioms
2. Example programs for each physics
3. Architecture documentation (this document!)
4. API documentation (Doxygen) for details

---

## Specific Code Improvements

### 1. Refactor PhysicsInterface

**Current:**
```cpp
class PhysicsInterface {
   virtual void assemble(Matrix& A, Vector& b, ...) = 0;
};
```

**Improved:**
```cpp
class PhysicsInterface {
   virtual BilinearForm& get_bilinear_form() = 0;
   virtual LinearForm& get_linear_form() = 0;
   virtual void set_time(double t) = 0;
   virtual void update_coefficients() = 0;
};
```

**Benefits:**
- Physics class doesn't control assembly strategy
- Forms handle assembly level
- Easier to test (can assemble forms independently)

### 2. Refactor JouleHeating

**Current:**
```cpp
class PhysicsJouleHeating {
   void assemble(...) {
      // Everything hardcoded here
   }
};
```

**Improved:**
```cpp
class PhysicsJouleHeating {
   // Electrical problem
   BilinearForm a_elec;
   LinearForm b_elec;
   
   // Thermal problem  
   BilinearForm a_therm;
   LinearForm b_therm;
   
   // Coefficients
   std::unique_ptr<Coefficient> electrical_conductivity;
   std::unique_ptr<Coefficient> thermal_conductivity;
   std::unique_ptr<Coefficient> joule_heating_source;
   
   PhysicsJouleHeating(...) {
      // Setup forms
      a_elec.add_domain_integrator(
         new DiffusionIntegrator(electrical_conductivity.get()));
      
      a_therm.add_domain_integrator(
         new DiffusionIntegrator(thermal_conductivity.get()));
      b_therm.add_domain_integrator(
         new DomainLFIntegrator(joule_heating_source.get()));
      
      // Boundary conditions...
   }
   
   void solve_segregated() {
      for (int iter = 0; iter < max_iter; iter++) {
         a_elec.assemble(AssemblyLevel::PARTIAL);
         solver.solve(a_elec, b_elec, x_elec);
         
         update_joule_heating(x_elec);  // Update coefficient
         
         a_therm.assemble(AssemblyLevel::PARTIAL);
         solver.solve(a_therm, b_therm, x_therm);
         
         if (check_convergence()) break;
      }
   }
};
```

### 3. Add Assembly Level Support

```cpp
enum class AssemblyLevel {
   FULL,     // Start with this
   PARTIAL,  // Add this next for performance
   // ELEMENT and MATRIX_FREE later
};

class BilinearForm {
   AssemblyLevel assembly_level = AssemblyLevel::FULL;
   
   void set_assembly_level(AssemblyLevel level) {
      assembly_level = level;
      create_extension();  // Create appropriate extension
   }
   
   void assemble() {
      extension->assemble();  // Delegate to extension
   }
   
private:
   std::unique_ptr<BilinearFormExtension> extension;
   
   void create_extension() {
      switch (assembly_level) {
         case FULL: 
            extension = std::make_unique<FAExtension>(this);
            break;
         case PARTIAL:
            extension = std::make_unique<PAExtension>(this);
            break;
      }
   }
};
```

---

## Benchmarking Plan

### 1. Assembly Benchmarks

Compare assembly times for:
- Poisson problem (-Δu = f)
- 3D hex mesh, 100k elements
- Orders p = 1, 2, 3, 4

**Expected Results:**
```
                p=1    p=2    p=3    p=4
FULL           1.0x   2.3x   5.1x   9.8x
PARTIAL        0.8x   1.1x   1.5x   2.1x
Speedup        1.25x  2.1x   3.4x   4.7x
```

### 2. Solve Benchmarks

Total solve time for:
- Same problem
- Include assembly + linear solve

**Expected Results:**
```
                p=1    p=2    p=3    p=4
Our (FULL)     1.0x   2.1x   4.5x   8.9x
Our (PARTIAL)  0.9x   1.3x   2.0x   3.1x
MFEM (PA)      0.85x  1.2x   1.9x   2.9x
Overhead       ~5%    ~8%    ~5%    ~7%
```

Goal: Match MFEM within 10% for PA assembly.

### 3. Memory Benchmarks

Memory usage:
```
                p=1    p=2    p=3    p=4
FULL (MB)      120    280    520    890
PARTIAL (MB)   45     78     145    245
Reduction      2.7x   3.6x   3.6x   3.6x
```

---

## Conclusion

### Summary of Lessons

**✅ Adopt from MFEM:**
1. Assembly level abstraction (start simple, add complexity)
2. Integrator + Coefficient pattern (essential for flexibility)
3. Form classes as integrator aggregators
4. Operator algebra (uniform interface)
5. Restriction/prolongation operators (for PA)
6. Lazy initialization pattern
7. Marker arrays for selective assembly
8. Clear separation: mesh/FE space/DOF management

**⚠️ Improve Upon MFEM:**
1. Modern C++ (C++17): smart pointers, optionals, if constexpr
2. Templates for performance (innermost loops)
3. Modular headers (faster compilation)
4. Exception-based error handling
5. Less global state (testability)
6. Domain-specific optimizations
7. Cleaner code through better naming

**❌ Avoid from MFEM:**
1. Monolithic headers
2. Heavy use of raw pointers
3. Global singletons
4. Macro-based error handling
5. Over-generalization where not needed

### Next Steps

1. **Week 1-2:** Implement Coefficient system
2. **Week 3-4:** Implement Integrator interfaces
3. **Week 5-6:** Implement BilinearForm/LinearForm classes
4. **Week 7-8:** Refactor physics classes to use new architecture
5. **Week 9-12:** Implement Partial Assembly
6. **Month 4-6:** Optimize and benchmark

### Success Criteria

1. **Correctness:** Match MFEM results for standard problems (L2 error < 1e-10)
2. **Performance:** Within 10% of MFEM for PA assembly
3. **Flexibility:** Easy to add new physics (< 200 lines for new PDE)
4. **Maintainability:** New contributor can understand codebase in < 1 week
5. **Scalability:** Good MPI scaling (> 80% efficiency up to 100 ranks)

---

**End of Document**

This analysis provides a comprehensive roadmap for our FEM solver development, learning from MFEM's 15+ years of development while avoiding its pitfalls and leveraging modern C++ features for a cleaner, more maintainable, and potentially faster implementation.
