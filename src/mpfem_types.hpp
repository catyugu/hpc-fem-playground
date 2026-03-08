#ifndef MPFEM_MPFEM_TYPES_HPP
#define MPFEM_MPFEM_TYPES_HPP

#include "mfem.hpp"
#include "model/case_model.hpp"

#ifdef MFEM_USE_MPI
// Parallel types
using FemMesh = mfem::ParMesh;
using FemFECollection = mfem::H1_FECollection;
using FemFEspace = mfem::ParFiniteElementSpace;
using FemGridFunction = mfem::ParGridFunction;
using FemBilinearForm = mfem::ParBilinearForm;
using FemLinearForm = mfem::ParLinearForm;
using FemMatrix = mfem::HypreParMatrix;
using FemVector = mfem::Vector;
#else
// Serial types
using FemMesh = mfem::Mesh;
using FemFECollection = mfem::H1_FECollection;
using FemFEspace = mfem::FiniteElementSpace;
using FemGridFunction = mfem::GridFunction;
using FemBilinearForm = mfem::BilinearForm;
using FemLinearForm = mfem::LinearForm;
using FemMatrix = mfem::SparseMatrix;
using FemVector = mfem::Vector;
#endif

#endif