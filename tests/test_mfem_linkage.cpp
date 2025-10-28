/**
 * @file test_mfem_linkage.cpp
 * @brief Test that MFEM headers can be included and linked successfully
 * 
 * This is a basic smoke test to verify the build system is correctly
 * configured to find and link against MFEM.
 * 
 * Phase: 0, Step: 0.1
 */

#include "mfem.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
    // Initialize MPI if available
#ifdef MFEM_USE_MPI
    mfem::MPI_Session mpi(argc, argv);
    if (mpi.Root())
    {
        std::cout << "MFEM with MPI support detected" << std::endl;
    }
#else
    std::cout << "MFEM without MPI support" << std::endl;
#endif

    // Try to create a simple MFEM object to verify linkage
    mfem::Vector test_vector(3);
    test_vector = 1.0;
    
    double norm = test_vector.Norml2();
    constexpr double EXPECTED_NORM = 1.732050807568877; // sqrt(3)
    constexpr double TOLERANCE = 1.0e-12;
    
    if (std::abs(norm - EXPECTED_NORM) < TOLERANCE)
    {
#ifdef MFEM_USE_MPI
        if (mpi.Root())
#endif
        {
            std::cout << "Test passed: MFEM linkage successful" << std::endl;
        }
        return EXIT_SUCCESS;
    }
    else
    {
#ifdef MFEM_USE_MPI
        if (mpi.Root())
#endif
        {
            std::cerr << "Test failed: Expected norm " << EXPECTED_NORM 
                      << " but got " << norm << std::endl;
        }
        return EXIT_FAILURE;
    }
}
