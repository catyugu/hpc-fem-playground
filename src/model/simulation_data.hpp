#ifndef MPFEM_SIMULATION_DATA_HPP
#define MPFEM_SIMULATION_DATA_HPP

#include <vector>

namespace mpfem {

/**
 * @brief 3D coordinate for nodal output data.
 */
struct Coordinate3D {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @brief Single row from COMSOL result file.
 */
struct ComsolResultRow {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double potential = 0.0;
    double temperature = 0.0;
    double displacement = 0.0;
};

/**
 * @brief Coupled baseline simulation result fields.
 */
struct CoupledFieldResult {
    std::vector<Coordinate3D> coordinates;
    std::vector<double> electricPotential;
    std::vector<double> temperature;
    std::vector<double> displacement;
    std::vector<double> jouleSource;
};

} // namespace mpfem

#endif
