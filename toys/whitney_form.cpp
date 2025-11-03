#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <numeric>     // For std::iota
#include <stdexcept>
#include <algorithm>
#include <iomanip>     // For formatting output

// --- Type Definitions for Clarity ---
using MultiIndex = std::vector<int>;
using BasisKForm = std::vector<int>;
using Point = std::vector<double>;

// --- Point/Vector Math Helpers ---
Point operator-(const Point& a, const Point& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}
double dot(const Point& a, const Point& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double norm(const Point& a) {
    return std::sqrt(dot(a, a));
}

// Forward declarations
class Polynomial;
class DifferentialForm;

/**
 * @brief Represents a k-dimensional subsimplex (vertex, edge, face, or tet).
 */
struct Subsimplex {
    std::vector<int> vertex_indices;
    int dim;
    double volume; // k-dimensional volume (1 for vertex, length for edge, area for face)

    // Gets the BasisKForm associated with this subsimplex relative to v0
    // e.g., face {0, 1, 2} -> basis {1, 2} (dλ₁ ∧ dλ₂)
    BasisKForm getBasis() const {
        BasisKForm basis;
        for (int idx : vertex_indices) {
            if (idx != 0) {
                basis.push_back(idx);
            }
        }
        std::sort(basis.begin(), basis.end());
        return basis;
    }
};

/**
 * @brief Holds geometry for our reference tetrahedron.
 */
struct Tetrahedron {
    Point v0, v1, v2, v3;
    double volume;

    Tetrahedron(Point p0, Point p1, Point p2, Point p3)
        : v0(p0), v1(p1), v2(p2), v3(p3) {
        // Volume = 1/6 * |(v1-v0) . ((v2-v0) x (v3-v0))|
        volume = std::abs(dot(v1 - v0, cross(v2 - v0, v3 - v0))) / 6.0;
    }

    // Get barycentric coords for a point p
    std::vector<double> getBarycentricCoords(const Point& p) const {
        // This is simplified for our specific reference tetrahedron
        double x = p[0], y = p[1], z = p[2];
        return {1 - x - y - z, x, y, z}; // {λ0, λ1, λ2, λ3}
    }

    // Get an edge as a Subsimplex
    Subsimplex getEdge(int i, int j) const {
        std::vector<Point> v = {v0, v1, v2, v3};
        double length = norm(v[j] - v[i]);
        return {{i, j}, 1, length};
    }

    // Get a face as a Subsimplex
    Subsimplex getFace(int i, int j, int k) const {
        std::vector<Point> v = {v0, v1, v2, v3};
        double area = 0.5 * norm(cross(v[j] - v[i], v[k] - v[i]));
        return {{i, j, k}, 2, area};
    }
};

// --- Helper Functions ---
long long factorial(int n) {
    if (n < 0) return 1; // By convention for the formula
    long long res = 1;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}
std::string miToString(const MultiIndex& mi) { /* ... same as before ... */ 
    std::stringstream ss;
    for (size_t i = 0; i < mi.size(); ++i) {
        if (mi[i] > 0) {
            ss << "λ" << i;
            if (mi[i] > 1) ss << "^" << mi[i];
        }
    }
    std::string s = ss.str();
    return s.empty() ? "1" : s;
}
std::string basisToString(const BasisKForm& basis) { /* ... same as before ... */ 
    if (basis.empty()) return "1"; // For 0-forms
    std::stringstream ss;
    for (size_t i = 0; i < basis.size(); ++i) {
        ss << "dλ" << basis[i] << (i < basis.size() - 1 ? " ∧ " : "");
    }
    return ss.str();
}

////////////////////////////////////////////////////////////////////////////////
// CLASS: Polynomial (Represents a 0-Form)
////////////////////////////////////////////////////////////////////////////////
class Polynomial {
public:
    std::map<MultiIndex, double> terms;

    Polynomial() = default;
    Polynomial(double c) {
        if (c != 0.0) terms[{0, 0, 0, 0}] = c;
    }
    Polynomial(const MultiIndex& mi, double c) {
        if (c != 0.0) terms[mi] = c;
    }

    double evaluate(const std::vector<double>& lambdas) const { /* ... same as before ... */ 
        double result = 0.0;
        for (const auto& term : terms) {
            double term_val = term.second; // The coefficient
            const MultiIndex& mi = term.first;
            for (size_t i = 0; i < lambdas.size(); ++i) {
                term_val *= std::pow(lambdas[i], mi[i]);
            }
            result += term_val;
        }
        return result;
    }

    Polynomial operator+(const Polynomial& other) const { /* ... same as before ... */ 
        Polynomial result = *this;
        for (const auto& term : other.terms) {
            result.terms[term.first] += term.second;
        }
        return result;
    }

    Polynomial operator*(const Polynomial& other) const { /* ... same as before ... */ 
        Polynomial result;
        for (const auto& term1 : terms) {
            for (const auto& term2 : other.terms) {
                MultiIndex new_mi = {
                    term1.first[0] + term2.first[0],
                    term1.first[1] + term2.first[1],
                    term1.first[2] + term2.first[2],
                    term1.first[3] + term2.first[3]
                };
                result.terms[new_mi] += term1.second * term2.second;
            }
        }
        return result;
    }

    // Integrate over the full tetrahedron
    double integrate(const Tetrahedron& tet) const {
        // Use the subsimplex integration with all vertices
        Subsimplex s_tet = {{0, 1, 2, 3}, 3, tet.volume};
        return this->integrate(s_tet);
    }

    /**
     * @brief Integrates the polynomial over a k-subsimplex.
     * This is the core function for computing degrees of freedom.
     */
    double integrate(const Subsimplex& s) const {
        double integral = 0.0;
        int n_k = s.dim; // Dimension of the subsimplex

        for (const auto& term : terms) {
            const MultiIndex& mi_n = term.first; // Full 4-index
            double coeff = term.second;

            // 1. Restrict the polynomial to the subsimplex
            // If the polynomial has a λ_j term where j is not a vertex
            // of the subsimplex, that term is zero on the subsimplex.
            bool on_simplex = true;
            int sum_i_k = 0;
            std::vector<int> mi_k; // Multi-index restricted to the subsimplex
            
            for (int i = 0; i < 4; ++i) {
                bool is_vertex = false;
                for (int v_idx : s.vertex_indices) {
                    if (i == v_idx) {
                        is_vertex = true;
                        break;
                    }
                }

                if (mi_n[i] > 0 && !is_vertex) {
                    on_simplex = false;
                    break;
                }
                if (is_vertex) {
                    mi_k.push_back(mi_n[i]);
                    sum_i_k += mi_n[i];
                }
            }
            if (!on_simplex) continue; // This term is 0 on s

            // 2. Apply the barycentric integral formula
            // ∫_s (Π λ_j^ij) dV = ( (Π ij!) * k! ) / ( (Σij) + k )! * |s|
            long long num = factorial(n_k);
            for (int i_k : mi_k) {
                num *= factorial(i_k);
            }
            long long den = factorial(sum_i_k + n_k);
            
            integral += coeff * (static_cast<double>(num) / den) * s.volume;
        }
        return integral;
    }

    std::string toString() const { /* ... same as before ... */ 
        if (terms.empty()) return "0";
        std::stringstream ss;
        bool first = true;
        for (const auto& term : terms) {
            if (term.second == 0.0) continue;
            ss << (term.second > 0 ? (first ? "" : " + ") : " - ");
            ss << std::abs(term.second) << "*" << miToString(term.first);
            first = false;
        }
        return ss.str();
    }
};

Polynomial lambda(int i) { /* ... same as before ... */ 
    if (i < 0 || i > 3) throw std::out_of_range("Invalid lambda index");
    MultiIndex mi = {0, 0, 0, 0};
    mi[i] = 1;
    return Polynomial(mi, 1.0);
}

////////////////////////////////////////////////////////////////////////////////
// CLASS: DifferentialForm (Represents a k-Form)
////////////////////////////////////////////////////////////////////////////////
class DifferentialForm {
public:
    int k;
    std::map<BasisKForm, Polynomial> components;

    DifferentialForm(int k) : k(k) {}

    DifferentialForm operator+(const DifferentialForm& other) const { /* ... same as before ... */ 
        if (k != other.k) throw std::invalid_argument("Cannot add forms of different k");
        DifferentialForm result(k);
        result.components = components;
        for (const auto& comp : other.components) {
            result.components[comp.first] = result.components[comp.first] + comp.second;
        }
        return result;
    }
    DifferentialForm operator*(const Polynomial& p) const { /* ... same as before ... */ 
        DifferentialForm result(k);
        for (const auto& comp : components) {
            result.components[comp.first] = comp.second * p;
        }
        return result;
    }

    std::map<BasisKForm, double> evaluateCoefficients(const std::vector<double>& lambdas) const { /* ... same as before ... */ 
        std::map<BasisKForm, double> evaled_coeffs;
        for (const auto& comp : components) {
            evaled_coeffs[comp.first] = comp.second.evaluate(lambdas);
        }
        return evaled_coeffs;
    }

    /**
     * @brief Integrates the k-form over a k-dimensional subsimplex.
     * This calculates ∫_s ω
     * It assumes the basis (dλ)_σ is the "volume form" for s, which
     * is true for subsimplices s attached to vertex 0.
     */
    double integrate(const Subsimplex& s) const {
        if (k != s.dim) {
            throw std::invalid_argument("Cannot integrate k-form over subsimplex of different dimension");
        }
        if (k == 0) {
            // This is a 0-form (Polynomial) on a 0-subsimplex (Vertex)
            // This is just evaluation.
            Point v = (s.vertex_indices[0] == 0) ? Point{0,0,0} : 
                      (s.vertex_indices[0] == 1) ? Point{1,0,0} :
                      (s.vertex_indices[0] == 2) ? Point{0,1,0} : Point{0,0,1};
            return components.at({}).evaluate(Tetrahedron(v,v,v,v).getBarycentricCoords(v));
        }

        // Get the basis form {σ} corresponding to this subsimplex
        // This simple mapping only works for subsimplices attached to v0
        BasisKForm basis_s = s.getBasis();
        if (basis_s.empty() && k > 0) {
            throw std::invalid_argument("Integration on this subsimplex is not supported by this example");
        }

        // Find the polynomial coefficient 'a_σ' for this basis form
        Polynomial a_s;
        if (components.count(basis_s)) {
            a_s = components.at(basis_s);
        } else {
            return 0.0; // This form has no component on this basis
        }

        // Now, compute the integral ∫_s a_s * (dλ)_σ
        // From the paper, (dλ)_σ = ±(1 / k!|s|) * vol_s
        // We assume "+" for our reference simplex.
        // ∫_s a_s * (dλ)_σ = ∫_s a_s * (1 / (k!|s|)) vol_s
        //                   = (1 / (k!|s|)) * ∫_s a_s vol_s
        
        // a_s.integrate(s) computes ∫_s a_s vol_s
        double integral_a_s = a_s.integrate(s);
        
        return (1.0 / (factorial(k) * s.volume)) * integral_a_s;
    }

    std::string toString() const { /* ... same as before ... */ 
        if (components.empty()) return "0";
        std::stringstream ss;
        bool first = true;
        for (const auto& comp : components) {
            std::string p_str = comp.second.toString();
            if (p_str == "0") continue;
            
            ss << (first ? "" : " + ");
            ss << "(" << p_str << ") * [" << basisToString(comp.first) << "]";
            first = false;
        }
        return ss.str();
    }
};

DifferentialForm operator*(const Polynomial& p, const DifferentialForm& form) { /* ... same as before ... */ 
    return form * p;
}

std::pair<BasisKForm, int> wedgeBasis(const BasisKForm& a, const BasisKForm& b) { /* ... same as before ... */ 
    BasisKForm combined = a;
    combined.insert(combined.end(), b.begin(), b.end());

    int n = combined.size();
    if (n > 3) return {{}, 0}; // Can't have > 3-form in R³

    // Check for duplicates
    std::vector<int> sorted_check = combined;
    std::sort(sorted_check.begin(), sorted_check.end());
    for (size_t i = 0; i < sorted_check.size() - 1; ++i) {
        if (sorted_check[i] == sorted_check[i+1]) {
            return {{}, 0}; // e.g., dλ₁ ∧ dλ₁ = 0
        }
    }

    // Calculate sign from number of inversions to sort
    int sign = 1;
    BasisKForm sorted_form = combined;
    int swaps = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (sorted_form[j] > sorted_form[j+1]) {
                std::swap(sorted_form[j], sorted_form[j+1]);
                swaps++;
            }
        }
    }
    
    if (swaps % 2 != 0) sign = -1;

    return {sorted_form, sign};
}

DifferentialForm wedge(const DifferentialForm& formA, const DifferentialForm& formB) { /* ... same as before ... */ 
    int new_k = formA.k + formB.k;
    if (new_k > 3) throw std::out_of_range("Wedge product k > 3");

    DifferentialForm result(new_k);

    for (const auto& compA : formA.components) {
        for (const auto& compB : formB.components) {
            // 1. Compute new polynomial coefficient
            Polynomial new_poly = compA.second * compB.second;

            // 2. Compute new basis form and sign
            auto wedge_result = wedgeBasis(compA.first, compB.first);
            BasisKForm new_basis = wedge_result.first;
            int sign = wedge_result.second;

            if (sign == 0) continue; // This term is zero

            // 3. Add to the result, scaling by the sign
            DifferentialForm term_to_add(new_k);
            term_to_add.components[new_basis] = new_poly;
            result = result + (term_to_add * Polynomial(sign));
        }
    }
    return result;
}

// d_lambda(i) now correctly handles dλ₀
DifferentialForm d_lambda(int i) {
    if (i < 0 || i > 3) throw std::out_of_range("Invalid dλ index.");
    
    if (i == 0) {
        // dλ₀ = -dλ₁ - dλ₂ - dλ₃
        DifferentialForm d_l1 = d_lambda(1);
        DifferentialForm d_l2 = d_lambda(2);
        DifferentialForm d_l3 = d_lambda(3);
        return (d_l1 + d_l2 + d_l3) * Polynomial(-1.0);
    } else {
        // dλi is an independent basis vector
        DifferentialForm dli(1);
        dli.components[{i}] = Polynomial(1.0);
        return dli;
    }
}


////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION (Example Usage)
////////////////////////////////////////////////////////////////////////////////

int main() {
    // --- Setup ---
    Tetrahedron ref_tet({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    Subsimplex e_01 = ref_tet.getEdge(0, 1);
    Subsimplex f_012 = ref_tet.getFace(0, 1, 2);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "--- Reference Tetrahedron Volume: " << ref_tet.volume << std::endl;
    std::cout << "--- Edge e_01 Length: " << e_01.volume << std::endl;
    std::cout << "--- Face f_012 Area: " << f_012.volume << std::endl;

    // --- 1. HIGHER ORDER POLYNOMIALS (0-FORMS) ---
    std::cout << "\n--- 1. Higher-Order 0-Forms (r=3) ---" << std::endl;
    // P_rΛ⁰(T)
    // Construct p = 5 * λ₀ * λ₁ * λ₂
    Polynomial p_r3 = lambda(0) * lambda(1) * lambda(2) * Polynomial(5.0);
    
    std::cout << "Constructed Poly p = " << p_r3.toString() << std::endl;
    std::cout << "  (Degree is 1+1+1 = 3)" << std::endl;
    std::cout << "  Integrate ∫_T p dV: " << p_r3.integrate(ref_tet) << std::endl;

    // --- 2. INTEGRATION OVER SUBSIMPLICES ---
    
    // --- 2a. Integrate a 1-form over an Edge ---
    std::cout << "\n--- 2a. Integrate 1-Form over Edge e_01 ---" << std::endl;
    // Construct the Whitney 1-form for edge e_01: φ_01 = λ₀*dλ₁ - λ₁*dλ₀
    // This is a P_1Λ¹(T) form, but part of the P_1⁻Λ¹(T) basis
    DifferentialForm phi_01 = lambda(0) * d_lambda(1) + lambda(1) * d_lambda(0) * Polynomial(-1.0);
    
    std::cout << "Whitney 1-form φ_01 = " << phi_01.toString() << std::endl;
    // Note: The toString() shows it in the independent basis:
    // (λ0) * [dλ1] + (-1*λ1) * [dλ1] + (-1*λ1) * [dλ2] + (-1*λ1) * [dλ3]
    // = (λ0 - λ1) * [dλ1] - λ1 * [dλ2] - λ1 * [dλ3]

    // Integrate φ_01 over its own edge e_01
    // This integral should be 1/k! = 1/1! = 1 (or -1 depending on orientation)
    // "it follows from (4.2) that ∫_f φσ = ±1/k!"
    // Let's see what our code gets.
    double int_phi_01_on_e_01 = phi_01.integrate(e_01);
    std::cout << "  Integrate ∫_e_01 φ_01: " << int_phi_01_on_e_01 << std::endl;

    // --- 2b. Integrate a 2-form over a Face ---
    std::cout << "\n--- 2b. Integrate 2-Form over Face f_012 ---" << std::endl;
    // Construct the Whitney 2-form for face f_012:
    // φ_012 = λ₀(dλ₁∧dλ₂) - λ₁(dλ₀∧dλ₂) + λ₂(dλ₀∧dλ₁)
    DifferentialForm t0 = lambda(0) * wedge(d_lambda(1), d_lambda(2));
    DifferentialForm t1 = lambda(1) * wedge(d_lambda(0), d_lambda(2)) * Polynomial(-1.0);
    DifferentialForm t2 = lambda(2) * wedge(d_lambda(0), d_lambda(1));
    DifferentialForm phi_012 = t0 + t1 + t2;

    std::cout << "Whitney 2-form φ_012 = " << phi_012.toString() << std::endl;
    
    // Integrate φ_012 over its own face f_012
    // This integral should be ±1/k! = ±1/2! = ±0.5
    double int_phi_012_on_f_012 = phi_012.integrate(f_012);
    std::cout << "  Integrate ∫_f_012 φ_012: " << int_phi_012_on_f_012 << std::endl;

    return 0;
}