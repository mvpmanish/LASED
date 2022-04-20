#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


//----------------------TYPE AND CONSTANT DEFINITIONS--------------------------
// Use Eigen boost library
// Dynamic matrix of complex doubles type
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;
// Dynamic matrix of 8 columns to store States with (label, w, m, L, S, J = None, I = None, F = None)
typedef Eigen::Matrix<double, Eigen::Dynamic, 8> States;

using namespace std::complex_literals;

// Define constants
std::vector<double> q_decay = { -1, 0, 1 }; // Selection rules
constexpr double PI = 3.14159265358979323846264338327950288;
constexpr double C = 299792458;  // Speed of light in m/s
constexpr double H = 6.62607015e-34;  // Planck's constant in Js

//---------------------------HELPER FUNCTIONS----------------------------------
/**
 * Calculates the half-Rabi frequency. This is a rewritten version of halfRabinFreq in half_rabi_freq.py.
 * @brief Calculates the half - Rabi frequency in Grad / s.
 * @param intensity The intensity of the laser in mW / mm ^ 2.
 * @param lifetime The lifetime of the excited state to the ground state transition in nanoseconds.
 * @param wavelength The wavelength(in metres) corresponding to the resonant transition from the ground to excited state.
 * @return The half - Rabi frequency in Grad / s.
 */
double halfRabiFreq(double intensity, double lifetime, double wavelength, double rabi_scaling)
{
    double I = intensity * 1000;  // Convert mW / mm ^ 2 to W / m ^ 2
    double h = 6.6260715e-34;  // Plank's constant
    double c = 299792458;  // Speed of light
    double tau = lifetime * 1e-9;  // Convert ns lfetime to s
    if (rabi_scaling == 0){rabi_scaling = 1}
    return std::sqrt((3 * I * std::pow(wavelength, 3)) / (8 * PI * c * h * tau)) * 1e-9 * rabi_scaling;
}

//---------------------------CLASS DEFINITIONS-------------------------
/**
 * Class containing the solution of the equations of motion for a laser-atom system.
 */
class Solution
{
public:
    MatrixXc V; /**< The matrix of eigenvectors of the laser-atom system's matrix A.*/
    MatrixXc inv_V;  /**< The inverse of the matrix of eigenvectors.*/
    MatrixXc D;  /**< The diagonalised form of the laser-atom system's matrix A.*/
    /**
     * A parameterised constructor.
     */
    Solution(MatrixXc V_, MatrixXc inv_V_, MatrixXc D_)
    { // Constructor with parameters
        V = V_;
        inv_V = inv_V_;
        D = D_;
    }

    static Solution create(MatrixXc x, MatrixXc y, MatrixXc z)
    {
        return Solution(x, y, z);
    }

    MatrixXc timeEvolution(MatrixXc rho0, double time)
    {
        if (static_cast<size_t>(rho0.cols() > 1))
        {
            rho0 = rho0.transpose();
        }
        //MatrixXc rhot = V * ((time * D).array().exp().matrix().asDiagonal()) * inv_V * rho0;
        MatrixXc rhot = ((V * ((time * D).array().exp().matrix().asDiagonal())) * inv_V) * rho0;
        // D is a vector in matrix form, transform it in to array form to perform
        // componentwise operation, then transform it back to matrix form and diagonal form.
        return  rhot;
    }
    MatrixXc get_D()
    {
        return D;
    }
    MatrixXc get_V()
    {
        return V;
    }
    MatrixXc get_inv_V()
    {
        return inv_V;
    }
};


class SolveLaserAtomSystem
{
public:
    States E0;
    States G;
    double tau;
    std::vector<int> Q;
    double laser_wavelength;
    Eigen::MatrixXd couplingtable_q1;
    Eigen::MatrixXd couplingtable_q0;
    Eigen::MatrixXd couplingtable_qn1;
    double detuning;
    double laser_intensity = INFINITY;// change this later
    double laser_power = INFINITY;
    double tau_f = INFINITY;
    double tau_b = INFINITY;
    std::vector<double> rabi_scaling = { 0 };
    std::vector<double> rabi_factors = { 0 };
    std::vector<Eigen::MatrixXd> couplingtable;
    size_t boundary;
    double rabi;
    std::vector<int> Q_decay = { -1, 0, 1 };
    double atomic_velocity;
    MatrixXc A;

    SolveLaserAtomSystem
    (States E_,
        States G_,
        double tau_,
        std::vector<int> Q_,
        double laser_wavelength_,
        Eigen::MatrixXd couplingtable_q1_,
        Eigen::MatrixXd couplingtable_q0_,
        Eigen::MatrixXd couplingtable_qn1_,
        double detuning_,
        double laser_intensity_ = 0.0,
        double laser_power_ = 0.0,
        double tau_f_ = 0.0,
        double tau_b_ = 0.0,
        std::vector<double> rabi_scaling_ = { 0 },
        std::vector<double> rabi_factors_ = { 0 },
        double atomic_velocity = 0
    )
    {
        E0 = E_;
        G = G_;
        tau = tau_;
        Q = Q_;
        laser_wavelength = laser_wavelength_;
        couplingtable_q1 = couplingtable_q1_;
        couplingtable_q0 = couplingtable_q0_;
        couplingtable_qn1 = couplingtable_qn1_;
        detuning = detuning_;
        laser_intensity = laser_intensity_;
        laser_power = laser_power_;
        tau_f = tau_f_;
        tau_b = tau_b_;
        rabi_scaling = rabi_scaling_;
        rabi_factors = rabi_factors_;
        // make one coupling table from the three tables
        //couplingtable = { couplingtable_qn1, couplingtable_q0, couplingtable_q1 };
        boundary = static_cast<size_t>(G.rows());
        rabi = halfRabiFreq(laser_intensity, tau, laser_wavelength, rabi_scaling[0]);
    }
    //factory function
    static SolveLaserAtomSystem create(States E_,
        States G_,
        double tau_,
        std::vector<int> Q_,
        double laser_wavelength_,
        Eigen::MatrixXd couplingtable_q1_,
        Eigen::MatrixXd couplingtable_q0_,
        Eigen::MatrixXd couplingtable_qn1_,
        double detuning_,
        double laser_intensity_ = 0.0,
        double laser_power_ = 0.0,
        double tau_f_ = 0.0,
        double tau_b_ = 0.0,
        std::vector<double> rabi_scaling_ = { 0 },
        std::vector<double> rabi_factors_ = { 0 },
        double atomic_velocity = 0
    )
    {
        return SolveLaserAtomSystem(E_,
            G_,
            tau_,
            Q_,
            laser_wavelength_,
            couplingtable_q1_,
            couplingtable_q0_,
            couplingtable_qn1_,
            detuning_,
            laser_intensity_,
            laser_power_ ,
            tau_f_,
            tau_b_ ,
            rabi_scaling_ ,
            rabi_factors_ ,
            atomic_velocity
        );
    }
    // debug

    double get_tau()
    {
        return tau;
    }
    std::vector<int> get_Q()
    {
        return Q;
    }


    // helper functions

    double angularFreq(double wavelength)
    {
        return 2 * PI * C / wavelength * 1e-9;
    }

    double coupling(size_t e, size_t g, int q)
    {
        if (q == -1)
        {
            return couplingtable_qn1(e, g);
        }
        else if (q == 0)
        {
            return couplingtable_q0(e, g);
        }
        else if (q == 1)
        {
            return couplingtable_q1(e, g);
        }
        else
        {
            return 0;
        }
    }
    double E(size_t e, size_t g)
    {
        //change index
        e -= boundary;
        return E0(e, g);
    }
    double generalisedDecayConstant(size_t ep, size_t epp, size_t g)
        /*Calculates the branching ratio for the generalised decay constant.

        This ratio must be multiplied by 1 / tau to get the generalised decay constant in Grad / s.This is then used for evaluating vertical coherences.

        Parameters:
        ep(State) : Excited State object
        epp(State) : Excited State object
        g(State) : Ground State object
        G(list) : : List of all ground State objects
        tau(float) : Lifetime of state in ns
        Q_decay(list) : List of decay channel polarisations allowed by transition rules, usually[-1, 0, 1]

        Returns :
            Value of the branching ratio for the  generalised decay constant gamma_{ ep, epp, g }.
        */
    {
        double sum_decay_channels_epg = 0;
        double sum_decay_channels_eppg = 0;
        for (size_t gp = 0; gp < static_cast<size_t>(G.rows()); ++gp)
        {
            for (size_t q = 0; q < Q_decay.size(); ++q)
            {
                sum_decay_channels_epg += abs(coupling(ep, gp, Q_decay[q]) * coupling(ep, gp, Q_decay[q]));
                sum_decay_channels_eppg += abs(coupling(epp, gp, Q_decay[q]) * coupling(epp, gp, Q_decay[q]));
            }
        }

        // Calculate the decay constants separately
        double gamma_epg = 0;
        double gamma_eppg = 0;
        for (size_t q = 0; q < Q_decay.size(); ++q)
        {
            gamma_epg += abs(coupling(ep, g, Q_decay[q]) * coupling(ep, g, Q_decay[q]));
            gamma_eppg += abs(coupling(epp, g, Q_decay[q]) * coupling(epp, g, Q_decay[q]));
        }
        gamma_epg = gamma_epg / (sum_decay_channels_epg);
        gamma_eppg = gamma_eppg / (sum_decay_channels_eppg);

        // Calculate the generalised decay constant
        double gamma_epeppg = std::pow(gamma_epg * gamma_eppg, 0.5);
        // Calculate the sign : only one polarisation will result in non - zero coupling so can sum over all
        for (size_t q = 0; q < Q_decay.size(); ++q)
        {
            if (coupling(ep, g, Q_decay[q]) * coupling(epp, g, Q_decay[q]) < 0)
            {
                gamma_epeppg = -1 * gamma_epeppg;

            }
        }
        return gamma_epeppg;
    }

    double dopplerDelta(size_t e, size_t g, double w_q, double lambda_q, double v_z)
    {
        double atomic_velocity_detuning = 2 * PI * v_z / (lambda_q * 1e9);
        return w_q - atomic_velocity_detuning - E(e, 1) + G(g, 1);
    }

    std::complex<double> ggpp(size_t i, size_t j, size_t k, size_t l)
    {
        // i, j are g, gpp ; k, l could be e or g.
        // ggpp term
        if (k == i && l == j)
        {
            if (tau_b == 0)
            {
                return -1.0i * (G(i, 1) - G(j, 1)); // delta term
            }
            else
            {

                return -1.0i * (G(i, 1) - G(j, 1)) - 1 / tau_b;

            }
        }
        // ge term k = g, l = e
        else if (k == i && l >= boundary)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(l, j, Q[q_]) * 1.0i * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        // eg term
        else if (k >= boundary && l == j)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(k, i, Q[q_]) * (-1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        // ee term
        else if (k >= boundary && l >= boundary)
        {
            // change index k = ep, l = epp, i = g j = gpp
            double sum_decay_channels = 0;
            for (size_t q_ = 0; q_ < Q_decay.size(); ++q_)
            {
                for (size_t g_ = 0; g_ < static_cast<size_t>(G.rows()); ++g_)
                {
                    sum_decay_channels += std::abs(coupling(l, g_, Q_decay[q_]) * coupling(k, g_, Q_decay[q_]));
                }
            }
            if (sum_decay_channels != 0)
            {
                if (k == l)
                {
                    double result = 0;
                    for (size_t q_ = 0; q_ < Q_decay.size(); ++q_)
                    {
                        result += std::abs(coupling(k, j, Q_decay[q_]) * coupling(l, i, Q_decay[q_]));
                    }
                    result /= sum_decay_channels * tau;
                    return  result;
                }
                else
                {
                    double decay_const = generalisedDecayConstant(k, l, j) / tau;
                    return decay_const;
                }
            }
            else
            {
                return 0.0 + 0.0i;
            }
        }
        else
        {
            return 0.0 + 0.0i;
        }
        //else is zero
    }

    std::complex<double> eepp(size_t i, size_t j, size_t k, size_t l)
    {
        // i, j are e, epp ; k, l could be e or g.
        // eepp term

        if (k == i && l == j)
        {
            if (tau_f == 0.0)
            {
                return -1.0i * (E(i, 1) - E(j, 1)) - 1.0 / tau; // delta term and decay term
            }
            else
            {
                return -1.0i * (E(i, 1) - E(j, 1)) - 1.0 / tau - 1.0 / tau_f;
            }
        }
        // eg term
        else if (k == i && l < boundary)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(j, l, Q[q_]) * 1.0i * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        // gepp term
        else if (k < boundary && l == j)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(i, k, Q[q_]) * (-1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        //else is zero
        else
        {
            return 0;
        }
    }

    std::complex<double> ge(size_t i, size_t j, size_t k, size_t l)
    {
        // i, j are g, e ; k, l could be e or g.
        // ge term.
        if (k == i && l == j)
        {
            // this term also contains detuning
            std::complex<double> result = -1.0i * dopplerDelta(j, i, angularFreq(laser_wavelength),
                laser_wavelength, atomic_velocity) - 1 / (2 * tau);
            if (detuning != 0.0)
            {
                result += -1.0i * detuning; // detuning term
            }
            if (tau_f != 0.0)
            {
                result += -1.0 / (2.0* tau_f); // decay
            }
            if (tau_b != 0.0)
            {
                result += -1.0 / (2.0 * tau_b); // decay
            }
            return result;
        }
        // gg term
        else if (k == i && l < boundary)
        {

            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(j, l, Q[q_]) * (1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        // eep term
        else if (k >= boundary && l == j)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(k, i, Q[q_]) * (-1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        //else is zero
        else
        {
            return 0;
        }
    }
    std::complex<double> eg(size_t i, size_t j, size_t k, size_t l)
    {
        // i, j are e, g ; k, l could be e or g.
        // ge term.
        if (k == i && l == j)
        {
            // this term also contains detuning
            std::complex<double> result = 1.0i * dopplerDelta(i, j, angularFreq(laser_wavelength),
                laser_wavelength, atomic_velocity) - 1 / (2 * tau);
            if (detuning != 0.0)
            {
                result += 1.0i * detuning; // detuning term
            }
            if (tau_f != 0.0)
            {
                result += -1.0 / (2.0 * tau_f); // decay
            }
            if (tau_b != 0.0)
            {
                result += -1.0 / (2.0 * tau_b); // decay
            }
            return result;
        }
        // gg term
        else if (k < boundary && l == j)
        {

            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(i, k, Q[q_]) * (-1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        // eep term
        else if (k == i && l >= boundary)
        {
            std::complex<double> currSum = 0.0 + 0.0i;
            for (size_t q_ = 0; q_ < Q.size(); ++q_)
            {
                currSum += coupling(l, j, Q[q_]) * (1.0i) * rabi * rabi_factors[q_];
            }
            return currSum;
        }
        //else is zero
        else
        {
            return 0;
        }
    }


    // Generating A matrix

    MatrixXc AMatrix()
    {
        //size_t boundary = static_cast<size_t>(G.rows());
        size_t n = boundary + static_cast<size_t>(E0.rows());
        size_t N = n * n;
        MatrixXc A = MatrixXc::Zero(N, N);
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                for (size_t k = 0; k < n; ++k)
                {
                    for (size_t l = 0; l < n; ++l)
                    {
                        size_t I = i * n + j;
                        size_t J = k * n + l;
                        size_t I1 = j * n + i;


                        if (i < boundary && j < boundary)
                        {
                            A(I, J) = ggpp(i, j, k, l);
                        }
                        else if (i < boundary && j >= boundary)
                        {
                            A(I, J) = ge(i, j, k, l);
                            //A(I1, J) = std::conj(A(I, J));
                        }
                        else if (i >= boundary && j < boundary)
                        {
                            A(I, J) = eg(i, j, k, l);
                        }
                        else if (i >= boundary && j >= boundary)
                        {
                            A(I, J) = eepp(i, j, k, l);
                        }
                    }
                }
            }
        }
        return A;
    }

    // Generate solution
    Solution get_solution()
    {
        Eigen::ComplexEigenSolver<MatrixXc> ces;
        MatrixXc A_ = AMatrix();
        ces.compute(A_);
        auto V = ces.eigenvectors();
        MatrixXc inv_V = V.inverse();
        MatrixXc D = ces.eigenvalues(); //vector
        Solution solution(V, inv_V, D);
        return solution;
    }
    MatrixXc get_solution0()
    {
        Eigen::ComplexEigenSolver<MatrixXc> ces;
        MatrixXc A_ = AMatrix();
        ces.compute(A_);
        auto V = ces.eigenvectors();
        MatrixXc inv_V = V.inverse();
        MatrixXc D = ces.eigenvalues(); //vector
        Solution solution(V, inv_V, D);
        return A_;
    }

\

};

namespace py = pybind11;

PYBIND11_MODULE(CppLASED, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: python_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    //m.def("testmatrix", &nothing, R"pbdoc(
    //    test complex matrix to array conversion
    //)pbdoc");
    py::class_<Solution>(m, "Solution")
        //.def(py::init<const std::string&>())
        .def(py::init(&Solution::create))
        .def("get_D", &Solution::get_D)
        .def("get_V", &Solution::get_V)
        .def("get_inv_V", &Solution::get_inv_V)
        .def("timeEvolution", &Solution::timeEvolution);
        //.def("getName", &Pet::getName);

    py::class_<SolveLaserAtomSystem>(m, "SolveLaserAtomSystem")
        .def(py::init(&SolveLaserAtomSystem::create))
        .def("get_tau", &SolveLaserAtomSystem::get_tau)
        .def("get_Q", &SolveLaserAtomSystem::get_Q)
        //.def("ggpp", &SolveLaserAtomSystem::ggpp)
        .def("get_solution0", &SolveLaserAtomSystem::get_solution0)
        .def("get_solution", &SolveLaserAtomSystem::get_solution);

    //.def("getName", &Pet::getName);



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
