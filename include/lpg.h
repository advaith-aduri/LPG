#ifndef LPG_H
#define LPG_H

#include "LPG_global.h"
#include "utils.h"

#define RESOLUTION 1000000

class LPG_EXPORT LPG
{
    // Data
    ld_t n_core, n_clad, n_ext;
    ld_t r_core, r_clad;
    ld_t grating_period, average_index_change, refractive_index_contrast, length;
    ld_t Z0 = 377;
    // Methods
    // Helper Methods
    ld_t normalised_frequency(ld_t Wavelength, ld_t r_int, ld_t n_int, ld_t n_out);
    ld_t dispersion_lhs(ld_t U);
    ld_t dispersion_rhs(ld_t W);
    ld_t beta(ld_t Wavelength, ld_t N, ld_t R, ld_t U);

    // N cladding eff functions
    ld_t sigma_1(ld_t n_eff);
    ld_t sigma_2(ld_t n_eff);
    ld_t U1(ld_t Wavelength, ld_t n_eff);
    ld_t U2(ld_t Wavelength, ld_t n_eff);
    ld_t W3(ld_t Wavelength, ld_t n_eff);
    ld_t U21(ld_t Wavelength, ld_t n_eff);
    ld_t U32(ld_t Wavelength, ld_t n_eff);
    ld_t J(ld_t Wavelength, ld_t n_eff);
    ld_t K(ld_t Wavelength, ld_t n_eff);
    ld_t P(ld_t Wavelength, ld_t n_eff, ld_t r);
    ld_t Q(ld_t Wavelength, ld_t n_eff, ld_t r);
    ld_t R(ld_t Wavelength, ld_t n_eff, ld_t r);
    ld_t S(ld_t Wavelength, ld_t n_eff, ld_t r);
    ld_t epsilon_o(ld_t Wavelength, ld_t n_eff);
    ld_t epsilon_o1(ld_t Wavelength, ld_t n_eff);
    // Transmission Spectrum functions
    ld_t delta_cl_co(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff);
    ld_t K_co_co(ld_t Wavelength, ld_t n_core_eff);
    ld_t K_co_cl(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff);
    ld_t coupling_coefficient(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff);
public:
    // Methods
    // Constructor
    LPG();
    LPG(ld_t Ncore, ld_t Nclad, ld_t Next, ld_t Rcore, ld_t Rclad, ld_t GratingPeriod, ld_t AverageIndexChange, ld_t Length);

    // Setters
    void set_n_core(ld_t NCore);
    void set_n_clad(ld_t NClad);
    void set_n_ext(ld_t NExt);
    void set_r_core(ld_t RCore);
    void set_r_clad(ld_t RClad);
    void set_grating_period(ld_t GP);
    void set_average_index_change(ld_t AIC);
    void set_refractive_index_contrast(); // To be used only after all values are set

    // Effective Refractive Index : Need to move to private after testing
    ld_t n_core_eff(ld_t Wavelength);
    QVector<double> n_clad_eff(ld_t Wavelength);
    // Test Variables
    QVector<double> test_U, test_W, test_D_lhs, test_D_rhs, test_X_IP, test_Y_IP, test_Ucl, test_Ucl_D, test_T;
};

#endif // LPG_H
