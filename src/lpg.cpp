#include "lpg.h"

// Constructors
LPG::LPG()
{
}

LPG::LPG(ld_t Ncore, ld_t Nclad, ld_t Next, ld_t Rcore, ld_t Rclad, ld_t GratingPeriod, ld_t AverageIndexChange, ld_t Length)
{
    n_core = Ncore;
    n_clad = Nclad;
    n_ext = Next;
    r_core = Rcore;
    r_clad = Rclad;
    grating_period = GratingPeriod;
    average_index_change = AverageIndexChange;
    refractive_index_contrast = (r_core * r_core - r_clad * r_clad) / (2 * r_core * r_core);
    length = Length;
}

// Setters
void LPG::set_n_core(ld_t NCore)
{
    n_core = NCore;
}

void LPG::set_n_clad(ld_t NClad)
{
    n_clad = NClad;
}

void LPG::set_n_ext(ld_t NExt)
{
    n_ext = NExt;
}

void LPG::set_r_core(ld_t RCore)
{
    r_core = RCore;
}

void LPG::set_r_clad(ld_t RClad)
{
    r_clad = RClad;
}

void LPG::set_grating_period(ld_t GP)
{
    grating_period = GP;
}

void LPG::set_average_index_change(ld_t AIC)
{
    average_index_change = AIC;
}

void LPG::set_refractive_index_contrast()
{
    try {
        refractive_index_contrast = (r_core * r_core - r_clad * r_clad) / (2 * r_core * r_core);
    }  catch (...) {
        std::cout << "ERROR: Could not set value for 'refractive_index_contrast' \n Try calling this function after setting the values of all other variables in LPG class.\n Terminating Application...\n";
    }
}

// Helper Methods


ld_t LPG::normalised_frequency(ld_t Wavelength, ld_t r_int, ld_t n_int, ld_t n_out)
{
    /*
     * Returns normalised frequency for a given wavelength
     * Out:
     *      V = Normalised Frequency
     * In:
     *      Wavelength = Input Wavelength
     *      r_int = Radius of inner fibre
     *      n_int = refractive index of fibre
     *      n_out = refractive index of outer suface
     */
    ld_t V = (2 * M_PI * r_int / Wavelength) * sqrt(n_int * n_int - n_out * n_out);
    return V;
}

ld_t LPG::dispersion_lhs(ld_t U)
{
    try
    {
    return U * (utils::J1(U) / utils::J0(U));
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::dispersion_rhs(ld_t W)
{
    try
    {
        return W * (utils::K1(W) / utils::K0(W));
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::beta(ld_t Wavelength, ld_t N, ld_t R, ld_t U)
{
    ld_t t1, t2;
    t1 = (2 * M_PI * N / Wavelength);
    t2 = U / R;
    return sqrt(t1 * t1 - t2 * t2);
}

// Effective Refractive index
ld_t LPG::n_core_eff(ld_t Wavelength)
{
    ld_t V = LPG::normalised_frequency(Wavelength, r_core, n_core, n_clad);
    QVector<double> Uco, Wco, D_lhs, D_rhs, X_IP, Y_IP;
    double uco, wco, dlhs, drhs;
    for (int i = 1; i <= RESOLUTION; i++)
    {
        uco = i * V / RESOLUTION;
        wco = sqrt(V * V - uco * uco);
        dlhs = dispersion_lhs(uco);
        drhs = dispersion_rhs(wco);
        if (!isnan(dlhs) & !isnan(drhs))
        {
            Uco.append(uco);
            Wco.append(wco);
            D_lhs.append(dlhs);
            D_rhs.append(drhs);
        }
    }
    // test code for plotting
    test_U = Uco;
    test_W = Wco;
    test_D_lhs = D_lhs;
    test_D_rhs = D_rhs;

    // Finding intersection points
    std::vector<point_t> D_lhs_line = utils::QVec_to_Vec(Uco, D_lhs);
    std::vector<point_t> D_rhs_line = utils::QVec_to_Vec(Uco, D_rhs);

    std::vector<point_t> IP = utils::Intersection_point(D_lhs_line, D_rhs_line);
    utils::Vec_to_QVec(IP, &X_IP, &Y_IP);

    // test code for plotting
    test_X_IP = X_IP;
    test_Y_IP = Y_IP;



    // Calculating Effective Refractive index of the core
    ld_t UCO_min = X_IP.value(0);
    for (auto item: X_IP)
    {
        if (item < UCO_min)
        {
            UCO_min = item;
        }
    }
    ld_t B_co = beta(Wavelength, n_core, r_core, UCO_min);
    ld_t N_co_eff = (B_co * Wavelength) / (2 * M_PI);
    return N_co_eff;
}

QVector<double> LPG::n_clad_eff(ld_t Wavelength)
{
    QVector<double> Ucl, Wcl, D_lhs, D_rhs, X_IP, Y_IP, U_IP, D_lhs_test, N_clad_eff;
    ld_t V = LPG::normalised_frequency(Wavelength, r_clad, n_clad, n_ext);
    double ucl, wcl, dlhs, drhs;
    for (int i = 1; i <= RESOLUTION; i++)
    {
        ucl = i * V / RESOLUTION;
        wcl = sqrt(V * V - ucl * ucl);
        dlhs = LPG::dispersion_lhs(ucl);
        drhs = LPG::dispersion_rhs(wcl);
        if (!isnan(dlhs) & !isnan(drhs))
        {
            Ucl.append(ucl);
            Wcl.append(wcl);
            D_lhs.append(dlhs);
            D_rhs.append(drhs);
        }
    }
    // test code for plotting
    test_U = Ucl;
    test_W = Wcl;
    test_D_lhs = D_lhs;
    test_D_rhs = D_rhs;

    //Finding intersection points
    std::vector<point_t> D_lhs_line = utils::QVec_to_Vec(Ucl, D_lhs);
    std::vector<point_t> D_rhs_line = utils::QVec_to_Vec(Ucl, D_rhs);

    std::vector<point_t> IP = utils::Intersection_point(D_lhs_line, D_rhs_line);
    utils::Vec_to_QVec(IP, &X_IP, &Y_IP);

    // test code for plotting
    test_X_IP = X_IP;
    test_Y_IP = Y_IP;

    // Sorting IP and finding N clad eff
    std::sort(X_IP.begin(),X_IP.end());
    for (int i = 0; i < 50; i += 2)
    {
        U_IP.append(X_IP.value(i));
        for (auto item: IP)
        {
            if (item.get<0>() == X_IP.value(i))
            {
                D_lhs_test.append(item.get<1>());
            }
        }
    }
    test_Ucl = U_IP;
    test_Ucl_D = D_lhs_test;

    double n,b;
    for (auto u: U_IP)
    {
        b = sqrt((2 * M_PI * n_clad / Wavelength) * (2 * M_PI * n_clad / Wavelength) - (u / r_clad) * (u / r_clad));
        n = b * Wavelength / (2 * M_PI);
        N_clad_eff.append(n);
    }
    return N_clad_eff;
}

ld_t LPG::sigma_1(ld_t n_eff)
{
    return n_eff / Z0;
}

ld_t LPG::sigma_2(ld_t n_eff)
{
    return n_eff * Z0;
}

ld_t LPG::U1(ld_t Wavelength, ld_t n_eff)
{
    return (2 * M_PI / Wavelength) * sqrt(n_core * n_core - n_eff * n_eff);
}

ld_t LPG::U2(ld_t Wavelength, ld_t n_eff)
{
    return (2 * M_PI / Wavelength) * sqrt(n_clad * n_clad - n_eff * n_eff);
}

ld_t LPG::W3(ld_t Wavelength, ld_t n_eff)
{
    return (2 * M_PI / Wavelength) * sqrt(n_eff * n_eff - n_ext * n_ext);
}

ld_t LPG::U21(ld_t Wavelength, ld_t n_eff)
{
    ld_t u1 = LPG::U1(Wavelength, n_eff);
    ld_t u2 = LPG::U2(Wavelength, n_eff);
    return (1 / (u2 * u2)) - (1 / (u1 * u1));
}

ld_t LPG::U32(ld_t Wavelength, ld_t n_eff)
{
    ld_t u2 = LPG::U2(Wavelength, n_eff);
    ld_t w3 = LPG::W3(Wavelength, n_eff);
    return (1 / (w3 * w3)) - (1 / (u2 * u2));
}

ld_t LPG::J(ld_t Wavelength, ld_t n_eff)
{
    try {
        ld_t u1 = LPG::U1(Wavelength, n_eff);
        ld_t j = 0.5 * (utils::J0(u1 * r_core) - utils::J2(u1 * r_core)) / (u1 * utils::J1(u1 * r_core));
        return j;
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::K(ld_t Wavelength, ld_t n_eff)
{
    try {
        ld_t w3 = LPG::W3(Wavelength, n_eff);
        ld_t k = -0.5 * (utils::K0(w3 * r_clad) + utils::K2(w3 * r_clad)) / (w3 * utils::K1(w3 * r_clad));
        return k;
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::P(ld_t Wavelength, ld_t n_eff, ld_t r)
{
    try {
        ld_t u2 = LPG::U2(Wavelength, n_eff);
        ld_t p = utils::J1(u2 * r) * utils::N1(u2 * r_core) - utils::N1(u2 * r) * utils::J1(u2 * r_core);
        return p;
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::Q(ld_t Wavelength, ld_t n_eff, ld_t r)
{
    try {
        ld_t u2 = LPG::U2(Wavelength, n_eff);
        ld_t q = 0.5 * ((utils::J1(u2 * r) * (utils::N0(u2 * r_core) - utils::N2(u2 * r_core))) - (utils::N1(u2 * r) * (utils::J0(u2 * r_core) - utils::J2(u2 * r_core))));
        return q;
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::R(ld_t Wavelength, ld_t n_eff, ld_t r)
{
    try {
        ld_t u2 = LPG::U2(Wavelength, n_eff);
        ld_t r1 = utils::N1(u2 * r_core) * (utils::J0(u2 * r) - utils::J2(u2 * r));
        ld_t r2 = utils::J1(u2 * r_core) * (utils::N0(u2 * r) - utils::N2(u2 * r));
        return 0.5 * (r1 - r2);
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::S(ld_t Wavelength, ld_t n_eff, ld_t r)
{
    try {
        ld_t u2 = LPG::U2(Wavelength, n_eff);
        ld_t s1 = (utils::J0(u2 * r) - utils::J2(u2 * r)) * (utils::N0(u2 * r_core) - utils::N2(u2 * r_core));
        ld_t s2 = (utils::N0(u2 * r) - utils::N2(u2 * r)) * (utils::J0(u2 * r_core) - utils::J2(u2 * r_core));
        return 0.25 * (s1 - s2);
    }
    catch (...)
    {
        return std::nan("");
    }
}

ld_t LPG::epsilon_o(ld_t Wavelength, ld_t n_eff)
{
    // Variable initialisation
    ld_t sigma1 = LPG::sigma_1(n_eff);
    ld_t sigma2 = LPG::sigma_2(n_eff);
    // ld_t u1 = LPG::U1(Wavelength, n_eff);
    ld_t u2 = LPG::U2(Wavelength, n_eff);
    // ld_t w3 = LPG::W3(Wavelength, n_eff);
    ld_t u21 = LPG::U21(Wavelength, n_eff);
    ld_t u32 = LPG::U32(Wavelength, n_eff);
    ld_t j = LPG::J(Wavelength, n_eff);
    ld_t k = LPG::K(Wavelength, n_eff);
    ld_t p = LPG::P(Wavelength, n_eff, r_clad);
    ld_t q = LPG::Q(Wavelength, n_eff, r_clad);
    ld_t r = LPG::R(Wavelength, n_eff, r_clad);
    ld_t s = LPG::S(Wavelength, n_eff, r_clad);

    // Numerator terms
    ld_t n1 = u2 * (j * k - ((sigma1 * sigma2 * u21 * u32) / (n_clad * n_clad * r_core * r_clad))) * p; // negative sign due to sigma 1 and sigma 2 being complex numbers
    ld_t n2 = -(k * q) + (j * r) - (s / u2);
    // Denominator terms
    ld_t d1 = -u2 * ((u32 * j / (n_clad * n_clad * r_clad)) - (u21 * k / (n_core * n_core * r_core))) * p;
    ld_t d2 = u32 * q / (n_core * n_core * r_clad);
    ld_t d3 = u21 * r / (n_core * n_core * r_core);

    ld_t out = (-1 / sigma2) * (n1 + n2) / (d1 + d2 + d3);// -1/sigma2 due to sigma being complex
    return out;
}

ld_t LPG::epsilon_o1(ld_t Wavelength, ld_t n_eff)
{
    // Variable initialisation
    ld_t sigma1 = LPG::sigma_1(n_eff);
    ld_t sigma2 = LPG::sigma_2(n_eff);
    // ld_t u1 = LPG::U1(Wavelength, n_eff);
    ld_t u2 = LPG::U2(Wavelength, n_eff);
    // ld_t w3 = LPG::W3(Wavelength, n_eff);
    ld_t u21 = LPG::U21(Wavelength, n_eff);
    ld_t u32 = LPG::U32(Wavelength, n_eff);
    ld_t j = LPG::J(Wavelength, n_eff);
    ld_t k = LPG::K(Wavelength, n_eff);
    ld_t p = LPG::P(Wavelength, n_eff, r_clad);
    ld_t q = LPG::Q(Wavelength, n_eff, r_clad);
    ld_t r = LPG::R(Wavelength, n_eff, r_clad);
    ld_t s = LPG::S(Wavelength, n_eff, r_clad);

    // Numerator terms
    ld_t n1 = u2 * p * ((u32 * j / r_clad) - ((n_ext * n_ext * u21 * k) / (n_clad * n_clad * r_core)));
    ld_t n2 = - u32 * q / r_clad;
    ld_t n3 = - u21 * r / r_core;

    // Denominator terms
    ld_t d1 = u2 * p * ((n_ext * n_ext * j * k / (n_clad * n_clad)) - ((sigma1 * sigma2 * u21 * u32) / (n_core * n_core * r_core * r_clad)));
    ld_t d2 = - n_ext * n_ext * k * q / (n_core * n_core);
    ld_t d3 = j * r;
    ld_t d4 = - n_clad * n_clad * s / (n_core * n_core * u2);

    ld_t out = sigma1 * (n1 + n2 + n3) / (d1 + d2 + d3 + d4);
    return out;
}

ld_t LPG::delta_cl_co(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff)
{
    ld_t bco = n_core_eff * 2 * M_PI / Wavelength;
    ld_t bcl = n_clad_eff * 2 * M_PI / Wavelength;
    ld_t delta = 0.5 * (bco - bcl - (2 * M_PI / grating_period));
    return delta;
}


ld_t LPG::K_co_co(ld_t Wavelength, ld_t n_core_eff)
{
    ld_t b = (n_core_eff * n_core_eff - n_core * n_core) / (n_clad * n_clad - n_core * n_core);
    ld_t V = LPG::normalised_frequency(Wavelength, r_core, n_core, n_clad);
    ld_t d = (n_core - n_clad) / n_core;
    ld_t k = ((2 * M_PI * n_core * n_core * b) / (Wavelength * n_clad * sqrt(1 + 2 * b * d))) * (1 + (utils::J0(V * sqrt(1 - b)) * utils::J0(V * sqrt(1 - b))) / (utils::J1(V * sqrt(1 - b)) * utils::J1(V * sqrt(1 - b))));
    return k;
}

ld_t LPG::coupling_coefficient(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff)
{
    ld_t cc = LPG::delta_cl_co(Wavelength, n_core_eff, n_clad_eff) / (1  - (LPG::K_co_co(Wavelength, n_core_eff) / 2));
    return cc;
}

ld_t LPG::K_co_cl(ld_t Wavelength, ld_t n_core_eff, ld_t n_clad_eff)
{
    ld_t b = (n_core_eff * n_core_eff - n_core * n_core) / (n_clad * n_clad - n_core * n_core);
    ld_t V = LPG::normalised_frequency(Wavelength, r_core, n_core, n_clad);
    ld_t d = (n_core - n_clad) / n_core;
    ld_t u1 = r_core * sqrt((2 * M_PI * n_core / Wavelength) * (2 * M_PI * n_core / Wavelength) - (n_core_eff * 2 * M_PI / Wavelength) * (n_core_eff * 2 * M_PI / Wavelength));

    ld_t t1 = LPG::coupling_coefficient(Wavelength, n_core_eff, n_clad_eff) * 2 * M_PI / Wavelength;
    ld_t t2 = sqrt(M_PI * b / (Z0 * n_clad * sqrt(1 + 2 * b * d)));
    ld_t t3 = (n_core * n_core * u1) / (u1 * u1 - ((V * V * (1 - b)) / r_core * r_core));
    ld_t t4 = 1 + (LPG::sigma_2(n_clad_eff) * LPG::epsilon_o(Wavelength, n_clad_eff) / (n_core * n_core));
    ld_t t5 = (u1 * utils::J1(u1 * r_core) * utils::J0(V * sqrt(1 - b))) / (utils::J1(V * sqrt(1 - b)));
    ld_t t6 = (V * sqrt(1 - b) * utils::J0(u1 * r_core)) / r_core;
    ld_t t7 = sqrt(Z0 * b / (M_PI * n_clad_eff * sqrt(1 + 2 * b * d))) * (1 / (r_clad * utils::J1(V * sqrt(1 - b))));
    ld_t out = t1 * t2 * t3 * t4 * (t5 - t6) * t7;
    return out;
}







