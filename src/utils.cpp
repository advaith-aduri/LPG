#include "utils.h"


namespace utils
{
    // Conversion between QVector and vector<point_t>
    std::vector<point_t> QVec_to_Vec(QVector<double> X, QVector<double> Y)
    {
        std::vector<point_t> out;
        if (X.size() != Y.size())
        {
            std::cout << " X and Y are not of same size.\n Returning a empty vector.\n";
            return out;
        }
        for (int i = 0; i < X.size(); i++)
        {
            out.push_back(point_t(X.value(i), Y.value(i)));
        }
        return out;
    }

    void Vec_to_QVec(std::vector<point_t> vec, QVector<double> *X, QVector<double> *Y)
    {
        for (unsigned int i = 0; i < vec.size(); i++)
        {
            X->append(vec[i].get<0>());
            Y->append(vec[i].get<1>());
        }
    }

    // Function for finding points of intersection between 2 linestrings
    std::vector<point_t> Intersection_point(std::vector<point_t> a, std::vector<point_t> b)
    {
        linestring_t l1, l2;
        std::vector<point_t> intersection_points;

        for (auto item:a)
        {
            bg::append(l1, item);
        }

        for (auto item:b)
        {
            bg::append(l2, item);
        }

        bg::intersection(l1, l2, intersection_points);
        return intersection_points;
    }

    // Bessel functions of different kinds
    ld_t J0(ld_t in)
    {
        return bm::cyl_bessel_j(0, in);
    }

    ld_t J1(ld_t in)
    {
        return bm::cyl_bessel_j(1, in);
    }

    ld_t J2(ld_t in)
    {
        return bm::cyl_bessel_j(2, in);
    }

    ld_t K0(ld_t in)
    {
        return bm::cyl_bessel_k(0, in);
    }

    ld_t K1(ld_t in)
    {
        return bm::cyl_bessel_k(1, in);
    }

    ld_t K2(ld_t in)
    {
        return bm::cyl_bessel_k(2, in);
    }

    ld_t N0(ld_t in)
    {
        return bm::cyl_neumann(0, in);
    }

    ld_t N1(ld_t in)
    {
        return bm::cyl_neumann(1, in);
    }

    ld_t N2(ld_t in)
    {
        return bm::cyl_neumann(2, in);
    }

    // Helper Functions
    int lesser_index(double value, QVector<double> vec)
    {
        int i = 0;
        while (value > vec.value(i))
        {
            i++;
        }
        return i - 1;
    }

}
