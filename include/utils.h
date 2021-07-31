#ifndef UTILS_H
#define UTILS_H
// C++ Standard Libraries Include
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

// Boost Libraries Include
#include <boost/math/special_functions/bessel.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>

// QT Libraries Include
#include <QtCore/QVector>

// Defining namespaces and typedefs
namespace bg = boost::geometry;
namespace bm = boost::math;

typedef bg::model::point<float, 2, bg::cs::cartesian> point_t;
typedef bg::model::linestring<point_t> linestring_t;
typedef double ld_t;

namespace utils
{
    // Conversion between QVector and vector<point_t>
    std::vector<point_t> QVec_to_Vec(QVector<double> X, QVector<double> Y);
    void Vec_to_QVec(std::vector<point_t> vec, QVector<double> *X, QVector<double> *Y);
    // Calculates Intersection between two linestrings
    std::vector<point_t> Intersection_point(std::vector<point_t> a, std::vector<point_t> b);

    // Bessel functions
    ld_t J0(ld_t in);
    ld_t J1(ld_t in);
    ld_t J2(ld_t in);
    ld_t K0(ld_t in);
    ld_t K1(ld_t in);
    ld_t K2(ld_t in);
    ld_t N0(ld_t in);
    ld_t N1(ld_t in);
    ld_t N2(ld_t in);

    // Helper Functions
    int lesser_index(double value, QVector<double> vec);
}

#endif // UTILS_H
