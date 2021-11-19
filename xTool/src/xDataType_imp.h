#ifndef XDATATYPE__H
#error Do NOT include xDataType_imp.h alone
#endif

#include <complex>
#include <limits>

namespace xtool {
/// Specialization of xDatatype for float
template < >
class xDataType < float >
{
    public:
        static float zero() { return 0.f; }
        static float one() { return 1.f; }
        static float epsilonNumeric() { return std::numeric_limits < float >::epsilon(); }
        static std::string  stype() { return "real_single"; }
};

/// Specialization of xDatatype for double
template < >
class xDataType < double >
{
    public:
        static double zero() { return 0.; }
        static double one() { return 1.; }
        static double epsilonNumeric() { return std::numeric_limits < double >::epsilon(); }
        static std::string  stype() { return "real_double"; }

};

/// Specialization of xDatatype for double-complex
template < >
class xDataType < std::complex<double> >
{
    public:
    using cplx_t = std::complex<double>;
    static cplx_t zero() { return cplx_t(0.); }
    static cplx_t one() { return std::complex<double>(1., 0.); }
    static cplx_t epsilonNumeric() { return std::numeric_limits < cplx_t >::epsilon(); }
    static std::string  stype() { return "complex_double"; }
};
}
