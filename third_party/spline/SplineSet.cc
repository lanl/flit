/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2016                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Splines.hh"
#include <limits>
#include <cmath>

#ifdef __GCC__
#pragma GCC diagnostic ignored "-Wc++98-compat"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wimplicit-fallthrough"
#endif

/**
 *
 */

namespace Splines
{

using std::abs;
using std::sqrt;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//! spline constructor
SplineSet::SplineSet(string const & name)
    : _name(name)
    , baseValue(name + "_values")
    , basePointer(name + "_pointers")
    , _npts(0)
    , _nspl(0)
    , _X(nullptr)
    , _Y(nullptr)
    , _Yp(nullptr)
    , _Ypp(nullptr)
    , _Ymin(nullptr)
    , _Ymax(nullptr)
    , lastInterval(0)
{}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//! spline destructor
SplineSet::~SplineSet()
{
    baseValue.free();
    basePointer.free();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::info(ostream_type & s) const
{
    s << "SplineSet[" << name() << "] n.points = "
      << _npts << " n.splines = " << _nspl << '\n';
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        s << "\nSpline n." << i;
        switch(is_monotone[i])
        {
            case -2:
                s << " with NON monotone data\n";
                break;
            case -1:
                s << " is NOT monotone\n";
                break;
            case  0:
                s << " is monotone\n";
                break;
            case  1:
                s << " is strictly monotone\n";
                break;
            default:
                SPLINE_ASSERT(false, "SplineSet::info classification: " << is_monotone[i] << " not in range {-2,-1,0,1}");
        }
        splines[i]->info(s);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::dump_table(ostream_type & stream, integer num_points) const
{
    vector<real_type> vals;
    stream << 's';
    for(integer i = 0; i < numSplines(); ++i)
    {
        stream << '\t' << header(i);
    }
    stream << '\n';

    for(integer j = 0; j < num_points; ++j)
    {
        real_type s = xMin() + ((xMax() - xMin()) * j) / (num_points - 1);
        this->eval(s, vals);
        stream << s;
        for(integer i = 0; i < numSplines(); ++i)
        {
            stream << '\t' << vals[size_t(i)];
        }
        stream << '\n';
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

integer
SplineSet::getPosition(char const * hdr) const
{
    map<string, integer>::const_iterator it = header_to_position.find(hdr);
    SPLINE_ASSERT(it != header_to_position.end(), "SplineSet::getPosition(\"" << hdr << "\") not found!");
    return it->second;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::build(integer          nspl,
                 integer          npts,
                 char       const *headers[],
                 SplineType const stype[],
                 real_type  const X[],
                 real_type  const *Y[],
                 real_type  const *Yp[])
{
    SPLINE_ASSERT(nspl > 0,
                  "SplineSet::build expected positive nspl = " << nspl);
    SPLINE_ASSERT(npts > 1,
                  "SplineSet::build expected npts = " << npts <<
                  " greather than 1");
    _nspl = nspl;
    _npts = npts;
    // allocate memory
    splines.resize(size_t(_nspl));
    is_monotone.resize(size_t(_nspl));
    integer mem = npts;
    for(integer spl = 0; spl < nspl; ++spl)
    {
        switch(stype[size_t(spl)])
        {
            case QUINTIC_TYPE:
                mem += npts; // Y, Yp, Ypp
            case CUBIC_TYPE:
            case AKIMA_TYPE:
            case BESSEL_TYPE:
            case PCHIP_TYPE:
            case HERMITE_TYPE:
                mem += npts; // Y, Yp
            case CONSTANT_TYPE:
            case LINEAR_TYPE:
                mem += npts;
                break;
            case SPLINE_SET_TYPE:
            case SPLINE_VEC_TYPE:
                //default:
                SPLINE_ASSERT(false,
                              "SplineSet::build\nAt spline n. " << spl <<
                              " named " << headers[spl] <<
                              " cannot be done for type = " << stype[spl]);
        }
    }

    baseValue   . allocate(size_t(mem + 2 * nspl));
    basePointer . allocate(size_t(3 * nspl));

    _Y    = basePointer(size_t(_nspl));
    _Yp   = basePointer(size_t(_nspl));
    _Ypp  = basePointer(size_t(_nspl));
    _X    = baseValue(size_t(_npts));
    _Ymin = baseValue(size_t(_nspl));
    _Ymax = baseValue(size_t(_nspl));

    std::copy(X, X + npts, _X);
    for(size_t spl = 0; spl < size_t(nspl); ++spl)
    {
        real_type *& pY   = _Y[spl];
        real_type *& pYp  = _Yp[spl];
        real_type *& pYpp = _Ypp[spl];
        pY = baseValue(size_t(_npts));
        std::copy(Y[spl], Y[spl] + npts, pY);
        if(stype[spl] == CONSTANT_TYPE)
        {
            _Ymin[spl] = *std::min_element(pY, pY + npts - 1);
            _Ymax[spl] = *std::max_element(pY, pY + npts - 1);
        }
        else
        {
            _Ymin[spl] = *std::min_element(pY, pY + npts);
            _Ymax[spl] = *std::max_element(pY, pY + npts);
        }
        pYpp = pYp = nullptr;
        switch(stype[size_t(spl)])
        {
            case QUINTIC_TYPE:
                pYpp = baseValue(size_t(_npts));
            case CUBIC_TYPE:
            case AKIMA_TYPE:
            case BESSEL_TYPE:
            case PCHIP_TYPE:
            case HERMITE_TYPE:
                pYp = baseValue(size_t(_npts));
                if(stype[spl] == HERMITE_TYPE)
                {
                    SPLINE_ASSERT(Yp != nullptr && Yp[spl] != nullptr,
                                  "SplineSet::build\nAt spline n. " << spl <<
                                  " named " << headers[spl] <<
                                  "\nexpect to find derivative values");
                    std::copy(Yp[spl], Yp[spl] + npts, pYp);
                }
            case CONSTANT_TYPE:
            case LINEAR_TYPE:
            case SPLINE_SET_TYPE:
            case SPLINE_VEC_TYPE:
                //default:
                break;
        }
        string h = headers[spl];
        Spline * & s = splines[spl];

        is_monotone[spl] = -1;
        switch(stype[size_t(spl)])
        {
            case CONSTANT_TYPE:
                s = new ConstantSpline(h);
                static_cast<ConstantSpline*>(s)->reserve_external(_npts, _X, pY);
                static_cast<ConstantSpline*>(s)->npts = _npts;
                static_cast<ConstantSpline*>(s)->build();
                break;

            case LINEAR_TYPE:
                s = new LinearSpline(h);
                static_cast<LinearSpline*>(s)->reserve_external(_npts, _X, pY);
                static_cast<LinearSpline*>(s)->npts = _npts;
                static_cast<LinearSpline*>(s)->build();
                // check monotonicity of data
                {
                    integer flag = 1;
                    for(integer j = 1; j < _npts; ++j)
                    {
                        if(pY[j - 1] > pY[j])
                        {
                            flag = -1;    // non monotone data
                            break;
                        }
                        if(isZero(pY[j - 1] - pY[j]) && _X[j - 1] < _X[j])
                        {
                            flag = 0;    // non strict monotone
                        }
                    }
                    is_monotone[spl] = flag;
                }
                break;

            case CUBIC_TYPE:
                s = new CubicSpline(h);
                static_cast<CubicSpline*>(s)->reserve_external(_npts, _X, pY, pYp);
                static_cast<CubicSpline*>(s)->npts = _npts;
                static_cast<CubicSpline*>(s)->build();
                is_monotone[spl] = checkCubicSplineMonotonicity(_X, pY, pYp, _npts);
                break;

            case AKIMA_TYPE:
                s = new AkimaSpline(h);
                static_cast<AkimaSpline*>(s)->reserve_external(_npts, _X, pY, pYp);
                static_cast<AkimaSpline*>(s)->npts = _npts;
                static_cast<AkimaSpline*>(s)->build();
                is_monotone[spl] = checkCubicSplineMonotonicity(_X, pY, pYp, _npts);
                break;

            case BESSEL_TYPE:
                s = new BesselSpline(h);
                static_cast<BesselSpline*>(s)->reserve_external(_npts, _X, pY, pYp);
                static_cast<BesselSpline*>(s)->npts = _npts;
                static_cast<BesselSpline*>(s)->build();
                is_monotone[spl] = checkCubicSplineMonotonicity(_X, pY, pYp, _npts);
                break;

            case PCHIP_TYPE:
                s = new PchipSpline(h);
                static_cast<PchipSpline*>(s)->reserve_external(_npts, _X, pY, pYp);
                static_cast<PchipSpline*>(s)->npts = _npts;
                static_cast<PchipSpline*>(s)->build();
                is_monotone[spl] = checkCubicSplineMonotonicity(_X, pY, pYp, _npts);
                break;

            case HERMITE_TYPE:
                s = new HermiteSpline(h);
                static_cast<CubicSpline*>(s)->reserve_external(_npts, _X, pY, pYp);
                static_cast<CubicSpline*>(s)->npts = _npts;
                static_cast<CubicSpline*>(s)->build();
                is_monotone[spl] = checkCubicSplineMonotonicity(_X, pY, pYp, _npts);
                break;

            case QUINTIC_TYPE:
                s = new QuinticSpline(h);
                static_cast<QuinticSpline*>(s)->reserve_external(_npts, _X, pY, pYp, pYpp);
                static_cast<QuinticSpline*>(s)->npts = _npts;
                static_cast<QuinticSpline*>(s)->build();
                break;

            case SPLINE_SET_TYPE:
            case SPLINE_VEC_TYPE:
                //default:
                SPLINE_ASSERT(false,
                              "SplineSet::build\nAt spline n. " << spl << " named " << headers[spl] <<
                              "\n" << stype[size_t(spl)] <<
                              " not allowed as spline type\nin SplineSet::build for " << spl <<
                              "-th spline");
        }
        header_to_position[s->name()] = integer(spl);
    }

    baseValue   . must_be_empty("SplineSet::build, baseValue");
    basePointer . must_be_empty("SplineSet::build, basePointer");

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::getHeaders(vector<string> & h) const
{
    h.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        h[i] = splines[i]->name();
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval(real_type x, vector<real_type> & vals) const
{
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = (*splines[i])(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval(real_type x, real_type vals[], integer incy) const
{
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = (*splines[i])(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_D(real_type x, vector<real_type> & vals) const
{
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = splines[i]->D(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_D(real_type x, real_type vals[], integer incy) const
{
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = splines[i]->D(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_DD(real_type x, vector<real_type> & vals) const
{
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = splines[i]->DD(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_DD(real_type x, real_type vals[], integer incy) const
{
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = splines[i]->DD(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_DDD(real_type x, vector<real_type> & vals) const
{
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = splines[i]->DDD(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval_DDD(real_type x, real_type vals[], integer incy) const
{
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = splines[i]->DDD(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// vectorial values

Spline const *
SplineSet::intersect(integer     spl,
                     real_type   zeta,
                     real_type & x) const
{
    SPLINE_ASSERT(spl >= 0 && spl < _nspl,
                  "Spline n." << spl << " is not in SplineSet");
    SPLINE_ASSERT(is_monotone[size_t(spl)] > 0,
                  "Spline n." << spl <<
                  " is not monotone and can't be used as independent");
    Spline const * S = splines[size_t(spl)];
    // cerco intervallo intersezione
    real_type const * X = _Y[size_t(spl)];
    SPLINE_ASSERT(zeta >= X[0] && zeta <= X[size_t(_npts - 1)],
                  "SplineSet, evaluation at zeta = " << zeta <<
                  " is out of range: [" << X[0] <<
                  ", " << X[size_t(_npts - 1)] << "]");

    integer interval = integer(lower_bound(X, X + _npts, zeta) - X);
    if(interval > 0)
    {
        --interval;
    }
    if(isZero(X[size_t(interval)] - X[size_t(interval + 1)]))
    {
        ++interval;    // degenerate interval for duplicated nodes
    }
    if(interval >= _npts - 1)
    {
        interval = _npts - 2;
    }

    // compute intersection
    real_type a  = _X[size_t(interval)];
    real_type b  = _X[size_t(interval + 1)];
    real_type ya = X[size_t(interval)];
    real_type yb = X[size_t(interval + 1)];
    real_type DX = b - a;
    real_type DY = yb - ya;
    SPLINE_ASSERT(zeta >= ya && zeta <= yb,
                  "SplineSet, Bad interval [ " << ya << "," << yb << "] for zeta = " << zeta);
    SPLINE_ASSERT(a < b,
                  "SplineSet, Bad x interval [ " << a << "," << b << "]");
    if(S->type() == LINEAR_TYPE)
    {
        x = a + (b - a) * (zeta - ya) / (yb - ya);
    }
    else
    {
        real_type const * dX = _Yp[size_t(spl)];
        real_type        dya = dX[interval];
        real_type        dyb = dX[interval + 1];
        real_type coeffs[4] = { ya - zeta, dya, (3 * DY / DX - 2 * dya - dyb) / DX, (dyb + dya - 2 * DY / DX) / (DX * DX) };
        real_type real[3], imag[3];
        pair<int, int> icase = cubicRoots(coeffs, real, imag);
        SPLINE_ASSERT(icase.first > 0,
                      "SplineSet, No intersection found with independent spline at zeta = " << zeta);
        // cerca radice buona
        bool ok = false;
        for(integer i = 0; i < icase.first && !ok; ++i)
        {
            ok = real[i] >= 0 && real[i] <= DX;
            if(ok)
            {
                x = a + real[i];
            }
        }
        SPLINE_ASSERT(ok, "SplineSet, failed to find intersection with independent spline at zeta = " << zeta);
    }
    return S;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2(integer             spl,
                 real_type           zeta,
                 vector<real_type> & vals) const
{
    real_type x;
    intersect(spl, zeta, x);
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = (*splines[i])(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2(integer   spl,
                 real_type zeta,
                 real_type vals[],
                 integer   incy) const
{
    real_type x;
    intersect(spl, zeta, x);
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = (*splines[i])(x);
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_D(integer             spl,
                   real_type           zeta,
                   vector<real_type> & vals) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type ds = S->D(x);
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = splines[i]->D(x) / ds;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_D(integer   spl,
                   real_type zeta,
                   real_type vals[],
                   integer   incy) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type ds = S->D(x);
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = splines[i]->D(x) / ds;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_DD(integer             spl,
                    real_type           zeta,
                    vector<real_type> & vals) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type dt  = 1 / S->D(x);
    real_type dt2 = dt * dt;
    real_type ddt = -S->DD(x) * (dt * dt2);
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
    {
        vals[i] = splines[i]->DD(x) * dt2 + splines[i]->D(x) * ddt;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_DD(integer   spl,
                    real_type zeta,
                    real_type vals[],
                    integer   incy) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type dt  = 1 / S->D(x);
    real_type dt2 = dt * dt;
    real_type ddt = -S->DD(x) * (dt * dt2);
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
    {
        vals[ii] = splines[i]->DD(x) * dt2 + splines[i]->D(x) * ddt;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_DDD(integer             spl,
                     real_type           zeta,
                     vector<real_type> & vals) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type dt  = 1 / S->D(x);
    real_type dt3 = dt * dt * dt;
    real_type ddt = -S->DD(x) * dt3;
    real_type dddt = 3 * (ddt * ddt) / dt - S->DDD(x) * (dt * dt3);
    vals.resize(size_t(_nspl));
    for(size_t i = 0; i < size_t(_nspl); ++i)
        vals[i] = splines[i]->DDD(x) * dt3 +
                  3 * splines[i]->DD(x) * dt * ddt +
                  splines[i]->D(x) * dddt;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void
SplineSet::eval2_DDD(integer   spl,
                     real_type zeta,
                     real_type vals[],
                     integer   incy) const
{
    real_type x;
    Spline const * S = intersect(spl, zeta, x);
    real_type dt  = 1 / S->D(x);
    real_type dt3 = dt * dt * dt;
    real_type ddt = -S->DD(x) * dt3;
    real_type dddt = 3 * (ddt * ddt) / dt - S->DDD(x) * (dt * dt3);
    size_t ii = 0;
    for(size_t i = 0; i < size_t(_nspl); ++i, ii += size_t(incy))
        vals[ii] = splines[i]->DDD(x) * dt3 +
                   3 * splines[i]->DD(x) * dt * ddt +
                   splines[i]->D(x) * dddt;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2(real_type    zeta,
                 char const * indep,
                 char const * name) const
{
    vector<real_type> vals;
    eval2(getPosition(indep), zeta, vals);
    return vals[size_t(getPosition(name))];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_D(real_type    zeta,
                   char const * indep,
                   char const * name) const
{
    vector<real_type> vals;
    eval2_D(getPosition(indep), zeta, vals);
    return vals[size_t(getPosition(name))];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_DD(real_type    zeta,
                    char const * indep,
                    char const * name) const
{
    vector<real_type> vals;
    eval2_DD(getPosition(indep), zeta, vals);
    return vals[size_t(getPosition(name))];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_DDD(real_type    zeta,
                     char const * indep,
                     char const * name) const
{
    vector<real_type> vals;
    eval2_DDD(getPosition(indep), zeta, vals);
    return vals[size_t(getPosition(name))];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2(real_type zeta, integer indep, integer spl) const
{
    vector<real_type> vals;
    eval2(indep, zeta, vals);
    return vals[size_t(spl)];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_D(real_type zeta, integer indep, integer spl) const
{
    vector<real_type> vals;
    eval2_D(indep, zeta, vals);
    return vals[size_t(spl)];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_DD(real_type zeta, integer indep, integer spl) const
{
    vector<real_type> vals;
    eval2_DD(indep, zeta, vals);
    return vals[size_t(spl)];
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

real_type
SplineSet::eval2_DDD(real_type zeta, integer indep, integer spl) const
{
    vector<real_type> vals;
    eval2_DDD(indep, zeta, vals);
    return vals[size_t(spl)];
}
}
