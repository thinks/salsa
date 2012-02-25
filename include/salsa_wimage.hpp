//------------------------------------------------------------------------------
//
// Contributors: 
//             1) Tommy Hinks
//
//------------------------------------------------------------------------------

#ifndef SALSA_WIMAGE_HPP_INCLUDED
#define SALSA_WIMAGE_HPP_INCLUDED

#include "salsa_namespace.hpp"
#include "salsa_exception.hpp"
#include "salsa_image.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

//------------------------------------------------------------------------------

BEGIN_SALSA_NAMESPACE

template<typename P>
class wimage
{
public:

    typedef typename image<P> image_type;
    typedef typename image<P>::pixel_type pixel_type;

    explicit
    wimage(const double  wxmin, const double wxmax,
           const double  wymin, const double wymax, 
           const double  wdx,   const double wdy, 
           const P      &val = P())
        : _wxmin(std::min<double>(wxmin, wxmax))
        , _wxmax(std::max<double>(wxmin, wxmax))
        , _wdx(wdx)
        , _wxdim(wxmax - wxmin)
        , _wxinvdim(1./_wxdim)
        , _wymin(std::min<double>(wymin, wymax))
        , _wymax(std::max<double>(wymin, wymax))
        , _wdy(wdy)
        , _wydim(wymax - wymin)
        , _wyinvdim(1./_wydim)
        , _img(static_cast<std::size_t>(std::ceil(_wxdim/wdx)),
               static_cast<std::size_t>(std::ceil(_wydim/wdy)),
               val)
    {}

    // Default copy & assign.

    const image_type&
    image() const
    { return _img; }

    image_type&
    image()
    { return _img; }

    double 
    wxmin() const	
    { return _wxmin; }

    double 
    wxmax() const	
    { return _wxmax; }

    double 
    wdx() const	
    { return _wdx; }

    double 
    wxdim() const	
    { return _wxdim; }

    double 
    wymin() const	
    { return _wymin; }

    double 
    wymax() const	
    { return _wymax; }

    double 
    wdy() const	
    { return _wdy; }

    double 
    wydim() const	
    { return _wydim; }

    std::size_t
    mem_used() const
    { return (sizeof(wimage<P>) + _img.mem_used()); }

    const P& 
    operator[](const std::size_t i) const 
    { return _img[i]; }

    P&		 
    operator[](const std::size_t i)		  
    { return _img[i]; }

    const P& 
    operator()(const std::size_t x, const std::size_t y) const 
    { return _img(x, y); }

    P& 
    operator()(const std::size_t x, const std::size_t y)				
    { return _img(x, y); }

    const P& 
    operator()(const double wx, const double wy) const 
    { return _worldToPixel(wx, wy); }

    P& 
    operator()(const double wx, const double wy)				
    { return _worldToPixel(wx, wy); }

private:	// Member variables.

    double _wxmin;
    double _wxmax;
    double _wdx;
    double _wxdim;
    double _wxinvdim;
    double _wymin;
    double _wymax;
    double _wdy;
    double _wydim;
    double _wyinvdim;
    image_type _img;

private:


    double 
    _round(const double r) 
    { return (r > 0.0) ? std::floor(r + 0.5) : std::ceil(r - 0.5); }

    const P&
    _worldToPixel() const
    {
        return _img(
            static_cast<std::size_t>(
                _round((_img.width() - 1)*((wx - _wxmin)*_wxinvdim))),
            static_cast<std::size_t>(
                _round((_img.height() - 1)*((wy - _wymin)*_wyinvdim))));
    }

    P&
    _worldToPixel(const double wx, const double wy)
    {
        return _img(
            static_cast<std::size_t>(
                _round((_img.width() - 1)*((wx - _wxmin)*_wxinvdim))),
            static_cast<std::size_t>(
                _round((_img.height() - 1)*((wy - _wymin)*_wyinvdim))));
    }
};

//template <typename P>
//void 
//world2pixel(const imagew<P>&    wimg, 
//            const thx::vec2f64& w2, 
//            thx::vec2i32*       p)
//{
//    assert(0 != p);
//    assert(wimg.wmin().x <= w2.x && w2.x <= wimg.wmax().x);
//    assert(wimg.wmin().y <= w2.y && w2.y <= wimg.wmax().y);
//
//    p->x = thx::round_int32(
//            (wimg.img().width() - 1)*std::max<thx::float64>(
//                0.0, (w2.x - wimg.wmin().x)*wimg.inv_wdim().x));
//    p->y = thx::round_int32(
//            (wimg.img().height() - 1)*std::max<thx::float64>(
//                0.0, (w2.y - wimg.wmin().y)*wimg.inv_wdim().y));
//
//    assert(0 <= p->x && p->x < wimg.img().width());
//    assert(0 <= p->y && p->y < wimg.img().height());
//}

END_SALSA_NAMESPACE

//------------------------------------------------------------------------------

BEGIN_STD_NAMESPACE

template<typename P>
ostream&
operator<<(ostream& os, salsa::wimage<P>& rhs)
{
    os	<< "salsa::imagew<P>[0x" << &rhs << "]:\n"
        << "sizeof(P) : " << sizeof(P)				 << " [bytes]\n"
        << "Width     : " << rhs.image().width()	 << " [pixels]\n"
        << "Height    : " << rhs.image().height()	 << " [pixels]\n"
        << fixed << setprecision(3)
        << "W Min     : [" << rhs.wxmin() << ", " << rhs.wymin() << "] [m]\n"
        << "W Max     : [" << rhs.wxmax() << ", " << rhs.wymax() << "] [m]\n"
        << "W Dim     : [" << rhs.wxdim() << ", " << rhs.wydim() << "] [m]\n"
        << "W Dx      : [" << rhs.wdx()   << ", " << rhs.wdy()   << "] [m]\n"
        << "Mem Used  : " << 0.000001*rhs.mem_used() << " [MB]\n";
    return os;
}

END_STD_NAMESPACE

//------------------------------------------------------------------------------

#endif	// SALSA_WIMAGE_HPP_INCLUDED
