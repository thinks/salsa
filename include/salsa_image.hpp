//------------------------------------------------------------------------------
//
// Contributors: 
//             1) Tommy Hinks
//
//------------------------------------------------------------------------------

#ifndef SALSA_IMAGE_HPP_INCLUDED
#define SALSA_IMAGE_HPP_INCLUDED

#include "salsa_namespace.hpp"
#include "salsa_exception.hpp"
#include <vector>
#include <iostream>
#include <iomanip>

//------------------------------------------------------------------------------

BEGIN_SALSA_NAMESPACE

//! Basic image, just a 2D grid of pixels.
template<typename P>
class image
{	
public:

    typedef P pixel_type;

    // CTOR.
    explicit
    image(const std::size_t  width, 
          const std::size_t  height, 
          const P           &val    = P())
        : _width(width)
        , _height(height)
        , _pixels(width*height, val) // May throw.
    {}

    // Default copy & assign.

    void
    clear(const P &val = P())
    {
        int i;
        const int n = static_cast<int>(_pixels.size());
#pragma omp parallel for
        for (i = 0; i < n; ++i) {
            _pixels[i] = val;
        }
    }

    std::size_t            
    width()	const 
    { return _width; }

    std::size_t            
    height() const 
    { return _height; }
    
    //! In [pixels]
    std::size_t 
    size() const 
    { return _pixels.size(); }
    
    //! Approximate.
    std::size_t
    mem_used() const
    { return (sizeof(image<P>) + _pixels.capacity()*sizeof(P)); }

    // Operator access.
    //
    const P& 
    operator[](const std::size_t i) const 
    { return _pixels[i]; }

    P& 
    operator[](const std::size_t i) 
    { return _pixels[i]; }

    const P& 
    operator()(const std::size_t x, const std::size_t y) const 
    { return _pixels[_linear_index(x, y)]; }

    P& 
    operator()(const std::size_t x, const std::size_t y) 
    { return _pixels[_linear_index(x, y)]; }

private:	// Member variables.

    std::size_t    _width;	//!< Width.
    std::size_t    _height;	//!< Height.
    std::vector<P> _pixels;	//!< Data.

private:	// Implementation helper functions.

    std::size_t
    _linear_index(const std::size_t x, const std::size_t y) const
    { return (x + y*_width); }
};	

END_SALSA_NAMESPACE

//------------------------------------------------------------------------------

BEGIN_STD_NAMESPACE

template<typename P>
ostream&
operator<<(ostream& os, salsa::image<P>& rhs)
{
    os	<< "salsa::image<P>[0x" << &rhs << "]:\n"
        << "sizeof(P) : " << sizeof(P)				 << " [bytes]\n"
        << "Width     : " << rhs.width()			 << " [pixels]\n"
        << "Height    : " << rhs.height()			 << " [pixels]\n"
        << fixed << setprecision(3)
        << "Mem Used  : " << 0.000001*rhs.mem_used() << " [MB]\n";
    return os;
}

END_STD_NAMESPACE

//------------------------------------------------------------------------------

#endif	// SALSA_IMAGE_HPP_INCLUDED
