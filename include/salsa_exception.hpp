//------------------------------------------------------------------------------
//
// Contributors: 
//             1) Tommy Hinks
//
//------------------------------------------------------------------------------

#ifndef SALSA_EXCEPTION_HPP_INCLUDED
#define SALSA_EXCEPTION_HPP_INCLUDED

#include "salsa_namespace.hpp"
#include <sstream>
#include <exception>

//------------------------------------------------------------------------------

// Useful macro.
#define SALSA_THROW(MSG) {       \
    ::std::stringstream ss;      \
    ss << "salsa: ";             \
    ss << MSG;                   \
    throw salsa::base(ss.str()); \
}                                

//------------------------------------------------------------------------------

BEGIN_SALSA_NAMESPACE

class base : public std::exception
{
public:

    //! CTOR.
    explicit
    base(const std::string &msg)
        : std::exception()
        , _msg(msg) 
    {}

    //! Copy CTOR.
    base(const base &rhs)
        : std::exception(rhs)
        , _msg(rhs._msg)
    {}

    //! Assign.
    base&
    operator=(const base &rhs)
    {
        std::exception::operator=(rhs);
        _msg = rhs._msg;
        return *this;
    }

    //! DTOR.
    virtual 
    ~base()
    {}

    virtual const char* 
    what() const
    { return _msg.c_str(); }

private:    // Member variables.

    std::string _msg;
};

END_SALSA_NAMESPACE

//------------------------------------------------------------------------------

#endif	// SALSA_EXCEPTION_HPP_INCLUDED
