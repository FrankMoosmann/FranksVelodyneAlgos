/* CVM Class Library */
/* Copyright (C), Sergei Nikolaev, 1992-2006, http://cvmlib.com */

#ifndef _CVM_H
#define _CVM_H

#if defined (__INTEL_COMPILER)
#   pragma warning(disable:1744)
#   pragma warning(disable:1125)
#   pragma warning(disable:1195)
#   pragma warning(disable:383)
#   pragma warning(disable:810)
#   pragma warning(disable:981)
#   pragma warning(disable:1418)
#   pragma warning(disable:171)
#   pragma warning(disable:1684)
#   pragma warning(disable:1599)
#endif

#   if !defined (CVM_NO_MT)
#       define CVM_MT
#       if !defined (_PTHREADS)
#           define _PTHREADS
#       endif
#   endif

// MSVC++ 6.0 and higher settings
#if defined (_MSC_VER)
#   pragma once
#   define WIN32_LEAN_AND_MEAN        // Exclude rarely-used stuff from Windows headers
#   define _WIN32_WINNT 0x500       // at least Win2000 is required for InitializeCriticalSectionAndSpinCount
#   include <windows.h>
#   include <process.h>
#   include <time.h>
#   if (_MSC_VER < 1400)
#       error "Please use stable version 5.2 for older MSVC versions"
#   endif
#   if !defined (__INTEL_COMPILER) || !defined(_WIN64)
#       define CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
#   endif
#   pragma warning(disable:4290)
#   pragma warning(push)
#   pragma warning(disable:4250)
#   pragma warning(disable:4251)
#   pragma warning(disable:4311)
#   if defined (CVM_FLOAT)
#       pragma warning(disable:4244)
#   endif
#   if defined (SRC_EXPORTS) && !defined (CVM_EXPORTS)
#       define CVM_EXPORTS
#   endif
#   include <hash_map>
#   if defined (CVM_USES_STLPORT)
#       define CVM_BLOCKS_MAP std::hash_map
#   else
#       define CVM_BLOCKS_MAP stdext::hash_map
#   endif
#   ifdef CVM_STATIC
#       define CVM_API
#   else
#       ifdef CVM_EXPORTS
#           define CVM_API __declspec(dllexport)
#       else
#           define CVM_API __declspec(dllimport)
#       endif
#   endif
#   if defined (_MT) && !defined (CVM_NO_MT)
#       define CVM_MT
#   endif

#   include <limits>
    typedef __int64 CVM_LONGEST_INT;
    typedef unsigned __int64 CVM_LONGEST_UINT;

#   define CVM_VSNPRINTF vsnprintf_s
#   define CVM_VSNPRINTF_S_DEFINED
#   define CVM_STRCPY_S_DEFINED

// GCC settings
#elif defined (__GNUC__)
#   ifdef __MINGW32__               // Dev-C++ under Win32 assumed here
#       define WIN32_LEAN_AND_MEAN
#       include <windows.h>
#   else
#       include <semaphore.h>       // Unix
#   endif
#   ifdef __stdcall
#       undef __stdcall
#   endif
#   define __stdcall
#   define CVM_API
#   define CVM_STDEXT stdext
#   define CVM_BLOCKS_MAP std::map
typedef long long CVM_LONGEST_INT;
typedef unsigned long long CVM_LONGEST_UINT;
#define CVM_VSNPRINTF vsnprintf

#else
#   error "Unsupported compiler"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <cstring>
#include <limits>

// fix for missing __builtin_clog functions
#if defined (__INTEL_COMPILER)
#   if defined (_GLIBCXX_USE_C99_COMPLEX)
#       undef _GLIBCXX_USE_C99_COMPLEX
#       define _GLIBCXX_USE_C99_COMPLEX 0
#   endif
#endif

#include <complex>
#include <algorithm>
#include <exception>
#include <new>

#if defined (STLPORT)
#   define CVM_USES_STLPORT
#endif

#if defined (_DEBUG) || defined (DEBUG)
#   define CVM_DEBUG
#   include <assert.h>
#   define CVM_ASSERT(p,n) _cvm_assert(p,n);
#else
#   define CVM_ASSERT(p,n)
#endif

// errors codes
#define CVM_OK                           0
#define CVM_OUTOFMEMORY                  1
#define CVM_OUTOFRANGE                   2
#define CVM_OUTOFRANGE1                  3
#define CVM_OUTOFRANGE2                  4
#define CVM_WRONGSIZE                    5
#define CVM_SIZESMISMATCH                6
#define CVM_WRONGMKLARG                  7
#define CVM_WRONGMKLARG2                 8
#define CVM_SINGULARMATRIX               9
#define CVM_NOTPOSITIVEDEFINITE          10
#define CVM_WRONGCHOLESKYFACTOR          11
#define CVM_WRONGBUNCHKAUFMANFACTOR      12
#define CVM_NOTPOSITIVEDIAG              13
#define CVM_CONVERGENCE_ERROR            14
#define CVM_DIVISIONBYZERO               15
#define CVM_SEMAPHOREERROR               16
#define CVM_READ_ONLY_ACCESS             17
#define CVM_SUBMATRIXACCESSERROR         18
#define CVM_SUBMATRIXNOTAVAILABLE        19
#define CVM_MATRIXNOTSYMMETRIC           20
#define CVM_MATRIXNOTHERMITIAN           21
#define CVM_BREAKS_HERMITIANITY          22
#define CVM_METHODNOTAVAILABLE           23
#define CVM_NOTIMPLEMENTED               24

#define CVM_THE_LAST_ERROR_CODE          25                    // use it in derived classes
#define CVM_MATRIX_ELEMENT_SEPARATOR     " "
#define CVM_EOL                          std::endl

typedef unsigned char tbyte;                                   // memory allocation quantum

#ifdef CVM_NO_NAMESPACE
#   define CVM_NAMESPACE_BEG
#   define CVM_NAMESPACE_END
#else
#   define CVM_NAMESPACE_BEG namespace cvm {
#   define CVM_NAMESPACE_END }
#endif


CVM_NAMESPACE_BEG

template <typename T>               class basic_array;
template <typename TR, typename TC> class Array;
template <typename TR, typename TC> class Matrix;
template <typename TR, typename TC> class SqMatrix;
template <typename TR, typename TC> class BandMatrix;
template <typename TR>              class basic_rvector;
template <typename TR>              class basic_rmatrix;
template <typename TR>              class basic_srmatrix;
template <typename TR, typename TC> class basic_cvector;
template <typename TR, typename TC> class basic_cmatrix;
template <typename TR, typename TC> class basic_scmatrix;
template <typename TR>              class basic_srbmatrix;
template <typename TR, typename TC> class basic_scbmatrix;
template <typename TR>              class basic_srsmatrix;
template <typename TR, typename TC> class basic_schmatrix;
template <typename T,  typename TR> class type_proxy;

#if !defined (CVM_ALLOCATOR)
#   define CVM_ALLOCATOR std::allocator
#endif

class ErrMessages;
CVM_API ErrMessages& ErrMessagesInstance();

// error messages holder
class ErrMessages
{
    friend class cvmexception;
    // message string maps
    typedef std::map <int, std::string, std::less<int> > map_Msg;
    typedef map_Msg::iterator                            itr_Msg;
    typedef map_Msg::const_iterator                      citr_Msg;
    typedef std::pair <int, std::string>                 pair_Msg;

private:
    std::string msUnknown;
    map_Msg mmMsg;

    CVM_API ErrMessages();
    CVM_API const std::string& _get (int nException);
    CVM_API bool _add (int nNewCause, const char* szNewMessage);

    friend CVM_API ErrMessages& ErrMessagesInstance();

public:
    ~ErrMessages() {}
};

class cvmexception : public std::exception
{
protected:
    int mnCause;
    mutable char mszMsg[256];

public:
    cvmexception()
        : mnCause (CVM_OK)
    {
        mszMsg[0] = '\0';
    }

    explicit cvmexception (int nCause, ...)
        : mnCause (nCause)
    {
        va_list argList;
        va_start (argList, nCause);
#if defined (CVM_VSNPRINTF_S_DEFINED)
        const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg), sizeof(mszMsg) - 1, ErrMessagesInstance()._get(mnCause).c_str(), argList);
#else
        const int nLength = CVM_VSNPRINTF (mszMsg, sizeof(mszMsg) - 1, ErrMessagesInstance()._get(mnCause).c_str(), argList);
#endif
        va_end (argList);
        if (nLength >= (int) sizeof(mszMsg))
        {
            mszMsg[sizeof(mszMsg) - 1] = '\0';
        }
    }

    cvmexception (const cvmexception& e)
        : std::exception(e), mnCause (e.mnCause)
    {
#if defined (CVM_STRCPY_S_DEFINED)
        strcpy_s (mszMsg, sizeof(mszMsg), e.mszMsg);
#else
        strcpy (mszMsg, e.mszMsg);
#endif
    }

    virtual ~cvmexception() throw() {}

    int cause() const
    {
        return mnCause;
    }

    virtual const char* what() const throw()
    {
        return mszMsg;
    }

    static int getNextCause ()
    {
        return CVM_THE_LAST_ERROR_CODE;
    }

    static bool add (int nNewCause, const char* szNewMessage)
    {
        return ErrMessagesInstance()._add (nNewCause, szNewMessage);
    }
};

// exported utilities - were recently refactored from specialization because of Borland's understanding of ANSI C++
template <typename T>
CVM_API void __copy (int nSize, const T* pFrom, int nFromIncr, T* pTo, int nToIncr);
template <typename T>
CVM_API void __swap (int nSize, T* p1, int n1Incr, T* p2, int n2Incr);

template <typename TC, typename TR>
CVM_API TR _real (const TC& mT);
template <typename TC, typename TR>
CVM_API TR _imag (const TC& mT);

template<typename TR>
CVM_API TR __dot (const TR* mpD, int mnSize, int mnIncr, const TR* pD, int nIncr);
template<typename TC>
CVM_API TC __dotu (const TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr);
template<typename TC>
CVM_API TC __dotc (const TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr);

template <typename TR, typename TC>
CVM_API TR __norm (const TC*  pD, int nSize, int nIncr);

template <typename TC>
CVM_API int __idamax (const TC* pD, int nSize, int nIncr);
template <typename TC>
CVM_API int __idamin (const TC* pD, int nSize, int nIncr);

template <typename TC>
CVM_API void __add (TC* mpD, int mnSize, int mnIncr, const TC* pv, int nIncr);
template <typename TC>
CVM_API void __subtract (TC* mpD, int mnSize, int mnIncr, const TC* pv, int nIncr);
template <typename TR, typename TC>
CVM_API void __scal (TC* mpD, int mnSize, int mnIncr, TR dScal);
template<typename TC>
CVM_API void __conj (TC* mpD, int mnSize, int mnIncr);

template <typename TR, typename TC>
CVM_API void __copy_real (TC* mpD, int mnSize, int mnIncr, const TR* pRe, int nReIncr);
template <typename TR, typename TC>
CVM_API void __copy_imag (TC* mpD, int mnSize, int mnIncr, const TR* pRe, int nReIncr);
template <typename TR, typename TC>
CVM_API void __copy2 (TC* mpD, int mnSize, int mnIncr, const TR* pRe, const TR* pIm, int nReIncr = 1, int nImIncr = 1);

template <typename TC, typename TM, typename TV>
CVM_API void __gemv (bool bLeft, const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TC, typename TM, typename TV>
CVM_API void __gbmv (bool bLeft, const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TR, typename TM, typename TV>
CVM_API void __symv (const TM& m, TR dAlpha, const TV& v, TR dBeta, TV& vRes);
template <typename TC, typename TM, typename TV>
CVM_API void __shmv (const TM& m, TC dAlpha, const TV& v, TC dBeta, TV& vRes);
template <typename TC, typename TM>
CVM_API void __gemm (const TM& ml, bool bTrans1, const TM& mr, bool bTrans2, TC dAlpha, TM& mRes, TC dBeta);
template <typename TR, typename TSM, typename TM>
CVM_API void __symm (bool bLeft, const TSM& ml, const TM& mr, TR dAlpha, TM& mRes, TR dBeta);
template <typename TC, typename TSM, typename TM>
CVM_API void __hemm (bool bLeft, const TSM& ml, const TM& mr, TC dAlpha, TM& mRes, TC dBeta);

template <typename TC, typename TV>
CVM_API void __polynom (TC* mpD, int ldP, int mnM, const TC* pD, int ldA, const TV& v);
template <typename T>
CVM_API void __inv (T& m, const T& mArg) throw (cvmexception);
template <typename T, typename TR>
CVM_API void __exp (T& m, const T& mArg, TR tol) throw (cvmexception);
template <typename T, typename TR>
CVM_API void __exp_symm (T& m, const T& mArg, TR tol) throw (cvmexception);
template <typename TR, typename TC, typename TRM>
CVM_API void __solve (const TRM& m, int nrhs, const TC* pB, int ldB, TC* pX, int ldX, TR& dErr, const TC* pLU, const int* pPivots) throw (cvmexception);
template <typename TC, typename TM, typename TSM>
CVM_API void __svd (TC* pD, int nSize, int nIncr, const TM& mArg, TSM* mU, TSM* mVH) throw (cvmexception);
template <typename TR, typename TM, typename TX>
CVM_API void __pinv (TX& mX, const TM& mArg, TR threshold) throw (cvmexception);
template <typename TV, typename TSM, typename TSCM>
CVM_API void __eig (TV& vRes, const TSM& mArg, TSCM* mEigVect, bool bRightVect) throw (cvmexception);
template<typename TR, typename TM>
CVM_API void __cond_num (const TM& mArg, TR& dCondNum) throw (cvmexception);
template <typename TM>
void __low_up (TM& m, int* nPivots) throw (cvmexception);

template <typename TR>
CVM_API void __randomize (TR* mpD, int mnSize, int mnIncr, TR dFrom, TR dTo);
template <typename TC, typename TR>
CVM_API void __randomize_real (TC* mpD, int mnSize, int mnIncr, TR dFrom, TR dTo);
template <typename TC, typename TR>
CVM_API void __randomize_imag (TC* mpD, int mnSize, int mnIncr, TR dFrom, TR dTo);

template <typename TR, typename TM, typename TV>
CVM_API void __ger (TM& m, const TV& vCol, const TV& vRow, TR dAlpha);
template <typename TC, typename TM, typename TV>
CVM_API void __geru (TM& m, const TV& vCol, const TV& vRow, TC cAlpha);
template <typename TC, typename TM, typename TV>
CVM_API void __gerc (TM& m, const TV& vCol, const TV& vRow, TC cAlpha);

template <typename TC, typename TSM>
CVM_API void __syrk (bool bTransp, TC alpha, int k, const TC* pA, int ldA, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __syr2k (bool bTransp, TC alpha, int k, const TC* pA, int ldA, const TC* pB, int ldB, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __herk (bool bTransp, TC alpha, int k, const TC* pA, int ldA, TC beta, TSM& m);
template <typename TC, typename TSM>
CVM_API void __her2k (bool bTransp, TC alpha, int k, const TC* pA, int ldA, const TC* pB, int ldB, TC beta, TSM& m);

template <typename TM>
CVM_API int __cholesky (TM& m);
template <typename TM>
CVM_API void __bunch_kaufman (TM& m, int* nPivots) throw (cvmexception);
template<typename TR, typename TM, typename TV>
CVM_API void __poequ (const TM& m, TV& vScalings, TR& dCond, TR& dMax);


template <typename T> inline
const T& _cvm_min (const T& x, const T& y)
{
    return x < y ? x : y;
}

template <typename T> inline
const T& _cvm_max (const T& x, const T& y)
{
    return x > y ? x : y;
}

template<class TR>
inline const TR* __get_real_p (const std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_re;
#else
    return reinterpret_cast<const TR*>(c);
#endif
}

template<class TR>
inline const TR* __get_imag_p (const std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_im;
#else
    return reinterpret_cast<const TR*>(c) + 1;
#endif
}

template<class TR>
inline TR* __get_real_p (std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_re;
#else
    return reinterpret_cast<TR*>(c);
#endif
}

template<class TR>
inline TR* __get_imag_p (std::complex<TR>* c)
{
#if defined (CVM_USES_STLPORT)
    return &c->_M_im;
#else
    return reinterpret_cast<TR*>(c) + 1;
#endif
}

template<class TR> inline const TR* __get_real_p (const TR* d) { return d; }
template<class TR> inline const TR* __get_imag_p (const TR* d) { return d; }
template<class TR> inline       TR* __get_real_p (      TR* d) { return d; }
template<class TR> inline       TR* __get_imag_p (      TR* d) { return d; }

// characters for fortran subroutines
class Chars
{
protected:
    static char mchars[15];

public:
    static const char* pT () {return mchars;}
    static const char* pN () {return mchars + 1;}
    static const char* pU () {return mchars + 2;}
    static const char* pL () {return mchars + 3;}
    static const char* pP () {return mchars + 4;}
    static const char* pQ () {return mchars + 5;}
    static const char* pB () {return mchars + 6;}
    static const char* pE () {return mchars + 7;}
    static const char* pR () {return mchars + 8;}
    static const char* pA () {return mchars + 9;}
    static const char* pS () {return mchars + 10;}
    static const char* pV () {return mchars + 11;}
    static const char* pO () {return mchars + 12;}
    static const char* pI () {return mchars + 13;}
    static const char* pC () {return mchars + 14;}
};

template <class TR>
inline const TR& basic_cvmMachMin()
{
    static const TR _min = (std::numeric_limits<TR>::min)();
    return _min;
}

template <class TR>
inline const TR& basic_cvmMachSp()
{
    static const TR _eps = (std::numeric_limits<TR>::epsilon)();
    return _eps;
}


// class of memory blocks
class CVM_API MemoryBlocks
{
    struct BlockProperty
    {
        size_t mnSize;
        int mnRefCount;
        BlockProperty (size_t nSize, int nRefCount) : mnSize (nSize), mnRefCount (nRefCount) {}
    };

#ifdef CVM_USE_POOL_MANAGER
    typedef std::map<tbyte*, BlockProperty, std::less<tbyte*> > map_Blocks;     // pointer -> {size, refcount}
    typedef map_Blocks::iterator itr_Blocks;

    typedef std::multimap<size_t, tbyte*> map_FreeBs;                           // size -> pointer
    typedef map_FreeBs::iterator itr_FreeBs;

    typedef std::map<tbyte*, itr_FreeBs, std::less<tbyte*> > map_FreeIt;        // pointer -> iterator to FreeBs
    typedef map_FreeIt::iterator itr_FreeIt;

    map_FreeBs mFreeBs;                                                         // currently free blocks by sizes
    map_FreeIt mFreeIt;                                                         // currently free blocks iterators by pointers
#else
    typedef CVM_BLOCKS_MAP<CVM_LONGEST_UINT, BlockProperty> map_Blocks;         // pointer -> refcount
    typedef map_Blocks::iterator itr_Blocks;
#endif

    map_Blocks mBlocks;                                                         // currently occupied or freed blocks by pointers

public:
#ifdef CVM_USE_POOL_MANAGER
    void    AddBlock     (tbyte* pBlock, size_t nBytes, bool bOccupied);
    tbyte*  GetFreeBlock (size_t nBytes);
    void    AddPair      (tbyte* pBlock, size_t nBytes, size_t nRest);
#   ifdef CVM_DEBUG
    void    Assert       (const void* pvBlock, size_t nBytes);
#   endif
#else
    void AddNew (tbyte* pBlock, size_t nBytes)     // for just allocated only, i.e. non-const
    {
        mBlocks.insert (std::pair<CVM_LONGEST_UINT, BlockProperty>
                       (reinterpret_cast<CVM_LONGEST_UINT>(pBlock), BlockProperty(nBytes, 1)));
    }
#endif

    tbyte*  AddRef       (const tbyte* pBlock);
    int     FreeBlock    (tbyte* pBlock);
};


// memory pool class
class CVM_API MemoryPool
{
#ifdef CVM_USE_POOL_MANAGER
    typedef std::list<tbyte*> list_blocks;

    struct DeletePtr {
        template<class T>
        void operator () (T* p) const
        {
            ::delete[] p;
        }
    };

    list_blocks  mOutBlocks;                                                    // outer memory blocks
#endif

    MemoryBlocks mMemoryBlocks;                                                 // currently existing blocks and their statuses

public:
    MemoryPool();
   ~MemoryPool();

    tbyte* Malloc (size_t nBytes) throw (cvmexception);
    tbyte* AddRef (const tbyte* pD);                                            // increases a reference counter
    int    Free   (tbyte*& pToFree) throw (cvmexception);                       // decreases a reference counter and
                                                                                // returns memory back to the pool if the counter is zeroed
#ifdef CVM_USE_POOL_MANAGER
#   ifdef CVM_DEBUG
    void Assert (const void* pvBlock, size_t nBytes)                            // synchronized outside
    {
        mMemoryBlocks.Assert (pvBlock, nBytes);
    }
#   endif
    void Clear();                                                               // destroys all outer blocks in reverse order
#endif
};


CVM_API tbyte* _cvmMalloc  (size_t nBytes) throw (cvmexception);
CVM_API tbyte* _cvmAddRef  (const tbyte* pD);
CVM_API int    _cvmFree    (tbyte*& pD);
CVM_API void   _cvm_assert (const void* pvBlock, size_t nBytes);
CVM_API void   cvmExit();

template <typename T>
inline T* cvmMalloc (size_t nEls) throw (cvmexception)
{
    return reinterpret_cast<T*>(_cvmMalloc (nEls * sizeof (T)));
}
template <typename T>
inline T* cvmAddRef (const T* pD)
{
    return reinterpret_cast<T*>(_cvmAddRef (reinterpret_cast<const tbyte*>(pD)));
}
template <typename T>
inline int cvmFree (T*& pD)
{
    return _cvmFree (reinterpret_cast<tbyte*&>(pD));
}



// inline utilities
template <typename T>
inline void CleanMemory (T* p, int nEls)
{
    memset (p, 0, nEls * sizeof(T));
}

inline float _abs (const float& v)
{
    return static_cast<float>(fabs (v));
}
inline double _abs (const double& v)
{
    return fabs(v);
}
inline float _abs (const std::complex<float>& v)
{
    return static_cast<float>(sqrt (v.real() * v.real() + v.imag() * v.imag()));
}
inline double _abs (const std::complex<double>& v)
{
    return sqrt (v.real() * v.real() + v.imag() * v.imag());
}
inline int _abs (const int& v)
{
    return abs(v);
}
inline float _sqrt (const float& v)
{
    return static_cast<float>(sqrt (v));
}
inline double _sqrt (const double& v)
{
    return sqrt(v);
}
inline std::complex<float> _conjugate (const std::complex<float>& v)
{
    return std::complex<float> (v.real(), - v.imag());
}
inline std::complex<double> _conjugate (const std::complex<double>& v)
{
    return std::complex<double> (v.real(), - v.imag());
}
template <typename TR>
inline bool _conjugated (const std::complex<TR>& v1, const std::complex<TR>& v2)
{
    return _abs (v1.real() - v2.real()) <= basic_cvmMachMin<TR>() &&
           _abs (v1.imag() + v2.imag()) <= basic_cvmMachMin<TR>();
}

template <typename T, typename TR>
inline std::ostream& operator << (std::ostream& os, const type_proxy<T,TR>& mOut)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    os.imbue (std::locale::empty());
#endif
    os << static_cast<T>(mOut);
    return os;
}

template <typename T, typename TR>
inline std::istream& operator >> (std::istream& is, type_proxy<T,TR>& v)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    is.imbue (std::locale::empty());
#endif
    T t;
    is >> t;
    v = t;
    return is;
}

template <typename T>
inline std::istream& operator >> (std::istream& is, basic_array<T>& aIn)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    is.imbue (std::locale::empty());
#endif
    CVM_ASSERT(aIn.mpD, aIn.mnSize * sizeof(T))
    for (int i = 0; i < aIn.mnSize && is.good(); ++i)
    {
        is >> aIn.mpD[i];
    }
    return is;
}

template <typename T>
inline std::ostream& operator << (std::ostream& os, const basic_array<T>& aOut)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    os.imbue (std::locale::empty());
#endif
    CVM_ASSERT(aOut.mpD, aOut.mnSize * sizeof(T))
    for (int i = 0; i < aOut.mnSize && os.good(); ++i)
    {
        os << aOut.mpD[i] << CVM_MATRIX_ELEMENT_SEPARATOR;
    }
    os << CVM_EOL;
    return os;
}

template <typename TR, typename TC>
inline std::istream& operator >> (std::istream& is, Array<TR,TC>& aIn)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    is.imbue (std::locale::empty());
#endif
    const int nSize = aIn.mnSize * aIn.mnIncr;
    CVM_ASSERT(aIn.mpD, ((aIn.mnSize - 1) * aIn.mnIncr + 1) * sizeof(TC))
    for (int i = 0; i < nSize && is.good(); i += aIn.mnIncr)
    {
        is >> aIn.mpD[i];
    }
    return is;
}

template <typename TR, typename TC>
inline std::ostream& operator << (std::ostream& os, const Array<TR,TC>& aOut)
{
#if !defined (CVM_USES_STLPORT) && defined (_MSC_VER)
    os.imbue (std::locale::empty());
#endif
    const int nSize = aOut.mnSize * aOut.mnIncr;
    CVM_ASSERT(aOut.mpD, ((aOut.mnSize - 1) * aOut.mnIncr + 1) * sizeof(TC))

    for (int i = 0; i < nSize && os.good(); i += aOut.mnIncr)
    {
        os << aOut.mpD[i] << CVM_MATRIX_ELEMENT_SEPARATOR;
    }
    os << CVM_EOL;
    return os;
}

template <typename TR, typename TC>
inline std::ostream& operator << (std::ostream& os, const Matrix<TR,TC>& mOut)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    os.imbue (std::locale::empty());
#endif
    for (int i = 0; i < mOut.mnM; ++i)
    {
        for (int j = 0; j < mOut.mnN && os.good(); ++j)
        {
            os << mOut._ij_val (i, j) << CVM_MATRIX_ELEMENT_SEPARATOR;
        }
        os << CVM_EOL;
    }
    return os;
}

template <typename TR, typename TC>
inline std::istream& operator >> (std::istream& is, Matrix<TR,TC>& mIn)
{
#if defined (_MSC_VER) && !defined (CVM_USES_STLPORT)
    is.imbue (std::locale::empty());
#endif
    TC v;
    for (int i = 0; i < mIn.mnM; ++i)
    {
        for (int j = 0; j < mIn.mnN && is.good(); ++j)
        {
            is >> v;
            mIn._ij_proxy_val (i, j) = v;
        }
    }
    return is;
}

template <typename TR, typename TC, typename RM, typename RBM>
inline void _copy_b_matrix (RM& m, RBM& mb, bool bLeftToRight)
{
    const int nM   = mb.msize();
    const int nN   = mb.nsize();
    const int nKL  = mb.lsize();
    const int nKU  = mb.usize();
    const int nCol = 1 + nKL + nKU;
    int nS, nShiftL, nShiftR;
    TC* pL;
    TC* pR;

    for (int i = 0; i < nN; ++i)
    {
        nS = nCol;
        nShiftL = 0;
        nShiftR = 0;
        if (i < nKU)
        {
            nShiftR = nKU - i;
            nS -= nShiftR;
        }
        else
        {
            nShiftL = i - nKU;
        }
        if (nN - i <= nKL)
        {
            nS -= nKL + 1 - (nN - i);
        }

        pL = m.get()  + i * nM + nShiftL;
        pR = mb.get() + i * nCol + nShiftR;

        __copy<TC> (nS,
                    bLeftToRight ? pL : pR,
                    1,
                    bLeftToRight ? pR : pL,
                    1);
    }
}


template <typename TR, typename TC>
inline void _sum (TC* pD, int nSize, int nIncr,
                  const TC* p1, int nIncr1, const TC* p2, int nIncr2)           // pD = a1 + a2
{
    if (pD == p1)
    {
        if (pD == p2)
        {
            static const TR two(2.);
            __scal<TR, TC> (pD, nSize, nIncr, two);
        }
        else
        {
            __add<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
    else
    {
        if (pD == p2)
        {
            __add<TC> (pD, nSize, nIncr, p1, nIncr1);
        }
        else
        {
            __copy<TC> (nSize, p1, nIncr1, pD, nIncr);
            __add<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
}

template <typename TR, typename TC>
inline void _diff (TC* pD, int nSize, int nIncr, 
                   const TC* p1, int nIncr1, const TC* p2, int nIncr2)          // this = v1 - v2
{
    if (pD == p1)
    {
        if (pD == p2)
        {
            static const TR zero(0.);
            __scal<TR, TC> (pD, nSize, nIncr, zero);
        }
        else
        {
            __subtract<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
    else
    {
        if (pD == p2)
        {
            static const TR mone(-1.);
            __subtract<TC> (pD, nSize, nIncr, p1, nIncr1);
            __scal<TR, TC> (pD, nSize, nIncr, mone);
        }
        else
        {
            __copy<TC> (nSize, p1, nIncr1, pD, nIncr);
            __subtract<TC> (pD, nSize, nIncr, p2, nIncr2);
        }
    }
}

template <typename TR, typename TC>
inline void _incr (TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr)    // mpD = mpD + pD
{
    if (mpD == pD)
    {
        static const TR two(2.);
        __scal<TR, TC> (mpD, mnSize, mnIncr, two);
    }
    else
    {
        __add<TC> (mpD, mnSize, mnIncr, pD, nIncr);
    }
}

template <typename TR, typename TC>
inline void _decr (TC* mpD, int mnSize, int mnIncr, const TC* pD, int nIncr)    // mpD = mpD - pD
{
    if (mpD == pD)
    {
        static const TR zero(0.);
        __scal<TR, TC> (mpD, mnSize, mnIncr, zero);
    }
    else
    {
        __subtract<TC> (mpD, mnSize, mnIncr, pD, nIncr);
    }
}

template <typename TR, typename TC>
void _set_real (TC* mpD, int mnSize, int mnIncr, TR d)                          // fills real part
{
    const int nIncr2 = 2 * mnIncr;
    const int nSize = mnSize * nIncr2;
    TR* pD = __get_real_p<TR>(mpD);

    for (int i = 0; i < nSize; i += nIncr2)
    {
        CVM_ASSERT(pD, (i + 1) * sizeof(TR))
        pD[i] = d;
    }
}

template <typename TR, typename TC>
void _set_imag (TC* mpD, int mnSize, int mnIncr, TR d)                          // fills imaginary part
{
    const int nIncr2 = 2 * mnIncr;
    const int nSize = mnSize * nIncr2;
    TR* pD = __get_imag_p<TR>(mpD);

    for (int i = 0; i < nSize; i += nIncr2)
    {
        CVM_ASSERT(pD, (i + 1) * sizeof(TR))
        pD[i] = d;
    }
}


// this class provides read-write differentiation for a particular type
template <typename T, typename TR>
class type_proxy
{
    typedef type_proxy<T,TR> P;
    
protected:
    T&   mT;
    bool mbReadOnly;

public:
    type_proxy (T& ref, bool read_only) : mT(ref), mbReadOnly(read_only)
    {
    }

    type_proxy (const T& ref, bool read_only = true) :
            mT(const_cast<T&>(ref)), mbReadOnly(read_only)                      // read only by definition
    {
    }

    type_proxy (const type_proxy& p) : mT(p.mT), mbReadOnly(p.mbReadOnly)
    {
    }

    operator T () const
    {
        return mT;
    }
/*
    doesn't work - masks the operator above
    operator const T& () const
    {
        return mT;
    }

    operator T& ()
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return mT;
    }
*/
    T val() const
    {
        return mT;
    }

    T& get() throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return mT;
    }
    const T& get() const
    {
        return mT;
    }

    type_proxy& operator = (const type_proxy& p)
    {
        mT = p.mT;
        mbReadOnly = p.mbReadOnly;
        return *this;
    }
    type_proxy& operator = (const T& v) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT = v;
        return *this;
    }

    T* operator & () throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        return &mT;
    }
    const T* operator & () const
    {
        return &mT;
    }

    template <typename U>
    T& operator += (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT += T(u);
        return mT;
    }
    template <typename U>
    T& operator -= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT -= T(u);
        return mT;
    }
    template <typename U>
    T& operator *= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT *= T(u);
        return mT;
    }
    template <typename U>
    T& operator /= (const U& u) throw (cvmexception)
    {
        if (mbReadOnly) throw cvmexception (CVM_READ_ONLY_ACCESS);
        mT /= T(u);
        return mT;
    }

    template <typename U>
    T operator + (const U& u) const
    {
        return mT + T(u);
    }
    template <typename U>
    T operator - (const U& u) const
    {
        return mT - T(u);
    }
    template <typename U>
    T operator * (const U& u) const
    {
        return mT * T(u);
    }
    template <typename U>
    T operator / (const U& u) const
    {
        return mT / T(u);
    }

    T operator + (const P& u) const
    {
        return mT + u.mT;
    }
    T operator - (const P& u) const
    {
        return mT - u.mT;
    }
    T operator * (const P& u) const
    {
        return mT * u.mT;
    }
    T operator / (const P& u) const
    {
        return mT / u.mT;
    }

    T operator - () const
    {
        return - mT;
    }

    // specialized for std::complex<treal> only. link error would be received in other case
    TR real() const
    {
        return _real<T,TR>(mT);
    }
    TR imag() const
    {
        return _imag<T,TR>(mT);
    }
};


template <typename T, typename TR>
inline std::complex<TR> operator + (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u + p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator - (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u - p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator * (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u * p.val();
}
template <typename T, typename TR>
inline std::complex<TR> operator / (std::complex<TR> u, const type_proxy<T,TR>& p)
{
    return u / p.val();
}

template <typename U, typename T, typename TR>
inline T operator + (U u, const type_proxy<T,TR>& p)
{
    return T(u) + p.val();
}
template <typename U, typename T, typename TR>
inline T operator - (U u, const type_proxy<T,TR>& p)
{
    return T(u) - p.val();
}
template <typename U, typename T, typename TR>
inline T operator * (U u, const type_proxy<T,TR>& p)
{
    return T(u) * p.val();
}
template <typename U, typename T, typename TR>
inline T operator / (U u, const type_proxy<T,TR>& p)
{
    return T(u) / p.val();
}


// random numbers helper class
template <typename T>
class Randomizer
{
    T mdUMax;

    Randomizer ()
        : mdUMax (static_cast<T>(RAND_MAX))
    {
        srand (static_cast<unsigned int>(time(NULL)));
    }

    T _get (T dFrom, T dTo)
    {
        const T dMin = _cvm_min<T>(dFrom, dTo);
        const T dMax = _cvm_max<T>(dFrom, dTo);
        unsigned int nr = static_cast<unsigned int>(rand());
        return dMin + static_cast<T>(nr) * (dMax - dMin) / mdUMax;
    }

public:
    ~Randomizer () {}

    static T get (T dFrom, T dTo)
    {
        static Randomizer r;
        return r._get (dFrom, dTo);
    }
};


template <class T>
inline CVM_ALLOCATOR<T>& AllocatorInstance()
{
    static CVM_ALLOCATOR<T> _A;
    return _A;
}

template <typename T>
class basic_array
{
protected:
    int mnSize;                                                                 // number of array elements allocated
    T*  mpD;                                                                    // data pointer

public:
    // STL type definitions
    typedef T value_type;
    typedef value_type* pointer;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    basic_array()
        : mnSize(0), mpD(NULL)
    {
    }

    explicit basic_array (int nSize, bool bZeroMemory = true)
        : mnSize(nSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        if (bZeroMemory) CleanMemory<T> (mpD, mnSize);
    }

    basic_array (const T* p, int nSize)
        : mnSize(nSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, p, 1, mpD, 1);
    }

    basic_array (const T* first, const T* last)
        : mnSize(int(last - first)), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, first, 1, mpD, 1);
    }

    basic_array (const basic_array& a)
        : mnSize(a.mnSize), mpD (cvmMalloc<T>(size_t(mnSize)))
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        __copy<T> (mnSize, a.mpD, 1, mpD, 1);
    }

    virtual ~basic_array()
    {
        cvmFree<T>(mpD);
    }

    int size() const
    {
        return mnSize;
    }

    T* get()
    {
        return mpD;
    }

    const T* get() const
    {
        return mpD;
    }

    operator T* ()
    {
        return mpD;
    }

    operator const T* () const
    {
        return mpD;
    }

    // 1-based indexing operators
    T& operator () (int nFI) throw (cvmexception)                               // element access (returns l-value)
    {
        return this -> at(size_type(nFI - 1));
    }

    T operator () (int nFI) const throw (cvmexception)                          // element access (does not return l-value)
    {
        return this -> at(size_type(nFI - 1));
    }

    T& operator [] (size_type n) throw (cvmexception)                           // element access (returns l-value)
    {
        return this -> at(n - 1);
    }
    T operator [] (size_type n) const throw (cvmexception)                      // element access (does not return l-value)
    {
        return this -> at(n - 1);
    }
    T& operator [] (int n) throw (cvmexception)                                 // element access (returns l-value)
    {
        return this -> at(size_type(n - 1));
    }
    T operator [] (int n) const throw (cvmexception)                            // element access (does not return l-value)
    {
        return this -> at(size_type(n - 1));
    }

    basic_array& operator = (const basic_array& a) throw (cvmexception)
    {
        if (mnSize != a.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (a.mpD, 1);
        return *this;
    }

    basic_array& assign (const T* p)
    {
        this -> _assign (p, 1);
        return *this;
    }

    basic_array& set (T x)
    {
        this -> _set (x);                                                       // fills the content
        return *this;
    }

    basic_array& resize (int nNewSize) throw (cvmexception)
    {
        this -> _resize (nNewSize);
        return *this;
    }

    // std stuff
    iterator       begin()        {return this -> mpD;}
    const_iterator begin() const  {return this -> mpD;}
    iterator       end()          {return this -> mpD + mnSize;}
    const_iterator end()   const  {return this -> mpD + mnSize;}

    reverse_iterator rbegin()             {return reverse_iterator(end());}
    const_reverse_iterator rbegin() const {return const_reverse_iterator(end());}
    reverse_iterator rend()               {return reverse_iterator(begin());}
    const_reverse_iterator rend()   const {return const_reverse_iterator(begin());}

    size_type max_size() const    {return size_type(-1) / sizeof(T);}
    size_type capacity() const    {return size_type(mnSize);}
    bool empty() const            {return mnSize > 0;}

    reference front()             {return *begin();}
    const_reference front() const {return *begin();}
    reference back()              {return *(end() - 1);}
    const_reference back()  const {return *(end() - 1);}

    void reserve (size_type n) throw (cvmexception)
    {
        this -> _resize (int(n));
    }

    void assign (size_type n, const T& val) throw (cvmexception)
    {
        if (n > mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < n; ++i)
        {
            mpD[i] = val;
        }
    }

    void assign (const_iterator first, const_iterator last) throw (cvmexception)
    {
        const int n = last - first;
        if (n > mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < n; ++i)
        {
            mpD[i] = *(first + i);
        }
    }

    void resize (size_type nNewSize) throw (cvmexception)
    {
        this -> _resize (int(nNewSize));
    }

    void clear()
    {
        this -> _resize (0);
    }

    void swap (basic_array& v) throw (cvmexception)
    {
        if (mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        __swap<T> (mnSize, mpD, 1, v.mpD, 1);
    }

    // 0-based
    reference at (size_type n) throw (cvmexception)
    {
        return this -> _at(n);
    }

    const_reference at (size_type n) const throw (cvmexception)
    {
        return this -> _at(n);
    }

    // very very slow. provided for compatibility only
    void push_back (const T& x) throw (cvmexception)
    {
        this -> _resize (mnSize + 1);
        mpD[mnSize - 1] = x;
    }

    // very very slow. provided for compatibility only
    void pop_back () throw (cvmexception)
    {
        if (mnSize > 0) this -> _resize (mnSize - 1);
    }

    // very very slow. provided for compatibility only
    iterator insert (iterator position, const T& x) throw (cvmexception)
    {
        const int n = int (position - this -> begin());
        if (n > mnSize || n < 0) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        this -> _resize (mnSize + 1);
        for (int i = mnSize - 1; i > n; --i)
        {
            mpD[i] = mpD[i-1];
        }
        mpD[n] = x;
        return iterator (mpD + n);
    }

    // very very slow. provided for compatibility only
    iterator erase (iterator position) throw (cvmexception)
    {
        const int n = int (position - this -> begin());
        if (n > mnSize || n < 0) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = n; i < mnSize - 1; ++i)
        {
            mpD[i] = mpD[i+1];
        }
        if (mnSize > 0) this -> _resize (mnSize - 1);
        return iterator (mpD + n);
    }
    // end of stl stuff

protected:
    // 0-based
    virtual T& _at (size_type n) throw (cvmexception)
    {
        if (n >= size_type(mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, (n + 1) * sizeof(T))
        return mpD [n];
    }

    // 0-based
    virtual const T& _at (size_type n) const throw (cvmexception)
    {
        if (n >= size_type(mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(mpD, (n + 1) * sizeof(T))
        return mpD [n];
    }

    virtual void _assign (const T* p, int)
    {
        if (mpD != p)
        {
            memcpy (mpD, p, mnSize * sizeof(T));
        }
    }

    virtual void _set (T d)                                                     // fills the content
    {
        CVM_ASSERT(mpD, mnSize * sizeof(T))
        for (int i = 0; i < mnSize; ++i)
        {
            mpD[i] = d;
        }
    }

    virtual void _resize (int nNewSize) throw (cvmexception)
    {
        if (nNewSize < 0) throw cvmexception (CVM_WRONGSIZE, nNewSize);
        if (nNewSize == 0)                                                      // just let's free object memory in this case
        {
            cvmFree<T>(mpD);
            mnSize = 0;
        }
        else if (nNewSize != mnSize)
        {
            T* pD = cvmMalloc<T>(nNewSize);
            if (nNewSize > mnSize) CleanMemory<T> (pD, nNewSize);
            const int nMinSize = _cvm_min<int>(nNewSize, mnSize);

            if (nMinSize > 0)
            {
                __copy<T> (nMinSize, mpD, 1, pD, 1);
            }
            cvmFree<T>(mpD);
            mpD = pD;
            mnSize = nNewSize;
            CVM_ASSERT(mpD, mnSize * sizeof(T))
        }
    }

    friend std::istream& operator >> <> (std::istream& is, basic_array& aIn);
    friend std::ostream& operator << <> (std::ostream& os, const basic_array& aOut);
};

// abstract array of numbers allocatable in pool
template <typename TR, typename TC>
class Array : public basic_array<TC>
{
    typedef size_t size_type;
    typedef basic_array<TC> BasicArray;

protected:
    int mnIncr;                                                                 // distance between array members (default is 1)

    Array()
        : mnIncr (0)
    {
    }

    explicit Array (int nSize, bool bZeroMemory = true)
        : BasicArray(nSize, bZeroMemory), mnIncr (1)
    {
    }

    Array (TC* pD, int nSize, int nIncr)
        : mnIncr (nIncr)
    {
        if (nSize <= 0) throw cvmexception (CVM_WRONGSIZE, nSize);
        this -> mnSize = nSize;
        this -> mpD = cvmAddRef<TC> (pD);
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
    }

public:
    int incr() const
    {
        return mnIncr;
    }

    int indofmax() const                                                        // index of max. module, undefined for matrices
    {
        return this -> _indofmax();
    }

    int indofmin() const                                                        // index of max. module, undefined for matrices
    {
        return this -> _indofmin();
    }

    virtual TR norm() const                                                     // Euclid norm
    {
        return __norm<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual TR norminf() const
    {
        CVM_ASSERT(this -> mpD, ((this -> _indofmax() - 1) * this -> mnIncr + 1) * sizeof(TC))
        return _abs (this -> mpD[(this -> _indofmax() - 1) * this -> mnIncr]);
    }

    virtual TR norm1() const
    {
        TR dNorm(0.);
        const int nSize = this -> mnSize * this -> mnIncr;
        for (int i = 0; i < nSize; i += this -> mnIncr)
        {
            dNorm += _abs(this -> mpD[i]);
        }
        return dNorm;
    }

    virtual TR norm2() const
    {
        return this -> norm();
    }

    const int* _pincr() const
    {
        return &this -> mnIncr;
    }

    const int* _psize() const
    {
        return &this -> mnSize;
    }

    // to be redefined in classes with non-traditional storage, like band matrices etc.
    virtual TC* _pd()
    {
        return this -> mpD;
    }

    // to be redefined in classes with non-traditional storage, like band matrices etc.
    virtual const TC* _pd() const
    {
        return this -> mpD;
    }

    void _div (TR d) throw (cvmexception)                                       // this = this / d for real only
    {
        static const TR one(1.);
        if (_abs(d) <= basic_cvmMachMin<TR>()) throw cvmexception (CVM_DIVISIONBYZERO);
        this -> _scal (one / d);
    }

protected:
    // 0-based
    virtual TC& _at (size_type n) throw (cvmexception)
    {
        if (n >= size_type(this -> mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(this -> mpD, (n + 1) * sizeof(TC))
        return this -> mpD [n * this -> mnIncr];
    }

    // 0-based
    virtual const TC& _at (size_type n) const throw (cvmexception)
    {
        if (n >= size_type(this -> mnSize)) throw cvmexception (CVM_OUTOFRANGE, n);
        CVM_ASSERT(this -> mpD, (n + 1) * sizeof(TC))
        return this -> mpD [n * this -> mnIncr];
    }

    bool _equals (const Array& a) const                                         // compares array elements
    {
        bool bRes = false;
        if (this -> mnSize == a.size())
        {
            bRes = true;
            if (this -> mpD != a.get())
            {
                for (int i = 0; i < this -> mnSize; ++i)
                {
                    if (_abs (this -> mpD[i * this -> mnIncr] - a.mpD[i * a.mnIncr]) > basic_cvmMachMin<TR>())
                    {
                        bRes = false;
                        break;
                    }
                }
            }
        }
        return bRes;
    }

    void _normalize()                                                           // array normalizing
    {
        const TR dNorm = this -> norm();
        if (dNorm > basic_cvmMachMin<TR>())
        {
            static const TR one(1.);
            this -> _scal (one / dNorm);
        }
    }

    void _replace (const Array& a) throw (cvmexception)                         // this = a
    {
        cvmFree<TC> (this -> mpD);
        this -> mpD = cvmMalloc<TC>(a.mnSize);
        this -> mnSize = a.mnSize;
        this -> mnIncr = 1;
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
    }

// virtual methods
    virtual void _scal (TR d)
    {
        __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    virtual int _indofmax() const                                               // index of max. module, undefined for matrices
    {
        return __idamax<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual int _indofmin() const                                               // index of max. module, undefined for matrices
    {
        return __idamin<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual void _set (TC d)                                                    // fills the content
    {
        CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
        const int nSize = this -> mnSize * this -> mnIncr;
        for (int i = 0; i < nSize; i += this -> mnIncr)
        {
            this -> mpD[i] = d;
        }
    }

    virtual void _assign (const TC* pD, int nIncr)
    {
        if (this -> mpD != pD)
        {
            __copy<TC> (this -> mnSize, pD, nIncr, this -> mpD, this -> mnIncr);
        }
    }

    virtual void _assign_shifted (TC* pDshifted, const TC* pD, int nSize, int nIncr, int)
    {
        if (pDshifted != pD)
        {
            __copy<TC> (nSize, pD, nIncr, pDshifted, this -> mnIncr);
        }
    }

    virtual void _resize (int nNewSize) throw (cvmexception)
    {
        if (nNewSize < 0) throw cvmexception (CVM_WRONGSIZE, nNewSize);
        if (nNewSize == 0)                                                      // just let's free object memory in this case
        {
            cvmFree<TC>(this -> mpD);
            this -> mnSize = this -> mnIncr = 0;
        }
        else if (nNewSize != this -> mnSize)
        {
            TC* pD = cvmMalloc<TC>(nNewSize);
            if (nNewSize > this -> mnSize) CleanMemory<TC> (pD, nNewSize);
            const int nMinSize = _cvm_min<int>(nNewSize, this -> mnSize);

            if (nMinSize > 0)
            {
                __copy<TC> (nMinSize, this -> mpD, this -> mnIncr, pD, 1);
            }
            cvmFree<TC>(this -> mpD);
            this -> mpD    = pD;
            this -> mnSize = nNewSize;
            this -> mnIncr = 1;
            CVM_ASSERT(this -> mpD, ((this -> mnSize - 1) * this -> mnIncr + 1) * sizeof(TC))
        }
    }

    friend std::istream& operator >> <> (std::istream& is, Array& aIn);
    friend std::ostream& operator << <> (std::ostream& os, const Array& aOut);
};

// vector of real numbers
template <typename TR>
class basic_rvector : public Array<TR,TR>
{
    typedef std::complex<TR> TC;
    typedef Array<TR,TR> BaseArray;

public:
    basic_rvector()
    {
    }

    explicit basic_rvector (int nSize)
        : BaseArray (nSize)
    {
    }

    basic_rvector (int nSize, TR d)
        : BaseArray (nSize, false)
    {
        this -> _set (d);
    }

    // WARNING!
    // The following constructor does not allocate memory!
    // It just shares a memory allocated before.
    // It is intented to make possible the following syntax:
    //
    // basic_rmatrix m (10, 20);
    // basic_rvector v (20);
    // ...
    // m[1] = v;            // assigns v to the 1st row of m
    //
    //
    // And for example this code...
    //
    // basic_rmatrix m (10,20);
    // basic_rvector vRow = m[1];
    //
    // ...will also call THIS constructor, and memory will be shared!

    // If you need the code like this with memory allocation, use the following:
    //
    // basic_rmatrix m (10,20);
    // basic_rvector vRow (m.msize());
    // vRow = m[1];
    //
    basic_rvector (TR* pD, int nSize, int nIncr = 1)
        : BaseArray (pD, nSize, nIncr)
    {
    }

    basic_rvector (const basic_rvector& v)
        : BaseArray (v.size(), false)
    {
        __copy<TR> (this -> mnSize, v, v.incr(), this -> mpD, this -> mnIncr);
    }

    basic_rvector& operator = (const basic_rvector& v) throw (cvmexception)     // assignment (equal sizes)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (v, v.incr());
        return *this;
    }

    basic_rvector& assign (const TR* pD, int nIncr = 1)                         // foreign array assignment
    {
        this -> _assign (pD, nIncr);
        return *this;
    }

    basic_rvector& assign (int n, const basic_rvector& v) throw (cvmexception)  // subvector assignment
    {
        if (n <= 0 || n > this -> mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        --n;
        if (v.mnSize + n > this -> mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> mpD + this -> mnIncr * n, v, v.size(), v.incr(), 0);
        return *this;
    }

    basic_rvector& set (TR d)                                                   // sets the content
    {
        this -> _set (d);
        return *this;
    }

    basic_rvector& resize (int nNewSize) throw (cvmexception)
    {
        this -> _resize (nNewSize);
        return *this;
    }

    bool operator == (const basic_rvector& v) const
    {
        return this -> _equals (v);
    }

    bool operator != (const basic_rvector& v) const
    {
        return !(this -> operator == (v));
    }

    // vector replacement
    basic_rvector& operator << (const basic_rvector& v) throw (cvmexception)
    {
        this -> _replace (v);
        __copy<TR> (this -> mnSize, v.mpD, v.mnIncr, this -> mpD, this -> mnIncr);
        return *this;
    }

    // v1 + v2
    basic_rvector operator + (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vSum (*this);
        __add<TR>(vSum.mpD, vSum.mnSize, vSum.mnIncr, v._pd(), v.incr());
        return vSum;
    }

    // v1 - v2
    basic_rvector operator - (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vDiff (*this);
        __subtract<TR>(vDiff.mpD, vDiff.mnSize, vDiff.mnIncr, v._pd(), v.incr());
        return vDiff;
    }

    // this = v1 + v2
    basic_rvector& sum (const basic_rvector& v1, const basic_rvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _sum<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    // this = v1 + v2
    basic_rvector& diff (const basic_rvector& v1, const basic_rvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _diff<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    basic_rvector& operator += (const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _incr<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_rvector& operator -= (const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _decr<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_rvector operator - () const throw (cvmexception)
    {
        static const TR mone(-1.);
        basic_rvector vRes (*this);
        vRes._scal (mone);
        return vRes;
    }

    basic_rvector operator * (TR dMult) const throw (cvmexception)
    {
        basic_rvector vRes (*this);
        vRes._scal (dMult);
        return vRes;
    }

    basic_rvector operator / (TR dDiv) const throw (cvmexception)
    {
        basic_rvector vRes (*this);
        vRes._div (dDiv);
        return vRes;
    }

    // this = this * d
    basic_rvector& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    // this = this / d
    basic_rvector& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_rvector& normalize()
    {
        this -> _normalize();
        return *this;
    }

    TR operator * (const basic_rvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        return __dot<TR>(this -> mpD, this -> mnSize, this -> mnIncr, v.mpD, v.mnIncr);
    }

    basic_rvector operator * (const basic_rmatrix<TR>& m) const throw (cvmexception)
    {
        if (this -> mnSize != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rvector vRes (m.nsize());
        m._multiply (vRes, *this, true);
        return vRes;
    }

    // this = v * m
    basic_rvector& mult (const basic_rvector& v, const basic_rmatrix<TR>& m) throw (cvmexception)
    {
        if (this -> mnSize != m.nsize() || v.size() != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, true);
        return *this;
    }

    // this = m * v
    basic_rvector& mult (const basic_rmatrix<TR>& m, const basic_rvector& v) throw (cvmexception)
    {
        if (this -> mnSize != m.msize() || v.size() != m.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, false);
        return *this;
    }

    // v_col * v_row (rank-1 update)
    basic_rmatrix<TR> rank1update (const basic_rvector& v) const
    {
        basic_rmatrix<TR> mRes (this -> mnSize, v.mnSize);
        static const TR one(1.);
        __ger<TR, basic_rmatrix<TR>, basic_rvector>(mRes, *this, v, one);
        return mRes;
    }

    // linear solvers
    basic_rvector& solve (const basic_srmatrix<TR>& mA, const basic_rvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> mnSize || mA.msize() != vB.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_rvector& solve (const basic_srmatrix<TR>& mA, const basic_rvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, vB, dErr);
    }

    basic_rvector& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> size() || mA.msize() != vB.size() || mA.msize() != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_rvector& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots, 
                             const basic_rvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, vB, dErr);
    }

    // singular value decomposition
    basic_rvector& svd (const basic_rmatrix<TR>& m) throw (cvmexception)
    {
        m._svd (*this, NULL, NULL);
        return *this;
    }

    basic_rvector& svd (const basic_cmatrix<TR,TC>& m) throw (cvmexception)
    {
        m._svd (*this, NULL, NULL);
        return *this;
    }

    basic_rvector& svd (const basic_rmatrix<TR>& m, basic_srmatrix<TR>& mU,
                                                    basic_srmatrix<TR>& mVH) throw (cvmexception)
    {
        m._svd (*this, &mU, &mVH);
        return *this;
    }

    basic_rvector& svd (const basic_cmatrix<TR,TC>& m, basic_scmatrix<TR,TC>& mU,
                                                       basic_scmatrix<TR,TC>& mVH) throw (cvmexception)
    {
        m._svd (*this, &mU, &mVH);
        return *this;
    }

    // eigenvalues
    // we don't use _eig here since this is the special case - symmetric matrix
    basic_rvector& eig (const basic_srsmatrix<TR>& m) throw (cvmexception)
    {
        __eig<basic_rvector, basic_srsmatrix<TR>, basic_srmatrix<TR> > (*this, m, NULL, true);
        return *this;
    }

    basic_rvector& eig (const basic_srsmatrix<TR>& m, basic_srmatrix<TR>& mEigVect) throw (cvmexception)
    {
        __eig<basic_rvector, basic_srsmatrix<TR>, basic_srmatrix<TR> > (*this, m, &mEigVect, true);
        return *this;
    }

    // we don't use _eig here since this is the special case - hermitian matrix
    basic_rvector& eig (const basic_schmatrix<TR,TC>& m) throw (cvmexception)
    {
        __eig<basic_rvector, basic_schmatrix<TR,TC>, basic_scmatrix<TR,TC> > (*this, m, NULL, true);
        return *this;
    }

    basic_rvector& eig (const basic_schmatrix<TR,TC>& m, basic_scmatrix<TR,TC>& mEigVect) throw (cvmexception)
    {
        __eig<basic_rvector, basic_schmatrix<TR,TC>, basic_scmatrix<TR,TC> > (*this, m, &mEigVect, true);
        return *this;
    }

    // ?gemv routines perform a matrix-vector operation defined as
    // this = alpha*m*v + beta * this,
    basic_rvector& gemv (bool bLeft, const basic_rmatrix<TR>& m, TR dAlpha, const basic_rvector& v, TR dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gemv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }

    basic_rvector& gbmv (bool bLeft, const basic_srbmatrix<TR>& m, TR dAlpha, const basic_rvector& v, TR dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gbmv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }

    basic_rvector& randomize (TR dFrom, TR dTo)
    {
        __randomize<TR>(this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        return *this;
    }
};


// vector of complex numbers
template <typename TR, typename TC>
class basic_cvector : public Array<TR,TC>
{
    typedef Array<TR,TC> BaseArray;

public:
    basic_cvector()
    {
    }

    explicit basic_cvector (int nSize)
        : BaseArray (nSize)
    {
    }

    basic_cvector (int nSize, TC c)
        : BaseArray (nSize, false)
    {
        this -> _set (c);
    }

    // see warning for similar constructor of basic_rvector
    basic_cvector (TC* pD, int nSize, int nIncr = 1)
        : BaseArray (pD, nSize, nIncr)
    {
    }

    basic_cvector (const basic_cvector& v)
        : BaseArray(v.size(), false)
    {
        __copy<TC> (this -> mnSize, v, v.incr(), this -> mpD, this -> mnIncr);
    }

    basic_cvector (const TR* pRe, const TR* pIm, int nSize, int nIncrRe = 1, int nIncrIm = 1)
        : BaseArray (nSize, pRe == NULL || pIm == NULL)
    {
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, pRe, pIm, nIncrRe, nIncrIm);
    }

    basic_cvector (const basic_rvector<TR>& vRe, const basic_rvector<TR>& vIm)
        : BaseArray (vRe.size(), false)
    {
        if (vRe.size() != vIm.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, vRe, vIm, vRe.incr(), vIm.incr());
    }

    basic_cvector (const TR* pA, int nSize, bool bRealPart = true, int nIncr = 1)
        : BaseArray (nSize)
    {
        if (bRealPart)
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, pA, NULL, nIncr, 0);
        else
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, NULL, pA, 0, nIncr);
    }

    explicit basic_cvector (const basic_rvector<TR>& v, bool bRealPart = true)
        : BaseArray (v.size())
    {
        if (bRealPart)
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, v, NULL, v.incr(), 0);
        else
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, NULL, v, 0, v.incr());
    }

    // returns vector (real part) which CAN be l-value
    basic_rvector<TR> real()
    {
        return basic_rvector<TR>(__get_real_p<TR>(this -> mpD), this -> mnSize, this -> mnIncr * 2);
    }

    // returns vector (imaginary part) which CAN be l-value
    basic_rvector<TR> imag()
    {
        return basic_rvector<TR>(__get_imag_p<TR>(this -> mpD), this -> mnSize, this -> mnIncr * 2);
    }

    basic_cvector& operator = (const basic_cvector& v) throw (cvmexception)     // assignment (equal sizes)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (v, v.incr());
        return *this;
    }

    // foreign array assignment
    basic_cvector& assign (const TC* pD, int nIncr = 1)
    {
        this -> _assign (pD, nIncr);
        return *this;
    }

    basic_cvector& assign (int n, const basic_cvector& v) throw (cvmexception)  // subvector assignment
    {
        if (n <= 0 || n > this -> mnSize) throw cvmexception (CVM_OUTOFRANGE, n);
        --n;
        if (v.mnSize + n > this -> mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> mpD + this -> mnIncr * n, v, v.size(), v.incr(), 0);
        return *this;
    }

    // fills the content
    basic_cvector& set (TC c)
    {
        this -> _set (c);
        return *this;
    }

    // assigns real array
    basic_cvector& assign_real (const basic_rvector<TR>& vRe) throw (cvmexception)
    {
        if (this -> mnSize != vRe.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, vRe, vRe.incr());
        return *this;
    }

    // assigns imaginary array
    basic_cvector& assign_imag (const basic_rvector<TR>& vIm) throw (cvmexception)
    {
        if (this -> mnSize != vIm.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_imag<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, vIm, vIm.incr());
        return *this;
    }

    // fills real part
    basic_cvector& set_real (TR d)
    {
        this -> _set_real_number (d);
        return *this;
    }

    // fills imaginary part
    basic_cvector& set_imag (TR d)
    {
        this -> _set_imag_number (d);
        return *this;
    }

    basic_cvector& resize (int nNewSize) throw (cvmexception)
    {
        this -> _resize (nNewSize);
        return *this;
    }

    bool operator == (const basic_cvector& v) const
    {
        return this -> _equals (v);
    }

    bool operator != (const basic_cvector& v) const
    {
        return !(this -> operator == (v));
    }

    // vector replacement
    basic_cvector& operator << (const basic_cvector& v) throw (cvmexception)
    {
        this -> _replace (v);
        __copy<TC> (this -> mnSize, v.mpD, v.mnIncr, this -> mpD, this -> mnIncr);
        return *this;
    }

    basic_cvector operator + (const basic_cvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cvector vSum (*this);
        __add<TC>(vSum.mpD, vSum.mnSize, vSum.mnIncr, v._pd(), v.incr());
        return vSum;
    }

    basic_cvector operator - (const basic_cvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cvector vDiff (*this);
        __subtract<TC>(vDiff.mpD, vDiff.mnSize, vDiff.mnIncr, v._pd(), v.incr());
        return vDiff;
    }

    // this = v1 + v2
    basic_cvector& sum (const basic_cvector& v1, const basic_cvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _sum<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    // this = v1 + v2
    basic_cvector& diff (const basic_cvector& v1, const basic_cvector& v2) throw (cvmexception)
    {
        if (this -> mnSize != v1.mnSize || this -> mnSize != v2.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _diff<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, v1, v1.incr(), v2, v2.incr());
        return *this;
    }

    basic_cvector& operator += (const basic_cvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _incr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_cvector& operator -= (const basic_cvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        _decr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, v, v.incr());
        return *this;
    }

    basic_cvector operator - () const throw (cvmexception)
    {
        static const TR mone(-1.);
        basic_cvector vRes (*this);
        vRes._scal(mone);
        return vRes;
    }

    basic_cvector operator * (TR dMult) const
    {
        basic_cvector vRes (*this);
        vRes._scal (dMult);
        return vRes;
    }

    basic_cvector operator / (TR dDiv) const throw (cvmexception)
    {
        basic_cvector vRes (*this);
        vRes._div (dDiv);
        return vRes;
    }

    basic_cvector operator * (TC cMult) const
    {
        basic_cvector vRes (*this);
        vRes._scal (cMult);
        return vRes;
    }

    basic_cvector operator / (TC cDiv) const throw (cvmexception)
    {
        basic_cvector vRes (*this);
        vRes._div (cDiv);
        return vRes;
    }

    basic_cvector& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_cvector& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_cvector& operator *= (TC cMult)
    {
        this -> _scal (cMult);
        return *this;
    }

    basic_cvector& operator /= (TC cDiv)
    {
        this -> _div (cDiv);
        return *this;
    }

    basic_cvector& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // conjugated vector
    basic_cvector operator ~ () const throw (cvmexception)
    {
        basic_cvector vRes (*this);
        __conj<TC> (vRes.mpD, vRes.mnSize, vRes.mnIncr);
        return vRes;
    }

    basic_cvector& conj (const basic_cvector& v) throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign (v, v.incr());
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    basic_cvector& conj()
    {
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    TC operator * (const basic_cvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        return __dotu<TC>(this -> mpD, this -> mnSize, this -> mnIncr, v.mpD, v.mnIncr);
    }

    TC operator % (const basic_cvector& v) const throw (cvmexception)
    {
        if (this -> mnSize != v.mnSize) throw cvmexception (CVM_SIZESMISMATCH);
        return __dotc<TC>(this -> mpD, this -> mnSize, this -> mnIncr, v.mpD, v.mnIncr);
    }

    basic_cvector operator * (const basic_cmatrix<TR,TC>& m) const throw (cvmexception)
    {
        if (this -> mnSize != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cvector vRes (m.nsize());
        m._multiply (vRes, *this, true);
        return vRes;
    }

    // this = v * m
    basic_cvector& mult (const basic_cvector& v, const basic_cmatrix<TR,TC>& m) throw (cvmexception)
    {
        if (this -> mnSize != m.nsize() || v.size() != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, true);
        return *this;
    }

    // this = m * v
    basic_cvector& mult (const basic_cmatrix<TR,TC>& m, const basic_cvector& v) throw (cvmexception)
    {
        if (this -> mnSize != m.msize() || v.size() != m.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        m._multiply (*this, v, false);
        return *this;
    }

    // v_col * v_row (rank-1 update, unconjugated)
    basic_cmatrix<TR,TC> rank1update_u (const basic_cvector& v) const
    {
        basic_cmatrix<TR,TC> mRes (this -> mnSize, v.mnSize);
        static const TC one(1.);
        __geru<TC, basic_cmatrix<TR,TC>, basic_cvector>(mRes, *this, v, one);
        return mRes;
    }

    // v_col * v_row (rank-1 update, conjugated)
    basic_cmatrix<TR,TC> rank1update_c (const basic_cvector& v) const
    {
        basic_cmatrix<TR,TC> mRes (this -> mnSize, v.mnSize);
        static const TC one(1.);
        __gerc<TC, basic_cmatrix<TR,TC>, basic_cvector>(mRes, *this, v, one);
        return mRes;
    }

    // linear solvers
    basic_cvector& solve (const basic_scmatrix<TR,TC>& mA, const basic_cvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> size() || mA.msize() != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_cvector& solve (const basic_scmatrix<TR,TC>& mA, const basic_cvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, vB, dErr);
    }

    basic_cvector& solve_lu (const basic_scmatrix<TR,TC>& mA, const basic_scmatrix<TR,TC>& mLU, const int* pPivots,
                             const basic_cvector& vB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != this -> size() || mA.msize() != vB.size() || mA.msize() != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (vB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_cvector& solve_lu (const basic_scmatrix<TR,TC>& mA, const basic_scmatrix<TR,TC>& mLU, const int* pPivots,
                             const basic_cvector& vB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, vB, dErr);
    }

    // eigenvaluess
    basic_cvector& eig (const basic_srmatrix<TR>& m) throw (cvmexception)
    {
        m._eig(*this, NULL, true);
        return *this;
    }

    basic_cvector& eig (const basic_srmatrix<TR>& m, basic_scmatrix<TR,TC>& mEigVect, bool bRightVect = true) throw (cvmexception)
    {
        m._eig(*this, &mEigVect, bRightVect);
        return *this;
    }

    basic_cvector& eig (const basic_scmatrix<TR,TC>& m) throw (cvmexception)
    {
        m._eig(*this, NULL, true);
        return *this;
    }

    basic_cvector& eig (const basic_scmatrix<TR,TC>& m, basic_scmatrix<TR,TC>& mEigVect, bool bRightVect = true) throw (cvmexception)
    {
        m._eig(*this, &mEigVect, bRightVect);
        return *this;
    }

    // ?gemv routines perform a matrix-vector operation defined as
    // this := alpha*m*v + beta * this,
    basic_cvector& gemv (bool bLeft, const basic_cmatrix<TR,TC>& m, TC dAlpha, const basic_cvector& v, TC dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gemv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }

    basic_cvector& gbmv (bool bLeft, const basic_scbmatrix<TR,TC>& m, TC dAlpha, const basic_cvector& v, TC dBeta) throw (cvmexception)
    {
        if ((bLeft ? m.msize() != v.mnSize : m.nsize() != v.mnSize) ||
            (bLeft ? m.nsize() != this -> mnSize : m.msize() != this -> mnSize)) throw cvmexception (CVM_SIZESMISMATCH);
        m._gbmv (bLeft, dAlpha, v, dBeta, *this);
        return *this;
    }

    basic_cvector& randomize_real (TR dFrom, TR dTo)
    {
        __randomize_real<TC,TR> (this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        return *this;
    }

    basic_cvector& randomize_imag (TR dFrom, TR dTo)
    {
        __randomize_imag<TC,TR> (this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        return *this;
    }

protected:
    void _div (TC d) throw (cvmexception)
    {
        if (_abs(d) <= basic_cvmMachMin<TR>()) throw cvmexception (CVM_DIVISIONBYZERO);
        static const TC one(1.);
        this -> _scal (one / d);
    }

    void _scal (TC d)
    {
        __scal<TC, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    void _set_real_number (TR d)
    {
        _set_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    void _set_imag_number (TR d)
    {
        _set_imag<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }
};

// generalized matrix
template <typename TR, typename TC>
class Matrix : public Array<TR,TC>
{
    typedef Array<TR,TC> BaseArray;

protected:
    int mnM;                                                                    // number of rows
    int mnN;                                                                    // number of columns
    int mnLD;                                                                   // leading dimension

    Matrix()
        : mnM(0), mnN(0), mnLD(0)
    {
    }

    Matrix (int nM, int nN, int nLD, bool bZeroMemory)
        : BaseArray (nLD * nN, bZeroMemory), mnM(nM), mnN(nN), mnLD(nLD)
    {
        if (mnM  <= 0) throw cvmexception (CVM_WRONGSIZE, mnM);
        if (mnN  <= 0) throw cvmexception (CVM_WRONGSIZE, mnN);
        if (mnLD <= 0) throw cvmexception (CVM_WRONGSIZE, mnLD);
    }

    Matrix (TC* pD, int nM, int nN, int nLD, int nSize)         // for submatrices
        : BaseArray (pD, nSize, 1), mnM(nM), mnN(nN), mnLD(nLD)
    {
    }

    Matrix (const BaseArray& v, bool bBeColumn)                 // true = column, false = row
        : BaseArray (v.size()), mnM (bBeColumn ? v.size() : 1), mnN (bBeColumn ? 1 : v.size()), mnLD(mnM)
    {
        __copy<TC> (this -> mnSize, v, v.incr(), this -> mpD, this -> mnIncr);
    }

public:
    int msize() const
    {
        return mnM;
    }

    int nsize() const
    {
        return mnN;
    }

    int ld() const
    {
        return mnLD;
    }

    int rowofmax() const
    {
        return (this -> _indofmax() - 1) % mnM + 1;
    }

    int rowofmin() const
    {
        return (this -> _indofmin() - 1) % mnM + 1;
    }

    int colofmax() const
    {
        return (this -> _indofmax() - 1) / mnM + 1;
    }

    int colofmin() const
    {
        return (this -> _indofmin() - 1) / mnM + 1;
    }

    virtual TR norm1() const
    {
        int i, j, k;
        TR  rSum, rNorm(0.);

        for (j = 0; j < mnN; ++j)
        {
            rSum = TR(0.);

            k = j * mnLD;
            for (i = 0; i < mnM; ++i)
            {
                CVM_ASSERT(this -> mpD, (k + i + 1) * sizeof(TC))
                rSum += _abs (this -> mpD [k + i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    virtual TR norminf() const
    {
        int i, j;
        TR  rSum, rNorm(0.);

        for (i = 0; i < mnM; ++i)
        {
            rSum = TR(0.);

            for (j = 0; j < mnN; ++j)
            {
                CVM_ASSERT(this -> mpD, (j * this -> mnLD + i + 1) * sizeof(TC))
                rSum += _abs (this -> mpD [j * this -> mnLD + i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    const int* _pm() const
    {
        return &mnM;
    }

    const int* _pn() const
    {
        return &mnN;
    }

    const int* _pld() const
    {
        return &mnLD;
    }

    TC* _sub_pointer_nocheck (int row, int col)
    {
        return this -> _pd() + ((col - 1) * this -> ld() + row - 1);
    }

    TC* _sub_pointer (int row, int col, int height, int width) throw (cvmexception)
    {
        if (row    <= 0) throw cvmexception (CVM_WRONGSIZE, row);
        if (col    <= 0) throw cvmexception (CVM_WRONGSIZE, col);
        if (height <= 0) throw cvmexception (CVM_WRONGSIZE, height);
        if (width  <= 0) throw cvmexception (CVM_WRONGSIZE, width);
        if (row + height - 1 > mnM || col + width - 1 > mnN) throw cvmexception (CVM_SIZESMISMATCH);
        return _sub_pointer_nocheck (row, col);
    }

    virtual int _ldm() const
    {
        return this -> ld();
    }

    virtual const int* _pldm() const
    {
        return this -> _pld();
    }

    virtual bool _continuous () const
    {
        return mnM == mnLD;
    }

    void _check_ld() const
    {
        if (!this -> _continuous())
        {
            throw cvmexception (CVM_SUBMATRIXACCESSERROR);
        }
    }

    virtual void _scal (TR d)
    {
        if (this -> _continuous())
        {
            __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __scal<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, d);
        }
    }

protected:
    virtual int _indofmax() const
    {
        this -> _check_ld();
        return __idamax<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual int _indofmin() const
    {
        this -> _check_ld();
        return __idamin<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
    }

    virtual void _assign (const TC* pD, int)
    {
        if (this -> mpD != pD)
        {
            if (this -> _continuous())
            {
                __copy<TC> (this -> mnSize, pD, 1, this -> mpD, this -> mnIncr);
            }
            else for (int i = 0; i < this -> mnN; ++i)
            {
                __copy<TC> (this -> mnM, pD + this -> mnM * i, 1, this -> mpD + this -> mnLD * i, this -> mnIncr);
            }
        }
    }

    virtual void _assign_shifted (TC* pDshifted, const TC* pD, int nRows, int nCols, int nLD)  // reusing nSise and nIncr parameter
    {
        if (pDshifted != pD)
        {
            for (int i = 0; i < nCols; ++i)
            {
                __copy<TC> (nRows, pD + nLD * i, 1, pDshifted + this -> mnLD * i, this -> mnIncr);
            }
        }
    }

    virtual void _set (TC d)
    {
        CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))
        int i, j, k;
        for (j = 0; j < this -> mnN; ++j)
        {
            k = j * mnLD;
            for (i = 0; i < this -> mnM; ++i)
            {
                this -> mpD[k + i] = d;
            }
        }
    }

    virtual void _massign (const Matrix& m)
    {
        if (this -> mpD != m.mpD)
        {
            if (this -> _continuous() && m._continuous())
            {
                __copy<TC> (this -> mnSize, m._pd(), m.incr(), this -> mpD, this -> mnIncr);
            }
            else
            {
                const TC* p = m._pd();
                const int nLD = m._ldm();
                for (int i = 0; i < this -> mnN; ++i)
                {
                    __copy<TC> (this -> mnM, p + nLD * i, m.incr(), this -> mpD + this -> mnLD * i, this -> mnIncr);
                }
            }
        }
    }

    virtual void _resize (int nNewM, int nNewN) throw (cvmexception)
    {
        if (nNewM != mnM || nNewN != mnN)
        {
            this -> _check_ld();

            if (nNewM < 0) throw cvmexception (CVM_WRONGSIZE, nNewM);
            if (nNewN < 0) throw cvmexception (CVM_WRONGSIZE, nNewN);
            const int nNewSize = nNewM * nNewN;

            if (nNewSize == 0)
            {
                cvmFree<TC>(this -> mpD);
                this -> mnSize = this -> mnIncr = this -> mnM = this -> mnN = 0;
            }
            else
            {
                TC* pD = cvmMalloc<TC>(nNewSize);
                if (nNewSize > this -> mnSize) CleanMemory<TC> (pD, nNewSize);
                const int nMinM = _cvm_min<int>(nNewM, this -> mnM);
                const int nMinN = _cvm_min<int>(nNewN, this -> mnN);

                for (int i = 0; i < nMinN; ++i)
                {
                    __copy<TC> (nMinM, this -> mpD + i * this -> mnM, this -> mnIncr, pD + i * nNewM, 1);
                }
                cvmFree<TC>(this -> mpD);
                this -> mpD    = pD;
                this -> mnSize = nNewSize;
                CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))

                this -> mnM    = nNewM;
                this -> mnN    = nNewN;
                this -> mnLD   = nNewM;
                this -> mnIncr = 1;
            }
        }
    }

    virtual int _ld_for_replace () const
    {
        return this -> mnLD;
    }

    virtual int _size_for_replace () const
    {
        return this -> mnSize;
    }

    void _replace (const Matrix& m) throw (cvmexception)
    {
        this -> _check_ld();                                    // submatrix replacement is obviously not possible
        cvmFree<TC>(this -> mpD);
        this -> mnSize = m._size_for_replace();
        this -> mpD = cvmMalloc<TC>(this -> mnSize);
        this -> mnIncr = 1;
        CVM_ASSERT(this -> mpD, (this -> mnSize * sizeof(TC)))
        this -> mnM  = m.mnM;
        this -> mnN  = m.mnN;
        this -> mnLD = m._ld_for_replace();
    }

    void _transp (const Matrix& m)
    {
        int i;
        if (this -> mnM > this -> mnN) for (i = 0; i < this -> mnN; ++i)
        {
            __copy<TC> (m.nsize(), m.get() + i, m.ld(), this -> mpD + i * this -> mnLD, 1);
        }
        else for (i = 0; i < this -> mnM; ++i)
        {
            __copy<TC> (m.msize(), m.get() + i * m.ld(), 1, this -> mpD + i, this -> mnLD);
        }
    }

    virtual type_proxy<TC,TR> _ij_proxy_val (int i, int j)                      // zero based
    {
        CVM_ASSERT(this -> mpD, (this -> mnLD * j + i + 1) * sizeof(TC))
        return type_proxy<TC,TR>(this -> mpD [this -> mnLD * j + i], false);
    }

    virtual TC _ij_val (int i, int j) const                                     // zero based
    {
        CVM_ASSERT(this -> mpD, (this -> mnLD * j + 1) * sizeof(TC))
        return this -> mpD [this -> mnLD * j + i];
    }

    virtual void _swap_rows (int n1, int n2) throw (cvmexception)
    {
        if (n1 <= 0 || n1 > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, n1);
        if (n2 <= 0 || n2 > this -> mnM) throw cvmexception (CVM_OUTOFRANGE2, n2);
        if (n1 != n2)
        {
            __swap<TC> (this -> mnN, this -> mpD + n1 - 1, this -> mnLD, this -> mpD + n2 - 1, this -> mnLD);
        }
    }

    virtual void _swap_cols (int n1, int n2) throw (cvmexception)
    {
        if (n1 <= 0 || n1 > mnN) throw cvmexception (CVM_OUTOFRANGE1, n1);
        if (n2 <= 0 || n2 > mnN) throw cvmexception (CVM_OUTOFRANGE2, n2);
        if (n1 != n2)
        {
            __swap<TC> (this -> mnM, this -> mpD + (n1 - 1) * this -> mnLD, 1, this -> mpD + (n2 - 1) * this -> mnLD, 1);
        }
    }

    virtual const TC* _p (const Matrix& m) const
    {
        return m._pd();
    }

    virtual void _msum (const Matrix& m1, const Matrix& m2)
    {
        if (this -> _continuous() && m1._continuous() && m2._continuous())
        {
            _sum<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m1), m1.incr(), this -> _p(m2), m2.incr());
        }
        else for (int i = 0; i < mnN; ++i)
        {
            _sum<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m1) + m1._ldm() * i, m1.incr(), this -> _p(m2) + m2._ldm() * i, m2.incr());
        }
    }

    void _mdiff (const Matrix& m1, const Matrix& m2)
    {
        if (this -> _continuous() && m1._continuous() && m2._continuous())
        {
            _diff<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m1), m1.incr(), this -> _p(m2), m2.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _diff<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m1) + m1._ldm() * i, m1.incr(), this -> _p(m2) + m2._ldm() * i, m2.incr());
        }
    }

    virtual void _mincr (const Matrix& m)
    {
        if (this -> _continuous() && m._continuous())
        {
            _incr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m), m.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _incr<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m) + m._ldm() * i, m.incr());
        }
    }

    virtual void _mdecr (const Matrix& m)
    {
        if (this -> _continuous() && m._continuous())
        {
            _decr<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, this -> _p(m), m.incr());
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _decr<TR, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, this -> _p(m) + m._ldm() * i, m.incr());
        }
    }

    // matrix cleaning (we ALWAYS have mnIncr = 1 for matrices)
    virtual void _vanish()
    {
        CVM_ASSERT(this -> mpD, this -> mnSize * sizeof(TC))
        if (this -> _continuous())
        {
            memset (this -> mpD, 0, this -> mnSize * sizeof(TC));
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            memset (this -> mpD + this -> mnLD * i, 0, this -> mnM * sizeof(TC));
        }
    }

public:
    friend std::ostream& operator << <>(std::ostream& os, const Matrix<TR,TC>& mOut);
    friend std::istream& operator >> <>(std::istream& os, Matrix<TR,TC>& mIn);
};


// generalized square matrix
template <typename TR, typename TC>
class SqMatrix
{
protected:
    SqMatrix()
    {
    }

    virtual ~SqMatrix()
    {
    }

    virtual int _size    () const = 0;
    virtual int _msize   () const = 0;
    virtual int _nsize   () const = 0;
    virtual int _ld      () const = 0;
    virtual const TC* _p () const = 0;
    virtual TC* _p () = 0;

    // it differs from Matrix::_transp because in this case we can do it in-place.
    void _sq_transp()
    {
        const int mnM  = this -> _msize();
        const int mnLD = this -> _ld();
        TC* mpD = this -> _p();
        if (mnM > 1)
        {
            const int nM1 = mnLD + 1, nM1m = mnLD - 1, nM2m = mnM - 1;
            int i = 1, j = 1, m;
            for (;;)
            {
                m = mnM - i;
                __swap<TC> (m, mpD + j, 1, mpD + j + nM1m, mnLD);
                if (i >= nM2m)
                {
                    break;
                }
                ++i;
                j += nM1;
            }
        }
    }

    void _sq_plus_plus()                                        // plus identity
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnSize = this -> _size();
        const int nNext = this -> _msize() + 1;
        CVM_ASSERT(mpD, mnSize * sizeof(TC))
        for (int i = 0; i < mnSize; i += nNext)
        {
            mpD[i] += one;
        }
    }

    void _sq_minus_minus()                                      // minus identity
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnSize = this -> _size();
        const int nNext = this -> _msize() + 1;
        CVM_ASSERT(mpD, mnSize * sizeof(TC))
        for (int i = 0; i < mnSize; i += nNext)
        {
            mpD[i] -= one;
        }
    }

public:
    void _clean_low_triangle ()
    {
        const int mnM  = this -> _msize();
        const int mnLD = this -> _ld();
        TC* mpD = this -> _p();
        int n = 1;
        static const TR zero(0.);
        for (int i = 1; i < mnM; ++i)
        {
            __scal<TR,TC> (mpD + n, mnM - i, 1, zero);   // column by column
            n += mnLD + 1;
        }
    }
};

// matrix of real numbers
template <typename TR>
class basic_rmatrix : public Matrix <TR,TR>
{
    typedef std::complex<TR> TC;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef basic_rvector<TR> RVector;

    friend class basic_rvector<TR>; // _multiply

public:
    basic_rmatrix()
    {
    }

    basic_rmatrix (int nM, int nN)
        : BaseMatrix (nM, nN, nM, true)
    {
    }

    basic_rmatrix (TR* pD, int nM, int nN)
        : BaseMatrix (pD, nM, nN, nM, nM * nN)
    {
    }

    basic_rmatrix (const basic_rmatrix& m)
        : BaseMatrix (m.mnM, m.mnN, m.mnM, false)
    {
        this -> _massign(m);
    }

    explicit basic_rmatrix (const RVector& v, bool bBeColumn = true)            // true = column, false = row
        : BaseMatrix (v, bBeColumn)
    {
    }

    // submatrix constructor, 1-based
    basic_rmatrix (basic_rmatrix& m, int nRow, int nCol, int nHeight, int nWidth)
        : BaseMatrix (m._sub_pointer (nRow, nCol, nHeight, nWidth), nHeight, nWidth, m.ld(), nHeight * nWidth)
    {
        m._check_submatrix();
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    RVector operator () (int nI) throw (cvmexception)                          // returns COLUMN which IS an l-value
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    RVector operator [] (int nI) throw (cvmexception)                          // returns ROW which IS an l-value
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    const RVector operator () (int nI) const throw (cvmexception)              // returns column which IS NOT an l-value
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    const RVector operator [] (int nI) const throw (cvmexception)              // returns ROW which IS NOT an l-value
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    RVector diag (int nDiag) throw (cvmexception)                               // returns diagonal which IS l-value (shares memory) 
    {                                                                           // 0 - main, negative - low, positive - up
        return this -> _diag(nDiag);
    }

    const RVector diag (int nDiag) const throw (cvmexception)                   // returns diagonal which IS NOT l-value
    {                                                                           // 0 - main, negative - low, positive - up
        return this -> _diag(nDiag);
    }

    basic_rmatrix& operator = (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    basic_rmatrix& assign (const TR* pD)                                        // assigns foregn array (nIncr = 1)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_rmatrix& assign (int nRow, int nCol, const basic_rmatrix& m) throw (cvmexception)    // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.mnM + nRow - 1 > this -> mnM || m.mnN + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    basic_rmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_rmatrix& resize (int nNewM, int nNewN) throw (cvmexception)
    {
        this -> _resize (nNewM, nNewN);
        return *this;
    }

    bool operator == (const basic_rmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_rmatrix& m) const
    {
        return !operator == (m);
    }

    basic_rmatrix& operator << (const basic_rmatrix& m) throw (cvmexception)    // matrix replacement
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_rmatrix operator + (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_rmatrix operator - (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_rmatrix& sum (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_rmatrix& diff (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_rmatrix& operator += (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_rmatrix& operator -= (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_rmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_rmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    basic_rmatrix operator * (TR dMult) const
    {
        basic_rmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_rmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_rmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_rmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_rmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_rmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix
    basic_rmatrix operator ~ () const throw (cvmexception)
    {
        basic_rmatrix mRes (this -> mnN, this -> mnM);
        mRes._transp (*this);
        return mRes;
    }

    basic_rmatrix& transpose (const basic_rmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnN || this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        if (this -> mpD == m.mpD)
        {
            basic_rmatrix mTmp(m);
            this -> _transp (mTmp);
        }
        else
        {
            this -> _transp (m);
        }
        return *this;
    }

    // in-place transpose Matrix - changes dimensions
    basic_rmatrix& transpose() throw (cvmexception)
    {
        basic_rmatrix mTmp (*this);
        this -> _resize (this -> mnN, this -> mnM);
        this -> _transp (mTmp);
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    basic_rmatrix operator * (const basic_rmatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_rmatrix mRes (this -> mnM, m.mnN);
        mRes.mult (*this, m);
        return mRes;
    }

    // this = m1 * m2
    // overridden in srsmatrix, NOT in srmatrix
    basic_rmatrix& mult (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    // this = v_col * v_row (rank-1 update)
    basic_rmatrix& rank1update (const RVector& vCol, const RVector& vRow) throw (cvmexception)
    {
        static const TR one(1.);
        this -> _check_rank1update();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _vanish();
        __ger<TR, basic_rmatrix, RVector> (*this, vCol, vRow, one);
        return *this;
    }

    basic_rmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_rmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers
    basic_rmatrix& solve (const basic_srmatrix<TR>& mA, const basic_rmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize() || mA.msize() != this -> msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_rmatrix& solve (const basic_srmatrix<TR>& mA, const basic_rmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, mB, dErr);
    }

    basic_rmatrix& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize()  || mA.msize() != this -> msize() ||
            mA.msize() != mLU.msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_rmatrix& solve_lu (const basic_srmatrix<TR>& mA, const basic_srmatrix<TR>& mLU, const int* pPivots,
                             const basic_rmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, mB, dErr);
    }

    // singular value decomposition (values in decreasing order)
    RVector svd() const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vRes, NULL, NULL);
        return vRes;
    }

    RVector svd (basic_srmatrix<TR>& mU, basic_srmatrix<TR>& mVH) const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vRes, &mU, &mVH);
        return vRes;
    }

    // pseudo (generalized) inversion - const version
    basic_rmatrix pinv (TR threshold = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_rmatrix mAx(this -> mnN, this -> mnM);
        this -> _pinv (mAx, threshold);
        return mAx;
    }

    // pseudo (generalized) inversion - non-const version
    basic_rmatrix& pinv (const basic_rmatrix& mA, TR threshold = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        if (mA.msize() != this -> nsize() || mA.nsize() != this -> msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._pinv (*this, threshold);
        return *this;
    }

    int rank (TR eps = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        int nRank = 0;
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        vS.normalize();

        for (; nRank < vS.size(); ++nRank)
        {
            if (vS [nRank * this -> mnIncr + 1] < eps) break;
        }
        return nRank;
    }

    basic_rmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    // this += alpha * v_col * v_row (rank-1 update)
    basic_rmatrix& ger (TR alpha, const RVector& vCol, const RVector& vRow) throw (cvmexception)
    {
        this -> _check_ger();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __ger<TR, basic_rmatrix, RVector> (*this, vCol, vRow, alpha);
        return *this;
    }

    // this = alpha * m1 * m2 + beta * this
    basic_rmatrix& gemm (const basic_rmatrix& m1, bool bTrans1, const basic_rmatrix& m2, bool bTrans2, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        this -> _check_gemm();
        if (this -> mnM != (bTrans1 ? m1.mnN : m1.mnM) || 
            this -> mnN != (bTrans2 ? m2.mnM : m2.mnN) || 
            (bTrans1 ? m1.mnM : m1.mnN) != (bTrans2 ? m2.mnN : m2.mnM)) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _gemm (bTrans1, m1, bTrans2, m2, dAlpha, dBeta);
        return *this;
    }

    // this = alpha*a*b + beta*this of this = alpha*b*a + beta*this  where a is symmetric
    basic_rmatrix& symm (bool bLeft, const basic_srsmatrix<TR>& ms, const basic_rmatrix& m, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        this -> _check_symm();
        if (this -> mnM != m.mnM || this -> mnN != m.mnN ||
            ms.msize() != (bLeft ? this -> mnM : this -> mnN)) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _symm (bLeft, ms, m, dAlpha, dBeta);
        return *this;
    }

    basic_rmatrix& randomize (TR dFrom, TR dTo)
    {
        this -> _randomize (dFrom, dTo);
        return *this;
    }

    // 2-norm (maximum singular value)
    virtual TR norm2() const
    {
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        return vS[1];
    }

    virtual void _svd (RVector& vRes, basic_srmatrix<TR>* pmU, basic_srmatrix<TR>* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_rmatrix, basic_srmatrix<TR> > (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (basic_rmatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_rmatrix, basic_rmatrix> (mX, *this, threshold);
    }

    // ?gemm routines perform a matrix-matrix operation with general matrices. The operation is defined as
    // c := alpha*op(a)*op(b) + beta*c,
    // where: op(x) is one of op(x) = x or op(x) = x' or op(x) = conjg(x'),
    void _gemm (bool bTrans1, const basic_rmatrix& m1, bool bTrans2, const basic_rmatrix& m2, TR dAlpha, TR dBeta) throw (cvmexception)     // this = m1 * m2
    {
        basic_rmatrix mTmp1, mTmp2;
        const TR* pD1 = m1.get();
        const TR* pD2 = m2.get();
        if (this -> mpD == pD1) mTmp1 << m1;
        if (this -> mpD == pD2) mTmp2 << m2;
        __gemm<TR, basic_rmatrix> (this -> mpD == pD1 ? mTmp1 : m1, bTrans1, this -> mpD == pD2 ? mTmp2 : m2, bTrans2, dAlpha, *this, dBeta);
    }

    // this = alpha*a*b + beta*this or this = alpha*b*a + beta*this  where a is symmetric
    void _symm (bool bLeft, const basic_srsmatrix<TR>& ms, const basic_rmatrix& m, TR dAlpha, TR dBeta) throw (cvmexception)
    {
        basic_rmatrix mTmp;
        basic_srsmatrix<TR> msTmp;
        const TR* pD1 = ms.get();
        const TR* pD2 = m._pd();
        if (this -> mpD == pD1) msTmp << ms;
        if (this -> mpD == pD2) mTmp << m;
        __symm<TR, basic_srsmatrix<TR>, basic_rmatrix> (bLeft, this -> mpD == pD1 ? msTmp : ms, this -> mpD == pD2 ? mTmp : m, dAlpha, *this, dBeta);
    }

    virtual void _check_submatrix () {}

protected:
    // protected constructors for inherited stuff
    basic_rmatrix (int nM, int nN, int nLD, bool bZeroMemory)
        : BaseMatrix (nM, nN, nLD, bZeroMemory)
    {
    }

    basic_rmatrix (TR* pD, int nM, int nN, int nLD, int nSize)
        : BaseMatrix (pD, nM, nN, nLD, nSize)
    {
    }

    // returns diagonal which IS l-value (shares memory)
    // 0 - main, negative - low, positive - up
    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        int nShift = 0;
        int nSize = 0;
        if (nDiag >= 0)
        {
            if (nDiag >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nShift = nDiag * this -> mnLD;
            nSize = this -> mnN > this -> mnM ? (nDiag > this -> mnN - this -> mnM ? this -> mnN - nDiag : this -> mnM) : this -> mnN - nDiag;
        }
        else
        {
            nShift = - nDiag;
            if (nShift >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nSize = this -> mnM > this -> mnN ? (nShift > this -> mnM - this -> mnN ? this -> mnM - nShift : this -> mnN) : this -> mnM - nShift;
        }
        return RVector (this -> mpD + nShift, nSize, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_rmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // compares matrix elements (equal sizes assumed)
    bool _mequals (const basic_rmatrix& m) const
    {
        return ((*this) - m).norminf() <= basic_cvmMachMin<TR>();
    }

    // ?gemv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gemv calls all virtual methods inside
    void _gemv (bool bLeft, TR dAlpha, const RVector& v, TR dBeta, RVector& vRes) const
    {
        RVector vTmp;
        basic_rmatrix mTmp;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gemv<TR, basic_rmatrix, RVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha, 
                                                   vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    // 0-based, returns l-value sharing memory
    virtual RVector _row (int m)
    {
        return RVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns l-value sharing memory
    virtual RVector _col (int n)
    {
        return RVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    // 0-based
    virtual const RVector _row (int m) const
    {
        return RVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns l-value sharing memory
    virtual const RVector _col (int n) const
    {
        return RVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    virtual void _mult (const basic_rmatrix& m1, const basic_rmatrix& m2) throw (cvmexception)
    {
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        if (this -> mnM != m1.mnM || this -> mnN != m2.mnN || m1.mnN != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _gemm (false, m1, false, m2, one, zero);
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool bLeft) const
    {
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        this -> _gemv (bLeft, one, v, zero, vRes);
    }

    virtual void _randomize (TR dFrom, TR dTo)
    {
        if (this -> _continuous())
        {
            __randomize<TR> (this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __randomize<TR> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, dFrom, dTo);
        }
    }

    virtual void _check_ger() {}
    virtual void _check_rank1update() {}
    virtual void _check_gemm() {}
    virtual void _check_symm() {}
};


template <typename TR>
class basic_srmatrix : public basic_rmatrix<TR>, public SqMatrix <TR, TR>
{
    typedef std::complex<TR> TC;
    typedef basic_rvector<TR> RVector;
    typedef basic_cvector<TR,TC> CVector;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef SqMatrix<TR, TR> BaseSqMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;

public:
    basic_srmatrix()
    {
    }

    explicit basic_srmatrix (int nMN)
        : BaseRMatrix (nMN, nMN)
    {
    }

    basic_srmatrix (TR* pD, int nMN)
        : BaseRMatrix (pD, nMN, nMN)
    {
    }

    basic_srmatrix (const basic_srmatrix& m)
        : BaseRMatrix (m.mnM, m.mnN, m.mnM, false), BaseSqMatrix()
    {
        this -> _massign(m);
    }

    basic_srmatrix (const BaseRMatrix& m)
        : BaseRMatrix (m.msize(), m.nsize(), m.msize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
    }

    // diagonal square matrix constructor
    explicit basic_srmatrix (const RVector& v)
        : BaseRMatrix (v.size(), v.size(), v.size(), true)
    {
        __copy<TR> (this -> mnM, v, v.incr(), this -> mpD, this -> mnM + 1);
    }

    // submatrix constructor
    // 1-based
    basic_srmatrix (BaseRMatrix& m, int nRow, int nCol, int nSize)
        : BaseRMatrix (m, nRow, nCol, nSize, nSize)
    {
        m._check_submatrix();
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns COLUMN which CAN be l-value
    RVector operator () (int nFI) throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns ROW which CAN be l-value
    RVector operator [] (int nFI) throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns ROW which CAN NOT be l-value
    const RVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    basic_srmatrix& operator = (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    // assigns foregn array (nIncr = 1)
    basic_srmatrix& assign (const TR* pD)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_srmatrix& assign (int nRow, int nCol, const BaseRMatrix& m) throw (cvmexception)      // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.msize() + nRow - 1 > this -> mnM || m.nsize() + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.msize(), m.nsize(), m.ld());
        return *this;
    }

    // fills the content
    basic_srmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_srmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    basic_srmatrix& operator << (const basic_srmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srmatrix operator + (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srmatrix operator - (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_srmatrix& sum (const basic_srmatrix& m1, const basic_srmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_srmatrix& diff (const basic_srmatrix& m1, const basic_srmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srmatrix& operator += (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_srmatrix& operator -= (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_srmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srmatrix operator * (TR dMult) const
    {
        basic_srmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srmatrix& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    basic_srmatrix operator ~ () const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        return mRes.transpose();
    }

    basic_srmatrix& transpose (const basic_srmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return this -> transpose();
    }

    basic_srmatrix& transpose()
    {
        this -> _transp();
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        return this -> BaseRMatrix::operator * (v);
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        return this -> BaseRMatrix::operator * (m);
    }

    basic_srmatrix operator * (const basic_srmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srmatrix mRes (this -> mnM);
        mRes.mult (*this, m);
        return mRes;
    }

    basic_srmatrix& operator *= (const basic_srmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        const basic_srmatrix mTmp (*this);
        this -> mult (mTmp, m);
        return *this;
    }

    basic_srmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_srmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers Ax=b. Also return solution flavor.
    RVector solve (const RVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, NULL, NULL);
        return vX;
    }

    RVector solve (const RVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (vB, dErr);
    }

    BaseRMatrix solve (const BaseRMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseRMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, NULL, NULL);
        return mX;
    }

    BaseRMatrix solve (const BaseRMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mB, dErr);
    }

    RVector solve_lu (const basic_srmatrix& mLU, const int* pPivots, const RVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, mLU, pPivots);
        return vX;
    }

    RVector solve_lu (const basic_srmatrix& mLU, const int* pPivots, const RVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, vB, dErr);
    }

    BaseRMatrix solve_lu (const basic_srmatrix& mLU, const int* pPivots, const BaseRMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseRMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, mLU, pPivots);
        return mX;
    }

    BaseRMatrix solve_lu (const basic_srmatrix& mLU, const int* pPivots, const BaseRMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, mB, dErr);
    }

    // matrix determinant
    TR det() const throw (cvmexception)
    {
        return this -> _det();
    }

    // low-up factorization
    basic_srmatrix& low_up (const basic_srmatrix& m, int* nPivots) throw (cvmexception)
    {
        (*this) = m;
        this -> _low_up (nPivots);
        return *this;
    }

    basic_srmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_srmatrix mRes (*this);
        mRes._low_up (nPivots);
        return mRes;
    }

    // reciprocal of the condition number
    TR cond() const throw (cvmexception)
    {
        TR dCondNum(0.);
        __cond_num<TR, basic_srmatrix>(*this, dCondNum);        // universal method, no need to virtualize
        return dCondNum;
    }

    // matrix inversion
    basic_srmatrix& inv (const basic_srmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_srmatrix>(*this, mArg);                     // overridden in srsmatrix, no need to virtualize
        return *this;
    }

    basic_srmatrix inv() const throw (cvmexception)
    {
        basic_srmatrix mRes (this -> mnM);
        __inv<basic_srmatrix>(mRes, *this);                     // overridden in srsmatrix, no need to virtualize
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_srmatrix& exp (const basic_srmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp<basic_srmatrix, TR>(*this, mArg, tol);            // uses universal code inside - no need to virtualize
        return *this;
    }

    basic_srmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_srmatrix mRes (this -> mnM);
        __exp<basic_srmatrix, TR>(mRes, *this, tol);
        return mRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_srmatrix& polynom (const basic_srmatrix& m, const RVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), v.incr() > 1 ? v1 : v);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_srmatrix polynom (const RVector& v) const
    {
        basic_srmatrix mRes (this -> mnM);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (mRes.mpD, mRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, v.incr() > 1 ? v1 : v);
        return mRes;
    }

    // eigenvalues
    CVector eig (basic_scmatrix<TR,TC>& mEigVect, bool bRightVect = true) const throw (cvmexception)
    {
        CVector vEig(this -> mnM);
        this -> _eig (vEig, &mEigVect, bRightVect);
        return vEig;
    }

    CVector eig() const throw (cvmexception)
    {
        CVector vEig (this -> mnM);
        this -> _eig (vEig, NULL, true);
        return vEig;
    }

    // Cholesky factorization
    basic_srmatrix& cholesky (const basic_srsmatrix<TR>& m) throw (cvmexception)
    {
        this -> _check_cholesky();  // doesn't work for band matrices
        *this = m;
        int nOutInfo = __cholesky<basic_srmatrix> (*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        this -> _clean_low_triangle ();
        return *this;
    }

    // Bunch-Kaufman factorization
    basic_srmatrix& bunch_kaufman (const basic_srsmatrix<TR>& m, int* nPivots) throw (cvmexception)
    {
        this -> _check_bunch_kaufman();  // doesn't work for band matrices
        *this = m;
        __bunch_kaufman<basic_srmatrix> (*this, nPivots);
        return *this;
    }

    basic_srmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_srmatrix& randomize (TR dFrom, TR dTo)
    {
        this -> _randomize (dFrom, dTo);
        return *this;
    }

    virtual void _eig (CVector& vEig, basic_scmatrix<TR,TC>* mEigVect, bool bRightVect) const throw (cvmexception)
    {
        __eig<CVector, basic_srmatrix, basic_scmatrix<TR,TC> > (vEig, *this, mEigVect, bRightVect);
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        RVector vB1;
        RVector vX1;
        if (vB.incr() > 1) vB1 << vB;   // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TR, basic_srmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TR, basic_srmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

protected:
    // protected constructors for inherited stuff
    basic_srmatrix (int nMN, int nLD, bool bZeroMemory)
        : BaseRMatrix (nMN, nMN, nLD, bZeroMemory)
    {
    }

    basic_srmatrix (TR* pD, int nMN, int nLD, int nSize)
        : BaseRMatrix (pD, nMN, nMN, nLD, nSize)
    {
    }

    virtual int _size   () const {return this -> mnSize;}
    virtual int _msize  () const {return this -> mnM;}
    virtual int _nsize  () const {return this -> mnN;}
    virtual int _ld     () const {return this -> mnLD;}
    virtual const TR* _p() const {return this -> mpD;}
    virtual TR*       _p()       {return this -> mpD;}

    // returns diagonal which IS l-value (shares memory)
    // 0 - main, negative - low, positive - up
    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nD >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        return RVector (this -> mpD + (nDiag > 0 ? nDiag * this -> mnLD : nD), this -> mnM - nD, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_srmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // returns main diagonal of low_up factorization
    virtual RVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _transp()
    {
        this -> _sq_transp();
    }

    virtual void _plus_plus()
    {
        this -> _sq_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _sq_minus_minus();
    }

    virtual TR _det() const throw (cvmexception)
    {
        TR dDet(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                dDet = this -> _ij_val (0, 0);
                break;
            case 2:
                dDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TR one(1.);
                    basic_array<int> naPivots (this -> mnM);
                    RVector vUpDiag = this -> _low_up_diag (naPivots);

                    dDet = one;
                    for (int i = 1; i <= this -> mnM; ++i)
                    {
                        dDet *= vUpDiag[i];
                        if (i != naPivots[i]) dDet = - dDet;
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_SINGULARMATRIX) throw e;
                }
                break;
        }
        return dDet;
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_srmatrix>(*this, nPivots);
    }

    virtual void _check_cholesky() {}
    virtual void _check_bunch_kaufman() {}
};


// matrix of complex numbers
template <typename TR, typename TC>
class basic_cmatrix : public Matrix<TR,TC>
{
    typedef basic_rvector<TR> RVector;
    typedef basic_cvector<TR,TC> CVector;
    typedef Array<TR,TC> BaseArray;
    typedef Matrix<TR,TC> BaseMatrix;

    friend class basic_cvector<TR,TC>;     // _multiply

public:
    basic_cmatrix()
    {
    }

    basic_cmatrix (int nM, int nN)
        : BaseMatrix (nM, nN, nM, true)
    {
    }

    basic_cmatrix (TC* pD, int nM, int nN)
        : BaseMatrix (pD, nM, nN, nM, nM * nN)
    {
    }

    basic_cmatrix (const basic_cmatrix& m)
        : BaseMatrix (m.mnM, m.mnN, m.mnM, false)
    {
        this -> _massign(m);
    }

    // true = column, false = row
    explicit basic_cmatrix (const CVector& v, bool bBeColumn = true)
        : BaseMatrix (v, bBeColumn)
    {
    }

    explicit basic_cmatrix (const basic_rmatrix<TR>& m, bool bRealPart = true)
        : BaseMatrix (m.msize(), m.nsize(), m.msize(), true)
    {
        if (bRealPart)
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, m._pd(), NULL);
        else
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, NULL, m._pd());
    }

    basic_cmatrix (const TR* pRe, const TR* pIm, int nM, int nN)
        : BaseMatrix (nM, nN, nM, true)
    {
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, pRe, pIm, 1, 1);
    }

    basic_cmatrix (const basic_rmatrix<TR>& mRe, const basic_rmatrix<TR>& mIm)
        : BaseMatrix (mRe.msize(), mRe.nsize(), mRe.msize(), false)
    {
        if (mRe.msize() != mIm.msize() || mRe.nsize() != mIm.nsize() || mRe.ld() != mIm.ld()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe, mIm, mRe.incr(), mIm.incr());
    }

    // submatrix constructor
    // 1-based
    basic_cmatrix (basic_cmatrix& m, int nRow, int nCol, int nHeight, int nWidth)
        : BaseMatrix (m._sub_pointer (nRow, nCol, nHeight, nWidth), nHeight, nWidth, m.ld(), nHeight * nWidth)
    {
        m._check_submatrix();
    }

    type_proxy<TC,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TC operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns COLUMN which IS an l-value
    CVector operator () (int nI) throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    // returns COLUMN which IS NOT an l-value
    const CVector operator () (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    // returns ROW which IS an l-value
    CVector operator [] (int nI) throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    // returns ROW which IS NOT an l-value
    const CVector operator [] (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    // real part
    const basic_rmatrix<TR> real() const
    {
        basic_rmatrix<TR> mRet (this -> mnM, this -> mnN);
        __copy<TR> (this -> mnSize, __get_real_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    // imaginary part
    const basic_rmatrix<TR> imag() const
    {
        basic_rmatrix<TR> mRet (this -> mnM, this -> mnN);
        __copy<TR> (this -> mnSize, __get_imag_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    // returns diagonal which IS l-value (shares memory) 
    // 0 - main, negative - low, positive - up
    CVector diag (int nDiag) throw (cvmexception)
    {
        return this -> _diag(nDiag);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    const CVector diag (int nDiag) const throw (cvmexception)
    {
        return this -> _diag(nDiag);
    }

    basic_cmatrix& operator = (const basic_cmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize() || this -> mnN != m.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    // assigns a foregn array (nIncr = 1)
    basic_cmatrix& assign (const TC* pD)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_cmatrix& assign (int nRow, int nCol, const basic_cmatrix& m) throw (cvmexception)       // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.mnM + nRow - 1 > this -> mnM || m.mnN + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    basic_cmatrix& set (TC c)
    {
        this -> _set (c);
        return *this;
    }

    basic_cmatrix& assign_real (const basic_rmatrix<TR>& mRe) throw (cvmexception)
    {
        if (this -> mnM != mRe.msize() || this -> mnN != mRe.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe._pd(), mRe.incr());
        return *this;
    }

    basic_cmatrix& assign_imag (const basic_rmatrix<TR>& mIm) throw (cvmexception)
    {
        if (this -> mnM != mIm.msize() || this -> mnN != mIm.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_imag<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mIm._pd(), mIm.incr());
        return *this;
    }

    // fills real part
    basic_cmatrix& set_real (TR d)
    {
        this -> _set_real_number (d);
        return *this;
    }

    // fills imaginary part
    basic_cmatrix& set_imag (TR d)
    {
        this -> _set_imag_number (d);
        return *this;
    }

    basic_cmatrix& resize (int nNewM, int nNewN) throw (cvmexception)
    {
        this -> _resize (nNewM, nNewN);
        return *this;
    }

    bool operator == (const basic_cmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_cmatrix& m) const
    {
        return !operator == (m);
    }

    basic_cmatrix& operator << (const basic_cmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_cmatrix operator + (const basic_cmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_cmatrix operator - (const basic_cmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_cmatrix& sum (const basic_cmatrix& m1, const basic_cmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_cmatrix& diff (const basic_cmatrix& m1, const basic_cmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m1.mnN || this -> mnM != m2.mnM || this -> mnN != m2.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_cmatrix& operator += (const basic_cmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_cmatrix& operator -= (const basic_cmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_cmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_cmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    basic_cmatrix operator * (TR dMult) const
    {
        basic_cmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_cmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_cmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_cmatrix operator * (TC cMult) const
    {
        basic_cmatrix mRes (*this);
        mRes._scal (cMult);
        return mRes;
    }

    basic_cmatrix operator / (TC cDiv) const throw (cvmexception)
    {
        basic_cmatrix mRes (*this);
        mRes._div (cDiv);
        return mRes;
    }

    basic_cmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_cmatrix& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_cmatrix& operator *= (TC cMult)
    {
        this -> _scal (cMult);
        return *this;
    }

    basic_cmatrix& operator /= (TC cDiv)
    {
        this -> _div (cDiv);
        return *this;
    }

    basic_cmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // hermitian conjugated matrix
    basic_cmatrix operator ~ () const throw (cvmexception)
    {
        basic_cmatrix mRes (this -> mnN, this -> mnM);
        mRes._transp (*this);
        __conj<TC> (mRes.mpD, mRes.mnSize, mRes.mnIncr);
        return mRes;
    }

    basic_cmatrix& conj (const basic_cmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnN || this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        if (this -> mpD == m.mpD)
        {
            basic_cmatrix mTmp(m);
            this -> _transp (mTmp);
        }
        else
        {
            this -> _transp (m);
        }
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    basic_cmatrix& conj () throw (cvmexception)
    {
        basic_cmatrix mTmp (*this);
        this -> _resize (this -> mnN, this -> mnM);
        this -> _transp (mTmp);
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    CVector operator * (const CVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    basic_cmatrix operator * (const basic_cmatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_cmatrix mRes (this -> mnM, m.mnN);
        mRes.mult (*this, m);
        return mRes;
    }

    // this = m1 * m2
    basic_cmatrix& mult (const basic_cmatrix& m1, const basic_cmatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    // this = v_col * v_row (rank-1 update, unconjugated)
    basic_cmatrix& rank1update_u (const CVector& vCol, const CVector& vRow) throw (cvmexception)
    {
        static const TC one(1.);
        this -> _check_rank1update_u();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _vanish();
        __geru<TC, basic_cmatrix, CVector> (*this, vCol, vRow, one);
        return *this;
    }

    // this = v_col * conj (v_row) (rank-1 update, conjugated)
    basic_cmatrix& rank1update_c (const CVector& vCol, const CVector& vRow) throw (cvmexception)
    {
        static const TC one(1.);
        this -> _check_rank1update_c();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _vanish();
        __gerc<TC, basic_cmatrix, CVector> (*this, vCol, vRow, one);
        return *this;
    }

    basic_cmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_cmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers
    basic_cmatrix& solve (const basic_scmatrix<TR,TC>& mA, const basic_cmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize() || mA.msize() != this -> msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, NULL, NULL);
        return *this;
    }

    basic_cmatrix& solve (const basic_scmatrix<TR,TC>& mA, const basic_cmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mA, mB, dErr);
    }

    basic_cmatrix& solve_lu (const basic_scmatrix<TR,TC>& mA, const basic_scmatrix<TR,TC>& mLU, const int* pPivots,
                             const basic_cmatrix& mB, TR& dErr) throw (cvmexception)
    {
        if (mA.msize() != mB.msize() || mA.msize() != this -> msize() || mA.msize() != mLU.msize() || mB.nsize() != this -> nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._solve (mB, *this, dErr, mLU, pPivots);
        return *this;
    }

    basic_cmatrix& solve_lu (const basic_scmatrix<TR,TC>& mA, const basic_scmatrix<TR,TC>& mLU, const int* pPivots,
                             const basic_cmatrix& mB) throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mA, mLU, pPivots, mB, dErr);
    }

    RVector svd() const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this ->_svd (vRes, NULL, NULL);
        return vRes;
    }

    RVector svd (basic_scmatrix<TR,TC>& mU, basic_scmatrix<TR,TC>& mVH) const throw (cvmexception)
    {
        RVector vRes (_cvm_min<int>(this -> mnM, this -> mnN));
        this ->_svd (vRes, &mU, &mVH);
        return vRes;
    }

    // pseudo (generalized) inversion - const version
    basic_cmatrix pinv (TR threshold = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_cmatrix mAx(this -> mnN, this -> mnM);
        this -> _pinv (mAx, threshold);
        return mAx;
    }

    // pseudo (generalized) inversion - non-const version
    basic_cmatrix& pinv (const basic_cmatrix& mA, TR threshold = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        if (mA.msize() != this -> nsize() || mA.nsize() != this -> msize()) throw cvmexception (CVM_SIZESMISMATCH);
        mA._pinv (*this, threshold);
        return *this;
    }

    int rank (TR eps = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        int nRank = 0;
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        vS.normalize();

        for (; nRank < vS.size(); ++nRank)
        {
            if (vS [nRank * this -> mnIncr + 1] < eps) break;
        }
        return nRank;
    }

    basic_cmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    // this += alpha * v_col * v_row (unconjugated)
    basic_cmatrix& geru (TC alpha, const CVector& vCol, const CVector& vRow) throw (cvmexception)
    {
        this -> _check_geru();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __geru<TC, basic_cmatrix, CVector> (*this, vCol, vRow, alpha);
        return *this;
    }

    // this += alpha * v_col * v_row (unconjugated)
    basic_cmatrix& gerc (TC alpha, const CVector& vCol, const CVector& vRow) throw (cvmexception)
    {
        this -> _check_gerc();
        if (this -> mnM != vCol.size() || this -> mnN != vRow.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __gerc<TC, basic_cmatrix, CVector> (*this, vCol, vRow, alpha);
        return *this;
    }

    // this = alpha * m1 * m2 + beta * this
    basic_cmatrix& gemm (const basic_cmatrix& m1, bool bTrans1, const basic_cmatrix& m2, bool bTrans2, TC dAlpha, TC dBeta) throw (cvmexception)
    {
        this -> _check_gemm();
        if (this -> mnM != (bTrans1 ? m1.mnN : m1.mnM) || 
            this -> mnN != (bTrans2 ? m2.mnM : m2.mnN) || 
            (bTrans1 ? m1.mnM : m1.mnN) != (bTrans2 ? m2.mnN : m2.mnM)) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _gemm (bTrans1, m1, bTrans2, m2, dAlpha, dBeta);
        return *this;
    }

    // this = alpha*a*b + beta*this of this = alpha*b*a + beta*this  where a is hermitian
    basic_cmatrix& hemm (bool bLeft, const basic_schmatrix<TR,TC>& ms, const basic_cmatrix& m, TC dAlpha, TC dBeta) throw (cvmexception)
    {
        this -> _check_hemm();
        if (ms.msize() != (bLeft ? this -> mnM : this -> mnN) || 
            this -> mnM != m.mnM || this -> mnN != m.mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _hemm (bLeft, ms, m, dAlpha, dBeta);
        return *this;
    }

    basic_cmatrix& randomize_real (TR dFrom, TR dTo)
    {
        this -> _randomize_real (dFrom, dTo);
        return *this;
    }

    basic_cmatrix& randomize_imag (TR dFrom, TR dTo)
    {
        this -> _randomize_imag (dFrom, dTo);
        return *this;
    }

    // 2-norm (maximum singular value)
    virtual TR norm2() const
    {
        RVector vS (_cvm_min<int>(this -> mnM, this -> mnN));
        this -> _svd (vS, NULL, NULL);
        return vS[1];
    }

    // singular values in decreasing order
    virtual void _svd (RVector& vRes, basic_scmatrix<TR,TC>* pmU, basic_scmatrix<TR,TC>* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_cmatrix, basic_scmatrix<TR,TC> > (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (basic_cmatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_cmatrix, basic_cmatrix> (mX, *this, threshold);
    }

    virtual void _check_submatrix () {}

    virtual void _scal (TC d)
    {
        if (this -> _continuous())
        {
            __scal<TC, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __scal<TC, TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, d);
        }
    }

    void _div (TC d) throw (cvmexception)
    {
        if (_abs(d) <= basic_cvmMachMin<TR>()) throw cvmexception (CVM_DIVISIONBYZERO);
        static const TC one(1.);
        this -> _scal (one / d);
    }

    // ?gemm routines perform a matrix-matrix operation with general matrices. The operation is defined as
    // c := alpha*op(a)*op(b) + beta*c,
    // where: op(x) is one of op(x) = x or op(x) = x' or op(x) = conjg(x'),
    void _gemm (bool bTrans1, const basic_cmatrix& m1, bool bTrans2, const basic_cmatrix& m2, TC dAlpha, TC dBeta) throw (cvmexception)
    {
        basic_cmatrix mTmp1, mTmp2;
        const TC* pD1 = m1.get();
        const TC* pD2 = m2.get();
        if (this -> mpD == pD1) mTmp1 << m1;
        if (this -> mpD == pD2) mTmp2 << m2;
        __gemm<TC, basic_cmatrix> (this -> mpD == pD1 ? mTmp1 : m1, bTrans1, this -> mpD == pD2 ? mTmp2 : m2, bTrans2, dAlpha, *this, dBeta);
    }

    // this = alpha*a*b + beta*this or this = alpha*b*a + beta*this  where a is hermitian
    void _hemm (bool bLeft, const basic_schmatrix<TR,TC>& ms, const basic_cmatrix& m, TC dAlpha, TC dBeta) throw (cvmexception)
    {
        basic_cmatrix mTmp;
        basic_schmatrix<TR,TC> msTmp;
        const TC* pD1 = ms.get();
        const TC* pD2 = m._pd();
        if (this -> mpD == pD1) msTmp << ms;
        if (this -> mpD == pD2) mTmp << m;
        __hemm<TC, basic_schmatrix<TR,TC>, basic_cmatrix> (bLeft, this -> mpD == pD1 ? msTmp : ms, this -> mpD == pD2 ? mTmp : m, dAlpha, *this, dBeta);
    }

protected:
    // protected constructors for inherited stuff
    basic_cmatrix (int nM, int nN, int nLD, bool bZeroMemory)
        : BaseMatrix (nM, nN, nLD, bZeroMemory)
    {
    }

    basic_cmatrix (TC* pD, int nM, int nN, int nLD, int nSize)
        : BaseMatrix (pD, nM, nN, nLD, nSize)
    {
    }

    // returns diagonal which IS l-value (shares memory) 
    // 0 - main, negative - low, positive - up
    virtual CVector _diag (int nDiag) throw (cvmexception)
    {
        int nShift = 0;
        int nSize = 0;
        if (nDiag >= 0)
        {
            if (nDiag >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nShift = nDiag * this -> mnLD;
            nSize = this -> mnN > this -> mnM ? (nDiag > this -> mnN - this -> mnM ? this -> mnN - nDiag : this -> mnM) : this -> mnN - nDiag;
        }
        else
        {
            nShift = - nDiag;
            if (nShift >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
            nSize = this -> mnM > this -> mnN ? (nShift > this -> mnM - this -> mnN ? this -> mnM - nShift : this -> mnN) : this -> mnM - nShift;
        }
        return CVector (this -> mpD + nShift, nSize, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (shares memory) 
    // 0 - main, negative - low, positive - up
    virtual const CVector _diag (int nDiag) const throw (cvmexception)
    {
        CVector vRet (const_cast<basic_cmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // compares matrix elements (equal sizes assumed)
    bool _mequals (const basic_cmatrix& m) const
    {
        return ((*this) - m).norminf() <= basic_cvmMachMin<TR>();
    }

    // ?gemv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gemv calls all virtual methods inside
    void _gemv (bool bLeft, TC dAlpha, const CVector& v, TC dBeta, CVector& vRes) const
    {
        CVector vTmp;
        basic_cmatrix mTmp;
        const TC* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gemv<TC, basic_cmatrix, CVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha, 
                                                   vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    // 0-based, returns an l-value sharing memory
    virtual CVector _row (int m)
    {
        return CVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns an l-value sharing memory
    virtual CVector _col (int n)
    {
        return CVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    // 0-based, returns not an l-value (copy)
    virtual const CVector _row (int m) const
    {
        return CVector (this -> mpD + m, this -> mnN, this -> mnLD);
    }

    // 0-based, returns not an l-value (copy)
    virtual const CVector _col (int n) const
    {
        return CVector (this -> mpD + this -> mnLD * n, this -> mnM);
    }

    virtual void _multiply (CVector& vRes, const CVector& v, bool bLeft) const
    {
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        this -> _gemv (bLeft, one, v, zero, vRes);
    }

    virtual void _mult (const basic_cmatrix& m1, const basic_cmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnN != m2.mnN || m1.mnN != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        this -> _gemm (false, m1, false, m2, one, zero);
    }

    virtual void _set_real_number (TR d)
    {
        if (this -> _continuous())
        {
            _set_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _set_real<TR,TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, d);
        }
    }

    virtual void _set_imag_number (TR d)
    {
        if (this -> _continuous())
        {
            _set_imag<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            _set_imag<TR,TC> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, d);
        }
    }

    virtual void _randomize_real (TR dFrom, TR dTo)
    {
        if (this -> _continuous())
        {
            __randomize_real<TC,TR> (this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __randomize_real<TC,TR> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, dFrom, dTo);
        }
    }

    virtual void _randomize_imag (TR dFrom, TR dTo)
    {
        if (this -> _continuous())
        {
            __randomize_imag<TC,TR> (this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
        }
        else for (int i = 0; i < this -> mnN; ++i)
        {
            __randomize_imag<TC,TR> (this -> mpD + this -> mnLD * i, this -> mnM, this -> mnIncr, dFrom, dTo);
        }
    }

    virtual void _check_geru() {}
    virtual void _check_gerc() {}
    virtual void _check_rank1update_u() {}
    virtual void _check_rank1update_c() {}
    virtual void _check_gemm() {}
    virtual void _check_hemm() {}
};


// square matrix of complex numbers
template <typename TR, typename TC>
class basic_scmatrix : public basic_cmatrix<TR,TC>, public SqMatrix<TR,TC>
{
    typedef basic_cvector<TR,TC> CVector;
    typedef Array<TR,TC> BaseArray;
    typedef Matrix<TR,TC> BaseMatrix;
    typedef SqMatrix<TR, TC> BaseSqMatrix;
    typedef basic_cmatrix<TR,TC> BaseCMatrix;

public:
    basic_scmatrix()
    {
    }

    explicit basic_scmatrix (int nMN)
        : BaseCMatrix (nMN, nMN)
    {
    }

    basic_scmatrix (TC* pD, int nMN)
        : BaseCMatrix (pD, nMN, nMN)
    {
    }

    basic_scmatrix (const basic_scmatrix& m)
        : BaseCMatrix (m.mnM, m.mnN, m.mnM, false), BaseSqMatrix()
    {
        this -> _massign(m);
    }

    basic_scmatrix (const BaseCMatrix& m)
        : BaseCMatrix (m.msize(), m.nsize(), m.msize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
    }

    // diagonal square matrix constructor
    explicit basic_scmatrix (const CVector& v)
        : BaseCMatrix (v.size(), v.size(), v.size(), true)
    {
        __copy<TC> (this -> mnM, v, v.incr(), this -> mpD, this -> mnM + 1);
    }

    explicit basic_scmatrix (const basic_srmatrix<TR>& m, bool bRealPart = true)
        : BaseCMatrix (m.msize(), m.msize(), m.msize(), true)
    {
        if (bRealPart)
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, m._pd(), NULL);
        else
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, NULL, m._pd());
    }

    basic_scmatrix (const TR* pRe, const TR* pIm, int nMN)
        : BaseCMatrix (pRe, pIm, nMN, nMN)
    {
    }

    basic_scmatrix (const basic_srmatrix<TR>& mRe, const basic_srmatrix<TR>& mIm)
        : BaseCMatrix (mRe.msize(), mRe.nsize(), mRe.msize(), false)
    {
        if (mRe.msize() != mIm.msize() || mRe.nsize() != mIm.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe, mIm, mRe.incr(), mIm.incr());
    }

    // submatrix constructor
    // 1-based
    basic_scmatrix (BaseCMatrix& m, int nRow, int nCol, int nSize)
        : BaseCMatrix (m, nRow, nCol, nSize, nSize)
    {
        m._check_submatrix();
    }

    type_proxy<TC,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TC operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns COLUMN which CAN be l-value
    CVector operator () (int nFI) throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns column which CAN NOT be l-value
    const CVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns ROW which CAN be l-value
    CVector operator [] (int nFI) throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns ROW which CAN NOT be l-value
    const CVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // real part
    const basic_srmatrix<TR> real() const
    {
        basic_srmatrix<TR> mRet (this -> mnM);
        __copy<TR> (this -> mnSize, __get_real_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    // imaginary part
    const basic_srmatrix<TR> imag() const
    {
        basic_srmatrix<TR> mRet (this -> mnM);
        __copy<TR> (this -> mnSize, __get_imag_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    basic_scmatrix& operator = (const basic_scmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    // assigns foregn array (nIncr = 1)
    basic_scmatrix& assign (const TC* pD)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_scmatrix& assign (int nRow, int nCol, const BaseCMatrix& m) throw (cvmexception)         // submatrix assignment
    {
        if (nRow <= 0 || nRow > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRow);
        if (nCol <= 0 || nCol > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nCol);
        if (m.msize() + nRow - 1 > this -> mnM || m.nsize() + nCol - 1 > this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRow, nCol), m._pd(), m.msize(), m.nsize(), m.ld());
        return *this;
    }

    // fills the content
    basic_scmatrix& set (TC c)
    {
        this -> _set (c);
        return *this;
    }

    // assigns real matrix
    basic_scmatrix& assign_real (const basic_srmatrix<TR>& mRe) throw (cvmexception)
    {
        if (this -> mnM != mRe.msize() || this -> mnN != mRe.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe._pd(), mRe.incr());
        return *this;
    }

    // assigns imaginary matrix
    basic_scmatrix& assign_imag (const basic_srmatrix<TR>& mIm) throw (cvmexception)
    {
        if (this -> mnM != mIm.msize() || this -> mnN != mIm.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_imag<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mIm._pd(), mIm.incr());
        return *this;
    }

    // fills real part
    basic_scmatrix& set_real (TR d)
    {
        this -> _set_real_number (d);
        return *this;
    }

    // fills imaginary part
    basic_scmatrix& set_imag (TR d)
    {
        this -> _set_imag_number (d);
        return *this;
    }

    basic_scmatrix& resize (int nNewM) throw (cvmexception)
    {
        this -> _resize (nNewM, nNewM);
        return *this;
    }

    basic_scmatrix& operator << (const basic_scmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_scmatrix operator + (const basic_scmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_scmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_scmatrix operator - (const basic_scmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_scmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    // this = v1 + v2
    basic_scmatrix& sum (const basic_scmatrix& m1, const basic_scmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    // this = v1 + v2
    basic_scmatrix& diff (const basic_scmatrix& m1, const basic_scmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_scmatrix& operator += (const basic_scmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_scmatrix& operator -= (const basic_scmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_scmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_scmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_scmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_scmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_scmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_scmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_scmatrix operator * (TR dMult) const
    {
        basic_scmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_scmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_scmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_scmatrix operator * (TC cMult) const
    {
        basic_scmatrix mRes (*this);
        mRes._scal (cMult);
        return mRes;
    }

    basic_scmatrix operator / (TC cDiv) const throw (cvmexception)
    {
        basic_scmatrix mRes (*this);
        mRes._div (cDiv);
        return mRes;
    }

    basic_scmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_scmatrix& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_scmatrix& operator *= (TC cMult)
    {
        this -> _scal (cMult);
        return *this;
    }

    basic_scmatrix& operator /= (TC cDiv)
    {
        this -> _div (cDiv);
        return *this;
    }

    basic_scmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // hermitian conjugated Matrix
    basic_scmatrix operator ~ () const throw (cvmexception)
    {
        basic_scmatrix mRes (*this);
        mRes._sq_transp();
        __conj<TC> (mRes.mpD, mRes.mnSize, mRes.mnIncr);
        return mRes;
    }

    // this = hermitian conjugate (m)
    basic_scmatrix& conj (const basic_scmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return this -> conj();
    }

    // this = hermitian conjugate (this)
    basic_scmatrix& conj() throw (cvmexception)
    {
        this -> _transp();
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    CVector operator * (const CVector& v) const throw (cvmexception)
    {
        return this -> BaseCMatrix::operator * (v);
    }

    // special exclusion since matrix product is not commutative
    BaseCMatrix operator * (const BaseCMatrix& m) const throw (cvmexception)
    {
        return this -> BaseCMatrix::operator * (m);
    }

    basic_scmatrix operator * (const basic_scmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        basic_scmatrix mRes (this -> mnM);
        mRes.mult (*this, m);
        return mRes;
    }

    basic_scmatrix& operator *= (const basic_scmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        const basic_scmatrix mTmp (*this);
        this -> mult (mTmp, m);
        return *this;
    }

    basic_scmatrix& swap_rows (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_rows (n1, n2);
        return *this;
    }

    basic_scmatrix& swap_cols (int n1, int n2) throw (cvmexception)
    {
        this -> _swap_cols (n1, n2);
        return *this;
    }

    // linear solvers
    CVector solve (const CVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (vB, dErr);
    }

    BaseCMatrix solve (const BaseCMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve (mB, dErr);
    }

    CVector solve (const CVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, NULL, NULL);
        return vX;
    }

    BaseCMatrix solve (const BaseCMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseCMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, NULL, NULL);
        return mX;
    }

    CVector solve_lu (const basic_scmatrix& mLU, const int* pPivots, const CVector& vB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, vB, dErr);
    }

    BaseCMatrix solve_lu (const basic_scmatrix& mLU, const int* pPivots, const BaseCMatrix& mB) const throw (cvmexception)
    {
        static TR dErr(0.);
        return this -> solve_lu (mLU, pPivots, mB, dErr);
    }

    CVector solve_lu (const basic_scmatrix& mLU, const int* pPivots, const CVector& vB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != vB.size() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vX (this -> mnM);
        this -> _solve (vB, vX, dErr, mLU, pPivots);
        return vX;
    }

    BaseCMatrix solve_lu (const basic_scmatrix& mLU, const int* pPivots, const BaseCMatrix& mB, TR& dErr) const throw (cvmexception)
    {
        if (this -> mnM != mB.msize() || this -> mnM != mLU.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseCMatrix mX (mB.msize(), mB.nsize());
        this -> _solve (mB, mX, dErr, mLU, pPivots);
        return mX;
    }

    // matrix determinant
    TC det() const throw (cvmexception)
    {
        return this -> _det();
    }

    // low-up factorization
    basic_scmatrix& low_up (const basic_scmatrix& m, int* nPivots) throw (cvmexception)
    {
        (*this) = m;
        this -> _low_up(nPivots);
        return *this;
    }

    basic_scmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_scmatrix mRes (*this);
        mRes._low_up(nPivots);
        return mRes;
    }

    // reciprocal of the condition number
    TR cond() const throw (cvmexception)
    {
        TR dCondNum(0.);
        __cond_num<TR, basic_scmatrix>(*this, dCondNum);        // universal method, no need to virtualize
        return dCondNum;
    }

    // matrix inversion
    basic_scmatrix& inv (const basic_scmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_scmatrix>(*this, mArg);
        return *this;
    }

    basic_scmatrix inv() const throw (cvmexception)
    {
        basic_scmatrix mRes (this -> mnM);
        __inv<basic_scmatrix>(mRes, *this);
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_scmatrix& exp (const basic_scmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp<basic_scmatrix, TR>(*this, mArg, tol);
        return *this;
    }

    basic_scmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_scmatrix mRes (this -> mnM);
        __exp<basic_scmatrix, TR>(mRes, *this, tol);
        return mRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_scmatrix& polynom (const basic_scmatrix& m, const CVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TC, CVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), v.incr() > 1 ? v1 : v);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_scmatrix polynom (const CVector& v) const
    {
        basic_scmatrix mRes (this -> mnM);
        CVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TC, CVector> (mRes.mpD, mRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, v.incr() > 1 ? v1 : v);
        return mRes;
    }

    // eigenvalues
    CVector eig (basic_scmatrix& mEigVect, bool bRightVect = true) const throw (cvmexception)
    {
        CVector vEig (this -> mnM);
        this -> _eig (vEig, &mEigVect, bRightVect);
        return vEig;
    }

    CVector eig() const throw (cvmexception)
    {
        CVector vEig (this -> mnM);
        this -> _eig (vEig, NULL, false);
        return vEig;
    }

    // Cholesky factorization
    basic_scmatrix& cholesky (const basic_schmatrix<TR,TC>& m) throw (cvmexception)
    {
        this -> _check_cholesky();  // doesn't work for band matrices
        *this = m;
        int nOutInfo = __cholesky<basic_scmatrix> (*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        this -> _clean_low_triangle ();
        return *this;
    }

    // Bunch-Kaufman factorization
    basic_scmatrix& bunch_kaufman (const basic_schmatrix<TR,TC>& m, int* nPivots) throw (cvmexception)
    {
        this -> _check_bunch_kaufman();  // doesn't work for band matrices
        *this = m;
        __bunch_kaufman<basic_scmatrix> (*this, nPivots);
        return *this;
    }

    basic_scmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_scmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_scmatrix& randomize_real (TR dFrom, TR dTo)
    {
        this -> _randomize_real (dFrom, dTo);
        return *this;
    }

    basic_scmatrix& randomize_imag (TR dFrom, TR dTo)
    {
        this -> _randomize_imag (dFrom, dTo);
        return *this;
    }

    virtual void _eig (CVector& vEig, basic_scmatrix<TR,TC>* mEigVect, bool bRightVect) const throw (cvmexception)
    {
        __eig<CVector, basic_scmatrix, basic_scmatrix> (vEig, *this, mEigVect, bRightVect);
    }

    virtual void _solve (const CVector& vB, CVector& vX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        CVector vB1;
        CVector vX1;
        if (vB.incr() > 1) vB1 << vB;   // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TC, basic_scmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseCMatrix& mB, BaseCMatrix& mX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TC, basic_scmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

protected:
    // protected constructors for inherited stuff
    basic_scmatrix (int nMN, int nLD, bool bZeroMemory)
        : BaseCMatrix (nMN, nMN, nLD, bZeroMemory)
    {
    }

    basic_scmatrix (TC* pD, int nMN, int nLD, int nSize)
        : BaseCMatrix (pD, nMN, nMN, nLD, nSize)
    {
    }

    virtual int _size   () const {return this -> mnSize;}
    virtual int _msize  () const {return this -> mnM;}
    virtual int _nsize  () const {return this -> mnN;}
    virtual int _ld     () const {return this -> mnLD;}
    virtual const TC* _p() const {return this -> mpD;}
    virtual TC*       _p()       {return this -> mpD;}

    // returns diagonal which IS l-value (shares memory) 
    // 0 - main, negative - low, positive - up
    virtual CVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nD >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        return CVector (this -> mpD + (nDiag > 0 ? nDiag * this -> mnLD : nD), this -> mnM - nD, this -> mnLD + 1);
    }

    // returns diagonal which IS NOT l-value (creates a copy)
    // 0 - main, negative - low, positive - up
    virtual const CVector _diag (int nDiag) const throw (cvmexception)
    {
        CVector vRet (const_cast<basic_scmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    // returns main diagonal of low_up factorization
    virtual CVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _transp()
    {
        this -> _sq_transp();
    }

    virtual void _plus_plus()
    {
        this -> _sq_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _sq_minus_minus();
    }

    virtual TC _det() const throw (cvmexception)
    {
        TC cDet(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                cDet = this -> _ij_val (0, 0);
                break;
            case 2:
                cDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TC one(1.);
                    basic_array<int> naPivots (this -> mnM);
                    CVector vUpDiag = this -> _low_up_diag (naPivots);

                    cDet = one;
                    for (int i = 1; i <= this -> mnM; ++i)
                    {
                        cDet *= vUpDiag[i];
                        if (i != naPivots[i])
                        {
                            cDet = - cDet;
                        }
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_SINGULARMATRIX) throw e;
                }
                break;
        }
        return cDet;
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_scmatrix>(*this, nPivots);
    }

    virtual void _check_cholesky() {}
    virtual void _check_bunch_kaufman() {}
};


// square band matrix
// way of storage:
/*
L=1, U=2          L=2, U=3        L=0, U=0      L=1, U=0

-                 -               d             d
--                --               d            *d
d**               ---               d            *d
*d**              d***               d            *d
 *d**             *d***               d            *d
  *d**            **d***               d            *d
   *d**            **d**                d            _
    *d*             **d*
     *d              **d
      -               --
                       -
*/
template <typename TR, typename TC>
class BandMatrix
{
protected:
    int mnKL;                                                   // number of non-zero sub-diagonals
    int mnKU;                                                   // number of non-zero super-diagonals

    BandMatrix()
        : mnKL(0), mnKU(0)
    {
    }

    BandMatrix (int nKL, int nKU)
        : mnKL (nKL), mnKU (nKU)
    {
    }

    virtual ~BandMatrix()
    {
    }

    virtual int       _size  () const = 0;
    virtual int       _msize () const = 0;
    virtual int       _nsize () const = 0;
    virtual int       _ld    () const = 0;
    virtual const TC* _p     () const = 0;
    virtual       TC* _p     () = 0;
    virtual void      _set_p (TC*) = 0;
    virtual void      _set   (TC* pD, int nSize, int nM, int nN, int nIncr, int nLD) = 0;

    void _bresize (int nNewM) throw (cvmexception)
    {
        TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnSize = this -> _size();
        const int mnIncr = 1;
        if (nNewM != mnM)
        {
            if (nNewM < 0) throw cvmexception (CVM_WRONGSIZE, nNewM);
            const int nNewSize = nNewM * (1 + mnKL + mnKU);

            if (nNewSize == 0)
            {
                cvmFree<TC>(mpD);
                this -> _set (NULL, 0, 0, 0, 0, 0);
                mnKL = mnKU = 0;
            }
            else
            {
                TC* pD = cvmMalloc<TC>(nNewSize);
                if (nNewSize > mnSize) CleanMemory<TC> (pD, nNewSize);
                const int nMinSize = _cvm_min<int>(mnSize, nNewSize);
                __copy<TC> (nMinSize, mpD, mnIncr, pD, 1);
                CVM_ASSERT(pD, nNewSize * sizeof(TC))
                cvmFree<TC>(mpD);
                this -> _set (pD, nNewSize, nNewM, nNewM, 1, 1 + mnKL + mnKU);
            }
        }
    }

    void _check_dimensions () const throw (cvmexception)
    {
        if (mnKL < 0) throw cvmexception (CVM_WRONGSIZE, mnKL);
        if (mnKU < 0) throw cvmexception (CVM_WRONGSIZE, mnKU);
    }

    void _check_dimensions (const BandMatrix& m) const throw (cvmexception)
    {
        if (this -> _msize() != m._msize() || mnKL != m.mnKL || mnKU != m.mnKU) throw cvmexception (CVM_SIZESMISMATCH);
    }

    bool _bcontinuous () const
    {
        return this -> _msize() == 0 || this -> _ld() == 1 + mnKL + mnKU;
    }

    void _mbassign (const Matrix<TR,TC>& m)                               // m is a band matrix here
    {
        TC* mpD = this -> _p();
        const int mnSize = this -> _size();
        if (mpD != m.get())
        {
            __copy<TC> (mnSize, m, m.incr(), mpD, 1);
        }
    }

    void _b_for_each (TC d, bool bRandom, bool bReal, TR dFrom, TR dTo)         // fills the content, doesn't touch head and tail
    {
        int i, jc, ju, jl;
        const int nCol = 1 + mnKL + mnKU;
        TC* mpD = this -> _p();
        const int mnN = this -> _nsize();
        for (int j = 0; j < mnN; ++j)
        {
            jc = j * nCol;
            ju = mnKU - j;
            jl = ju + mnN;
            for (i = 0; i < nCol; ++i)
            {
                if (i >= ju && i < jl)
                {
                    CVM_ASSERT(mpD, (jc + i + 1) * sizeof(TC))
                    if (bRandom)
                    {
                        TR* pD = bReal ? __get_real_p<TR> (mpD + (jc + i)) : __get_imag_p<TR> (mpD + (jc + i));
                        *pD = Randomizer<TR>::get(dFrom, dTo);
                    }
                    else
                    {
                        mpD[jc + i] = d;
                    }
                }
            }
        }
    }

    void _bset (TC d)                                           // fills the content, doesn't touch head and tail
    {
        static const TR zero = TR(0.);
        this -> _b_for_each (d, false, false, zero, zero);
    }

    void _b_randomize (bool bReal, TR dFrom, TR dTo)            // fills the content, doesn't touch head and tail
    {
        static const TC zero(0.);
        this -> _b_for_each (zero, true, bReal, dFrom, dTo);
    }

    TR _bnorm1() const                                          // 1-norm
    {
        int i, j;
        int nLen   = 0;
        int nShift = 0;
        TR  rNorm  = TR(0.);
        TR  rSum;
        const int mnN = this -> _nsize();
        basic_array<TC> rv (this -> _msize());

        for (j = 0; j < mnN; ++j)
        {
            rSum = TR(0.);
            this -> _get_col (j, rv, 1, &nLen, &nShift);

            nLen += nShift;
            for (i = nShift + 1; i <= nLen; ++i)
            {
                rSum += _abs (rv[i]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    TR _bnorminf() const                                        // infinity-norm
    {
        int i, j;
        int nLen = 0;
        int nShift = 0;
        TR  rNorm = TR(0.);
        TR  rSum;
        const int mnM = this -> _msize();
        basic_array<TC> rv (this -> _nsize());

        for (i = 0; i < mnM; ++i)
        {
            rSum = TR(0.);
            this -> _get_row (i, rv, 1, &nLen, &nShift);

            nLen += nShift;
            for (j = nShift + 1; j <= nLen; ++j)
            {
                rSum += _abs (rv[j]);
            }

            if (rSum > rNorm)
            {
                rNorm = rSum;
            }
        }
        return rNorm;
    }

    void _bvanish()
    {
        this -> _bset (TC(0.));
    }

    type_proxy<TC,TR> _b_ij_proxy_val (int i, int j)            // zero based
    {
        static const TC zero = TC(0.);
        TC* mpD = this -> _p();
        const int nA = j - mnKU;
        CVM_ASSERT(mpD, (i + j * (1 + mnKL + mnKU) - nA + 1) * sizeof(TC))
        return (nA < 0 || i >= nA) && i <= mnKL + j ? type_proxy<TC,TR>(mpD [i + j * (1 + mnKL + mnKU) - nA], false) : 
                                                      type_proxy<TC,TR>(zero, true);
    }

    TC _b_ij_val (int i, int j) const                           // zero based
    {
        static const TC zero = TC(0.);
        const TC* mpD = this -> _p();
        const int nA = j - mnKU;
        CVM_ASSERT(mpD, (i + j * (1 + mnKL + mnKU) - nA + 1) * sizeof(TC))
        return (nA < 0 || i >= nA) && i <= mnKL + j ? mpD [i + j * (1 + mnKL + mnKU) - nA] : zero;
    }

    void _get_col (int i, TC* pCol, int nIncr, int* pnLen = NULL, int* pnShift = NULL) const
    {
        const TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnN = this -> _nsize();
        const int nCol = 1 + mnKL + mnKU;
        int nS         = nCol;
        int nShiftSrc  = 0;
        int nShiftDest = 0;

        if (i > mnKU)
        {
            nShiftDest = i - mnKU;
        }
        else
        {
            nShiftSrc = mnKU - i;
            nS -= nShiftSrc;
        }

        if (mnM - i <= mnKL)
        {
            nS -= mnKL + 1 - (mnN - i);
        }

        __copy<TC> (nS,
                   mpD + i * nCol + nShiftSrc,
                   1,
                   pCol + nShiftDest,
                   nIncr);

        if (pnLen) *pnLen = nS;
        if (pnShift) *pnShift = nShiftDest;
    }

    void _get_row (int i, TC* pCol, int nIncr, int* pnLen = NULL, int* pnShift = NULL) const
    {
        const TC* mpD = this -> _p();
        const int mnM = this -> _msize();
        const int mnN = this -> _nsize();
        const int nCol = mnKL + mnKU;
        int nS         = mnN;
        int nShiftSrc  = i + mnKU;
        int nShiftDest = 0;

        if (i > mnKL)
        {
            nShiftDest = i - mnKL;
            nShiftSrc += nShiftDest * nCol;
            nS -= nShiftDest;
        }
        if (mnN - i > mnKU)
        {
            nS -= (mnM - i) - mnKU - 1;
        }

        __copy<TC> (nS,
                   mpD + nShiftSrc,
                   nCol,
                   pCol + nShiftDest,
                   nIncr);

        if (pnLen) *pnLen = nS;
        if (pnShift) *pnShift = nShiftDest;
    }

    void _btransp() throw (cvmexception)
    {
        TC* mpD = this -> _p();
        const int mnN = this -> _nsize();
        if (mnKL > 0 || mnKU > 0)
        {
            const int nLU  = mnKL + mnKU;
            const int nCol = 1 + nLU;
            int nS, nShiftSrc;
            TC* pL;
            TC* pR;
            TC* pD = cvmMalloc<TC>(nCol * mnN);

            for (int i = 0; i < mnN; ++i)
            {
                nS = nCol;
                nShiftSrc = 0;

                if (i < mnKU)
                {
                    nShiftSrc = mnKU - i;
                    nS -= nShiftSrc;
                }
                if (mnN - i <= mnKL)
                {
                    nS -= mnKL + 1 - (mnN - i);
                }

                pL = mpD + i * nCol + nShiftSrc;
                pR = pD + (i > mnKU ? (i - mnKU + 1) * nCol - 1 : mnKL + i);

                __copy<TC> (nS, pL, 1, pR, nLU);
            }

            cvmFree<TC>(mpD);
            this -> _set_p (pD);
            std::swap(mnKL, mnKU);
        }
    }

    void _b_plus_plus()
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnN = this -> _nsize();
        const int nNext = 1 + mnKL + mnKU;
        const int nSize = nNext * mnN;
        CVM_ASSERT(mpD, nSize * sizeof(TC))
        for (int i = mnKU; i < nSize; i += nNext)
        {
            mpD[i] += one;
        }
    }

    void _b_minus_minus()
    {
        TC* mpD = this -> _p();
        static const TC one(1.);
        const int mnN = this -> _nsize();
        const int nNext = 1 + mnKL + mnKU;
        const int nSize = nNext * mnN;
        CVM_ASSERT(mpD, nSize * sizeof(TC))
        for (int i = mnKU; i < nSize; i += nNext)
        {
            mpD[i] -= one;
        }
    }

    void _b_replace (const BandMatrix& m) throw (cvmexception)  // matrix replacement, no assignment
    {
        TC* mpD = this -> _p();
        cvmFree<TC>(mpD);
        int mnSize = m._size();
        mpD = cvmMalloc<TC>(mnSize);
        CVM_ASSERT(mpD, (mnSize * sizeof(TC)))
        mnKL = m.mnKL;
        mnKU = m.mnKU;
        this -> _set (mpD, mnSize, m._msize(), m._nsize(), 1, m._ld());
    }

    void _resize_lu (int nNewKL, int nNewKU) throw (cvmexception)
    {
        if (nNewKL != mnKL || nNewKU != mnKU)
        {
            if (mnKL < 0) throw cvmexception (CVM_WRONGSIZE, mnKL);
            if (mnKU < 0) throw cvmexception (CVM_WRONGSIZE, mnKU);
            TC* mpD = this -> _p();
            const int mnM = this -> _msize();
            const int mnN = this -> _nsize();
            const int nOldLD = 1 + mnKL + mnKU;
            const int nNewLD = 1 + nNewKL + nNewKU;
            const int nMinKL = _cvm_min<int>(mnKL, nNewKL);
            const int nMinKU = _cvm_min<int>(mnKU, nNewKU);
            const int nNewSize = mnN * (1 + nNewKL + nNewKU);
            TC* pD = cvmMalloc<TC>(nNewSize);
            CVM_ASSERT(pD, nNewSize * sizeof(TC))
            CleanMemory<TC> (pD, nNewSize);

            for (int i = - nMinKU; i <= nMinKL; ++i)
            {
                __copy<TC> (mnN, mpD + (mnKU + i), nOldLD, pD + (nNewKU + i), nNewLD);
            }

            cvmFree<TC>(mpD);
            mnKL = nNewKL;
            mnKU = nNewKU;
            this -> _set (pD, nNewSize, mnM, mnN, 1, nNewLD);
        }
    }

public:
    int lsize() const
    {
        return mnKL;
    }

    int usize() const
    {
        return mnKU;
    }

    const int* _pl() const
    {
        return &mnKL;
    }

    const int* _pu() const
    {
        return &mnKU;
    }
};


template <typename TR>
class basic_srbmatrix : public basic_srmatrix<TR>, public BandMatrix<TR,TR>
{
    typedef std::complex<TR> TC;
    typedef basic_rvector<TR> RVector;
    typedef basic_cvector<TR,TC> CVector;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;
    typedef basic_srmatrix<TR> BaseSRMatrix;
    typedef BandMatrix<TR,TR> BaseBandMatrix;

    friend class basic_scbmatrix<TR,TC>;  // basic_scbmatrix constructor

protected:
    mutable BaseSRMatrix mSM;

    void _bake_SM() const
    {
        mSM.resize(this -> mnM);
        mSM.vanish();
        _copy_b_matrix<TR, TR, BaseSRMatrix, basic_srbmatrix>(mSM, * const_cast<basic_srbmatrix*>(this), false);
    }

public:
    basic_srbmatrix()
    {
    }

    explicit basic_srbmatrix (int nMN)
        : BaseSRMatrix (nMN, 1, true), BaseBandMatrix (0, 0)                    // diagonal matrix
    {
    }

    basic_srbmatrix (int nMN, int nKL, int nKU)
        : BaseSRMatrix (nMN, 1 + nKL + nKU, true), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_srbmatrix (TR* pD, int nMN, int nKL, int nKU)
        : BaseSRMatrix (pD, nMN, 1 + nKL + nKU, nMN * (1 + nKL + nKU)), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_srbmatrix (const basic_srbmatrix& m)
        : BaseSRMatrix (m.mnM, 1 + m.mnKL + m.mnKU, false), BaseBandMatrix (m.mnKL, m.mnKU)
    {
        this -> _massign(m);
    }

    basic_srbmatrix (const BaseRMatrix& m, int nKL, int nKU)
        : BaseSRMatrix (m.msize(), 1 + nKL + nKU, false), BaseBandMatrix (nKL, nKU)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _check_dimensions();
        _copy_b_matrix<TR, TR, BaseRMatrix, basic_srbmatrix> (const_cast<BaseRMatrix&>(m), *this, true);
    }

    // diagonal square matrix constructor
    explicit basic_srbmatrix (const RVector& v)
        : BaseSRMatrix (v.size(), 1, false), BaseBandMatrix (0, 0)
    {
        __copy<TR> (this -> mnM, v, v.incr(), this -> mpD, 1);
    }

    type_proxy<TR,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    // returns row which CAN NOT be l-value
    const RVector operator [] (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    basic_srbmatrix& operator = (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions(m);
        this -> _massign(m);
        return *this;
    }

    // assigns foregn array (nIncr = 1)
    basic_srbmatrix& assign (const TR* pD)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    // fills the content
    basic_srbmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    basic_srbmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    basic_srbmatrix& resize_lu (int nNewKL, int nNewKU) throw (cvmexception)
    {
        this -> _resize_lu (nNewKL, nNewKU);
        return *this;
    }

    bool operator == (const basic_srbmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> mnKL == m.mnKL && this -> mnKU == m.mnKU && this -> _mequals (m);
    }

    bool operator != (const basic_srbmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_srbmatrix& operator << (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_ld();                                    // submatrix replacement is obviously not possible
        this -> _b_replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srbmatrix operator + (const basic_srbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_srbmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srbmatrix operator - (const basic_srbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_srbmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_srbmatrix& sum (const basic_srbmatrix& m1, const basic_srbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_srbmatrix& diff (const basic_srbmatrix& m1, const basic_srbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srbmatrix& operator += (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mincr (m);
        return *this;
    }

    basic_srbmatrix& operator -= (const basic_srbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mdecr (m);
        return *this;
    }

    basic_srbmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srbmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srbmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srbmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srbmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srbmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srbmatrix operator * (TR dMult) const
    {
        basic_srbmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srbmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srbmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srbmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srbmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix
    basic_srbmatrix operator ~ () const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        return mRes.transpose();
    }

    // well, not the best possible algorithm, has to be optimized
    basic_srbmatrix& transpose (const basic_srbmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnKL != m.mnKU || this -> mnKU != m.mnKL) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srbmatrix mTmp(m);
        mTmp.transpose();
        __copy<TR> (this -> mnSize, mTmp, mTmp.incr(), this -> mpD, this -> mnIncr);
        return *this;
    }

    basic_srbmatrix& transpose() throw (cvmexception)
    {
        this -> _transp();
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    BaseSRMatrix operator * (const BaseSRMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    basic_srbmatrix operator * (const basic_srbmatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        m._bake_SM();
        return basic_srbmatrix (mSM * m.mSM, 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKL + m.mnKL), 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKU + m.mnKU));
    }

    // low-up factorization
    basic_srbmatrix& low_up (const basic_srbmatrix& m, int* nPivots) throw (cvmexception)
    {
        (*this) = m;
        this -> _low_up (nPivots);
        return *this;
    }

    basic_srbmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_srbmatrix mRes (*this);
        mRes._low_up (nPivots);
        return mRes;
    }

    basic_srbmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srbmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_srbmatrix& randomize (TR dFrom, TR dTo)
    {
        this -> _b_randomize (true, dFrom, dTo);
        return *this;
    }

    // ATTENTION!!! This is not a good idea to use the following function. It's provided for
    // overriding of basic_rmatrix<TR>::mult only
    // this = m1 * m2
    basic_srbmatrix& mult (const BaseRMatrix& m1, const BaseRMatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    virtual TR* _pd()                                           // redefinition of Array's function
    {
#ifdef CVM_DEBUG
        assert (false);     // it's abnormal to call this function, this is pointer to copy, not to an object. only const version is allowable
#endif
        _bake_SM();
        return mSM._pd();
    }

    virtual const TR* _pd() const                               // redefinition of Array's function
    {
        _bake_SM();
        return mSM._pd();
    }

    // Euclid norm - band matrices require special care because of tails
    virtual TR norm() const throw (cvmexception)
    {
        TR dNorm = TR(0.), d;
        int i;

        for (i = 0; i <= this -> mnKL; ++i)
        {
            const RVector& dV = const_cast<basic_srbmatrix*>(this) -> diag(-i);
            d = dV.norm();
            dNorm += d * d;
        }
        for (i = 1; i <= this -> mnKU; ++i)
        {
            const RVector& dV = const_cast<basic_srbmatrix*>(this) -> diag(i);
            d = dV.norm();
            dNorm += d * d;
        }

        return _sqrt(dNorm);
    }

    virtual TR norm1() const
    {
        return this -> _bnorm1();
    }

    virtual TR norminf() const
    {
        return this -> _bnorminf();
    }

    // singular values in decreasing order
    virtual void _svd (RVector& vRes, BaseSRMatrix* pmU, BaseSRMatrix* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_srbmatrix, BaseSRMatrix> (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (BaseRMatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_srbmatrix, BaseRMatrix> (mX, *this, threshold);
    }

    virtual void _eig (CVector& vEig, basic_scmatrix<TR,TC>* mEigVect, bool bRightVect) const throw (cvmexception)
    {
        _bake_SM();
        mSM._eig (vEig, mEigVect, bRightVect);
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        RVector vB1;
        RVector vX1;
        if (vB.incr() > 1) vB1 << vB;       // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TR, basic_srbmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TR, basic_srbmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

    // ?gbmv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gbmv calls all virtual methods inside
    void _gbmv (bool bLeft, TR dAlpha, const RVector& v, TR dBeta, RVector& vRes) const
    {
        RVector vTmp;
        basic_srbmatrix mTmp;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gbmv<TR, basic_srbmatrix, RVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha, 
                                                     vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    virtual void _check_submatrix ()  throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "srbmatrix");}

protected:
    virtual int       _size  () const { return this -> mnSize; }
    virtual int       _msize () const { return this -> mnM; }
    virtual int       _nsize () const { return this -> mnN; }
    virtual int       _ld    () const { return this -> mnLD; }
    virtual const TR* _p     () const { return this -> mpD; }
    virtual       TR* _p     ()       { return this -> mpD; }

    virtual void _set_p (TR* pD)
    {
        this -> mpD = pD;
    }

    virtual void _set (TR* pD, int nSize, int nM, int nN, int nIncr, int nLD)
    {
        this -> mpD = pD;
        this -> mnSize = nSize;
        this -> mnM = nM;
        this -> mnN = nN;
        this -> mnIncr = nIncr;
        this -> mnLD = nLD;
    }

    // for _msum _mdiff etc.
    virtual const TR* _p (const BaseMatrix& m) const
    {
        return m.get();
    }

    virtual int _ldm() const
    {
        return this -> mnM;
    }

    virtual const int* _pldm() const
    {
        return &this -> mnM;
    }

    // 0-based
    virtual RVector _row (int m)
    {
        RVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual RVector _col (int n)
    {
        RVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const RVector _row (int m) const
    {
        RVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const RVector _col (int n) const
    {
        RVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    virtual RVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nDiag < 0 && nD > this -> mnKL || nDiag > 0 && nDiag > this -> mnKU) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        const int nLU = this -> mnKL + this -> mnKU;
        return RVector (this -> mpD + this -> mnKU + (nDiag < 0 ? nD : nD * nLU), this -> mnM - nD, nLU + 1);
    }

    virtual const RVector _diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (const_cast<basic_srbmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    virtual void _swap_rows(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "srbmatrix");}
    virtual void _swap_cols(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "srbmatrix");}
    virtual void _check_ger()           throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "ger", "srbmatrix");}
    virtual void _check_rank1update()   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update", "srbmatrix");}
    virtual void _check_gemm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "srbmatrix");}
    virtual void _check_symm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "symm", "srbmatrix");}
    virtual void _check_cholesky()      throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "cholesky", "srbmatrix");}
    virtual void _check_bunch_kaufman() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "bunch_kaufman", "srbmatrix");}

    virtual void _resize (int nNewM, int) throw (cvmexception)
    {
        this -> _bresize(nNewM);
    }

    virtual bool _continuous () const
    {
        return this -> _bcontinuous ();
    }

    virtual void _massign (const BaseMatrix& m)
    {
        this -> _mbassign(m);
    }

    virtual void _set (TR d)
    {
        this -> _bset(d);
    }

    virtual void _vanish()
    {
        this -> _bvanish();
    }

    // zero based
    virtual type_proxy<TR,TR> _ij_proxy_val (int i, int j)
    {
        return this -> _b_ij_proxy_val (i, j);
    }

    // zero based
    virtual TR _ij_val (int i, int j) const
    {
        return this -> _b_ij_val (i, j);
    }

    virtual void _transp() throw (cvmexception)
    {
        this -> _btransp();
    }

    virtual void _plus_plus()
    {
        this -> _b_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _b_minus_minus();
    }

    virtual int _indofmax() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmax();
    }

    virtual int _indofmin() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmin();
    }

    // returns main diagonal of low_up factorization
    virtual RVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, d);                // zero tails are supposed here
    }

    virtual void _mult (const BaseRMatrix& m1, const BaseRMatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.msize() || this -> mnN != m2.nsize() || m1.nsize() != m2.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseSRMatrix mR (this -> mnM);
        mR.mult (m1, m2);
        this -> _resize_lu (this -> mnM - 1, this -> mnM - 1);
        _copy_b_matrix<TR, TR, BaseSRMatrix, basic_srbmatrix> (const_cast<BaseSRMatrix&>(mR), *this, true);
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool bLeft) const
    {
        static const TR zero(0.);
        static const TR one(1.);
        this -> _gbmv (bLeft, one, v, zero, vRes);
    }

    virtual void _randomize (TR dFrom, TR dTo)
    {
        __randomize<TR>(this -> mpD, this -> mnSize, this -> mnIncr, dFrom, dTo);
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_srbmatrix>(*this, nPivots);
    }

    virtual int _ld_for_replace () const
    {
        return this -> mnM;
    }

    virtual int _size_for_replace () const
    {
        return this -> mnM * this -> mnN;
    }
};


template <typename TR, typename TC>
class basic_scbmatrix : public basic_scmatrix<TR,TC>, public BandMatrix<TR,TC>
{
    typedef basic_cvector<TR,TC> CVector;
    typedef Array<TR,TC> BaseArray;
    typedef Matrix<TR,TC> BaseMatrix;
    typedef basic_cmatrix<TR,TC> BaseCMatrix;
    typedef basic_scmatrix<TR,TC> BaseSCMatrix;
    typedef BandMatrix<TR,TC> BaseBandMatrix;

protected:
    mutable BaseSCMatrix mSM;

    void _bake_SM() const
    {
        mSM.resize(this -> mnM);
        mSM.vanish();
        _copy_b_matrix<TR, TC, BaseSCMatrix, basic_scbmatrix>(mSM, * const_cast<basic_scbmatrix*>(this), false);
    }

public:
    virtual TC* _pd()                                           // redefinition of Array's function
    {
#ifdef CVM_DEBUG
        // it's abnormal to call this function, this is pointer to copy, not to an object. only const version is allowable
        assert (false);
#endif
        _bake_SM();
        return mSM._pd();
    }

    virtual const TC* _pd() const                               // redefinition of Array's function
    {
        _bake_SM();
        return mSM._pd();
    }

    basic_scbmatrix()
    {
    }

    // diagonal matrix
    explicit basic_scbmatrix (int nMN)
        : BaseSCMatrix (nMN, 1, true), BaseBandMatrix (0, 0)
    {
    }

    basic_scbmatrix (int nMN, int nKL, int nKU)
        : BaseSCMatrix (nMN, 1 + nKL + nKU, true), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_scbmatrix (TC* pD, int nMN, int nKL, int nKU)
        : BaseSCMatrix (pD, nMN, 1 + nKL + nKU, nMN * (1 + nKL + nKU)), BaseBandMatrix (nKL, nKU)
    {
        this -> _check_dimensions();
    }

    basic_scbmatrix (const basic_scbmatrix& m)
        : BaseSCMatrix (m.mnM, 1 + m.mnKL + m.mnKU, false), BaseBandMatrix (m.mnKL, m.mnKU)
    {
        this -> _massign(m);
    }

    basic_scbmatrix (const BaseCMatrix& m, int nKL, int nKU)
        : BaseSCMatrix (m.msize(), 1 + nKL + nKU, false), BaseBandMatrix (nKL, nKU)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _check_dimensions();
        _copy_b_matrix<TR, TC, BaseCMatrix, basic_scbmatrix> (const_cast<BaseCMatrix&>(m), *this, true);
    }

    // diagonal square matrix constructor
    explicit basic_scbmatrix (const CVector& v)
        : BaseSCMatrix (v.size(), 1, false), BaseBandMatrix (0, 0)
    {
        __copy<TC> (this -> mnM, v, v.incr(), this -> mpD, 1);
    }

    explicit basic_scbmatrix (const basic_srbmatrix<TR>& m, bool bRealPart = true)
        : BaseSCMatrix (m.msize(), m.ld(), true), BaseBandMatrix (m.lsize(), m.usize())
    {
        if (bRealPart)
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, m._p(), NULL);
        else
            __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, NULL, m._p());
    }

    basic_scbmatrix (const basic_srbmatrix<TR>& mRe, const basic_srbmatrix<TR>& mIm)
        : BaseSCMatrix (mRe.msize(), mRe.ld(), false), BaseBandMatrix (mRe.lsize(), mRe.usize())
    {
        if (mRe.msize() != mIm.msize() || mRe.nsize() != mIm.nsize() ||
            mRe.lsize() != mIm.lsize() || mRe.usize() != mIm.usize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe, mIm, mRe.incr(), mIm.incr());
    }

    type_proxy<TC,TR> operator () (int nIm, int nIn) throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_proxy_val (nIm - 1, nIn - 1);
    }

    TC operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const CVector operator () (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnN) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _col (nI - 1);
    }

    // returns row which CAN NOT be l-value
    const CVector operator [] (int nI) const throw (cvmexception)
    {
        if (nI <= 0 || nI > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nI);
        return this -> _row (nI - 1);
    }

    // real part
    const basic_srbmatrix<TR> real() const
    {
        basic_srbmatrix<TR> mRet (this -> mnM, this -> mnKL, this -> mnKU);
        __copy<TR> (this -> mnSize, __get_real_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    // imaginary part
    const basic_srbmatrix<TR> imag() const
    {
        basic_srbmatrix<TR> mRet (this -> mnM, this -> mnKL, this -> mnKU);
        __copy<TR> (this -> mnSize, __get_imag_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    basic_scbmatrix& operator = (const basic_scbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions(m);
        this -> _massign(m);
        return *this;
    }

    // assigns a foregn array (nIncr = 1)
    basic_scbmatrix& assign (const TC* pD)
    {
        this -> _assign (pD, 1);
        return *this;
    }

    basic_scbmatrix& set (TC c)
    {
        this -> _set (c);
        return *this;
    }

    // assigns real matrix
    basic_scbmatrix& assign_real (const basic_srbmatrix<TR>& mRe) throw (cvmexception)
    {
        if (this -> mnM != mRe.msize() || this -> mnKL != mRe.lsize() || this -> mnKU != mRe.usize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe, mRe.incr());
        return *this;
    }

    // assigns imaginary array
    basic_scbmatrix& assign_imag (const basic_srbmatrix<TR>& mIm) throw (cvmexception)
    {
        if (this -> mnM != mIm.msize() || this -> mnKL != mIm.lsize() || this -> mnKU != mIm.usize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_imag<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mIm, mIm.incr());
        return *this;
    }

    // fills real part
    basic_scbmatrix& set_real (TR d)
    {
        this -> _set_real_number (d);
        return *this;
    }

    // fills imaginary part
    basic_scbmatrix& set_imag (TR d)
    {
        this -> _set_imag_number (d);
        return *this;
    }

    basic_scbmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    basic_scbmatrix& resize_lu (int nNewKL, int nNewKU) throw (cvmexception)
    {
        this -> _resize_lu (nNewKL, nNewKU);
        return *this;
    }

    bool operator == (const basic_scbmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> mnKL == m.mnKL && this -> mnKU == m.mnKU && this -> _mequals (m);
    }

    bool operator != (const basic_scbmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_scbmatrix& operator << (const basic_scbmatrix& m) throw (cvmexception)
    {
        this -> _check_ld();                                    // submatrix replacement is obviously not possible
        this -> _b_replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_scbmatrix operator + (const basic_scbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_scbmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_scbmatrix operator - (const basic_scbmatrix& m) const throw (cvmexception)
    {
        this -> _check_dimensions (m);
        basic_scbmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_scbmatrix& sum (const basic_scbmatrix& m1, const basic_scbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_scbmatrix& diff (const basic_scbmatrix& m1, const basic_scbmatrix& m2) throw (cvmexception)
    {
        this -> _check_dimensions (m1);
        this -> _check_dimensions (m2);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_scbmatrix& operator += (const basic_scbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mincr (m);
        return *this;
    }

    basic_scbmatrix& operator -= (const basic_scbmatrix& m) throw (cvmexception)
    {
        this -> _check_dimensions (m);
        this -> _mdecr (m);
        return *this;
    }

    basic_scbmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_scbmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_scbmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_scbmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_scbmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_scbmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_scbmatrix operator * (TR dMult) const
    {
        basic_scbmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_scbmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_scbmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_scbmatrix operator * (TC cMult) const
    {
        basic_scbmatrix mRes (*this);
        mRes._scal (cMult);
        return mRes;
    }

    basic_scbmatrix operator / (TC cDiv) const throw (cvmexception)
    {
        basic_scbmatrix mRes (*this);
        mRes._div (cDiv);
        return mRes;
    }

    basic_scbmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_scbmatrix& operator /= (TR dDiv)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_scbmatrix& operator *= (TC cMult)
    {
        this -> _scal (cMult);
        return *this;
    }

    basic_scbmatrix& operator /= (TC cDiv) throw (cvmexception)
    {
        this -> _div (cDiv);
        return *this;
    }

    basic_scbmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix
    basic_scbmatrix operator ~ () const throw (cvmexception)
    {
        basic_scbmatrix mRes (*this);
        return mRes.conj();
    }

    // well, not the best possible algorithm, has to be optimized
    basic_scbmatrix& conj (const basic_scbmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM || this -> mnKL != m.mnKU || this -> mnKU != m.mnKL) throw cvmexception (CVM_SIZESMISMATCH);
        basic_scbmatrix mTmp(m);
        mTmp.conj();
        __copy<TC> (this -> mnSize, mTmp, mTmp.incr(), this -> mpD, this -> mnIncr);
        return *this;
    }

    basic_scbmatrix& conj() throw (cvmexception)
    {
        this -> _transp();
        __conj<TC> (this -> mpD, this -> mnSize, this -> mnIncr);
        return *this;
    }

    CVector operator * (const CVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseCMatrix operator * (const BaseCMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    BaseSCMatrix operator * (const BaseSCMatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        return mSM * m;
    }

    basic_scbmatrix operator * (const basic_scbmatrix& m) const throw (cvmexception)
    {
        _bake_SM();
        m._bake_SM();
        return basic_scbmatrix (mSM * m.mSM, 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKL + m.mnKL), 
                                _cvm_min<int>(this -> mnM - 1, this -> mnKU + m.mnKU));
    }

    // low-up factorization
    basic_scbmatrix& low_up (const basic_scbmatrix& m, int* nPivots) throw (cvmexception)
    {
        *this = m;
        this -> _low_up (nPivots);
        return *this;
    }

    basic_scbmatrix low_up (int* nPivots) const throw (cvmexception)
    {
        basic_scbmatrix mRes (*this);
        mRes._low_up (nPivots);
        return mRes;
    }

    basic_scbmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_scbmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_scbmatrix& randomize_real (TR dFrom, TR dTo)
    {
        this -> _b_randomize (true, dFrom, dTo);
        return *this;
    }

    basic_scbmatrix& randomize_imag (TR dFrom, TR dTo)
    {
        this -> _b_randomize (false, dFrom, dTo);
        return *this;
    }

    // ATTENTION!!! This is not a good idea to use the following function. It's provided for
    // overriding of basic_rmatrix<TR>::mult only
    basic_scbmatrix& mult (const BaseCMatrix& m1, const BaseCMatrix& m2) throw (cvmexception)
    {
        this -> _mult (m1, m2);
        return *this;
    }

    // Euclid norm - band matrices require special care because of tails
    virtual TR norm() const throw (cvmexception)
    {
        TR dNorm = TR(0.), d;
        int i;
        for (i = 0; i <= this -> mnKL; ++i)
        {
            const CVector& dV = const_cast<basic_scbmatrix*>(this) -> diag(-i);
            d = dV.norm();
            dNorm += d * d;
        }
        for (i = 1; i <= this -> mnKU; ++i)
        {
            const CVector& dV = const_cast<basic_scbmatrix*>(this) -> diag(i);
            d = dV.norm();
            dNorm += d * d;
        }
        return _sqrt(dNorm);
    }

    virtual TR norm1() const
    {
        return this -> _bnorm1();
    }

    virtual TR norminf() const
    {
        return this -> _bnorminf();
    }

    // singular values in decreasing order
    virtual void _svd (basic_rvector<TR>& vRes, BaseSCMatrix* pmU, BaseSCMatrix* pmVH) const throw (cvmexception)
    {
        if (pmU != NULL && pmVH != NULL && (this -> mnM != pmU -> msize() || this -> mnN != pmVH -> msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __svd<TR, basic_scbmatrix, BaseSCMatrix> (vRes, vRes.size(), vRes.incr(), *this, pmU, pmVH);
    }

    virtual void _pinv (BaseCMatrix& mX, TR threshold) const throw (cvmexception)
    {
        __pinv<TR, basic_scbmatrix, BaseCMatrix> (mX, *this, threshold);
    }

    virtual void _eig (CVector& vEig, basic_scmatrix<TR,TC>* mEigVect, bool bRightVect) const throw (cvmexception)
    {
        _bake_SM();
        mSM._eig (vEig, mEigVect, bRightVect);
    }

    virtual void _solve (const CVector& vB, CVector& vX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        vX = vB;
        CVector vB1;
        CVector vX1;
        if (vB.incr() > 1) vB1 << vB;       // to make sure incr = 1
        if (vX.incr() > 1) vX1 << vX;
        __solve<TR, TC, basic_scbmatrix> (*this, 1, vB.incr() > 1 ? vB1 : vB, vB.size(), vX.incr() > 1 ? vX1 : vX, vX.size(), dErr, pLU, pPivots);
        if (vX.incr() > 1) vX = vX1;
    }

    virtual void _solve (const BaseCMatrix& mB, BaseCMatrix& mX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        mX = mB;
        __solve<TR, TC, basic_scbmatrix> (*this, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
    }

    // ?gbmv routines perform a matrix-vector operation defined as
    // vRes = alpha*m*v + beta * vRes or vRes = alpha*v'*m + beta * vRes
    // not virtual since __gbmv calls all virtual methods inside
    void _gbmv (bool bLeft, TC dAlpha, const CVector& v, TC dBeta, CVector& vRes) const
    {
        CVector vTmp;
        basic_scbmatrix mTmp;
        const TC* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == this -> mpD) mTmp << *this;
        __gbmv<TC, basic_scbmatrix, CVector> (bLeft, vRes.get() == this -> mpD ? mTmp : *this, dAlpha,
                                                     vRes.get() == pDv ? vTmp : v, dBeta, vRes);
    }

    virtual void _check_submatrix ()    throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "scbmatrix");}

protected:
    virtual int       _size  () const { return this -> mnSize; }
    virtual int       _msize () const { return this -> mnM; }
    virtual int       _nsize () const { return this -> mnN; }
    virtual int       _ld    () const { return this -> mnLD; }
    virtual const TC* _p     () const { return this -> mpD; }
    virtual       TC* _p     ()       { return this -> mpD; }

    virtual void _set_p (TC* pD)
    {
        this -> mpD = pD;
    }

    virtual void _set (TC* pD, int nSize, int nM, int nN, int nIncr, int nLD)
    {
        this -> mpD = pD;
        this -> mnSize = nSize;
        this -> mnM = nM;
        this -> mnN = nN;
        this -> mnIncr = nIncr;
        this -> mnLD = nLD;
    }

    virtual const TC* _p (const BaseMatrix& m) const // for _msum _mdiff etc.
    {
        return m.get();
    }

    virtual int _ldm() const
    {
        return this -> mnM;
    }

    virtual const int* _pldm() const
    {
        return &this -> mnM;
    }

    // 0-based, returns l-value sharing memory
    virtual CVector _row (int m)
    {
        CVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based, returns l-value sharing memory
    virtual CVector _col (int n)
    {
        CVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const CVector _row (int m) const
    {
        CVector vRet (this -> mnM);
        this -> _get_row (m, vRet, vRet.incr());
        return vRet;
    }

    // 0-based
    virtual const CVector _col (int n) const
    {
        CVector vRet (this -> mnM);
        this -> _get_col (n, vRet, vRet.incr());
        return vRet;
    }

    virtual CVector _diag (int nDiag) throw (cvmexception)
    {
        const int nD = abs (nDiag);
        if (nDiag < 0 && nD > this -> mnKL || nDiag > 0 && nDiag > this -> mnKU) throw cvmexception (CVM_OUTOFRANGE, nDiag);
        const int nLU = this -> mnKL + this -> mnKU;
        return CVector (this -> mpD + this -> mnKU + (nDiag < 0 ? nD : nD * nLU), this -> mnM - nD, nLU + 1);
    }

    virtual const CVector _diag (int nDiag) const throw (cvmexception)
    {
        CVector vRet (const_cast<basic_scbmatrix*>(this) -> _diag(nDiag));
        return vRet;
    }

    virtual void _resize (int nNewM, int) throw (cvmexception)
    {
        this ->_bresize(nNewM);
    }

    virtual bool _continuous () const
    {
        return this -> _bcontinuous ();
    }

    virtual void _massign (const BaseMatrix& m)
    {
        this -> _mbassign(m);
    }

    virtual void _set (TC d)
    {
        this -> _bset(d);
    }

    virtual void _vanish()
    {
        this -> _bvanish();
    }

    // zero based
    virtual type_proxy<TC,TR> _ij_proxy_val (int i, int j)
    {
        return this -> _b_ij_proxy_val (i, j);
    }

    // zero based
    virtual TC _ij_val (int i, int j) const
    {
        return this -> _b_ij_val (i, j);
    }

    virtual void _transp() throw (cvmexception)
    {
        this -> _btransp();
    }

    virtual void _plus_plus()
    {
        this -> _b_plus_plus();
    }

    virtual void _minus_minus()
    {
        this -> _b_minus_minus();
    }

    virtual int _indofmax() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmax();
    }

    virtual int _indofmin() const
    {
        this -> _check_ld();
        _bake_SM();
        return mSM.indofmin();
    }

    virtual CVector _low_up_diag (basic_array<int>& naPivots) const throw (cvmexception)
    {
        return this -> low_up (naPivots).diag(0);
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);     // zero tails are supposed here
    }

    virtual void _scal (TC d)
    {
        __scal<TC, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);     // zero tails are supposed here
    }

    virtual void _mult (const BaseCMatrix& m1, const BaseCMatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.msize() || this -> mnN != m2.nsize() || m1.nsize() != m2.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        BaseSCMatrix mR (this -> mnM);
        mR.mult (m1, m2);
        this -> _resize_lu (this -> mnM - 1, this -> mnM - 1);
        _copy_b_matrix<TR, TC, BaseSCMatrix, basic_scbmatrix> (const_cast<BaseSCMatrix&>(mR), *this, true);
    }

    virtual void _multiply (CVector& vRes, const CVector& v, bool bLeft) const
    {
        static const TC zero = TC(0.);
        static const TC one = TC(1.);
        this -> _gbmv (bLeft, one, v, zero, vRes);
    }

    virtual void _low_up (int* nPivots) throw (cvmexception)
    {
        __low_up<basic_scbmatrix>(*this, nPivots);
    }

    virtual int _ld_for_replace () const
    {
        return this -> mnM;
    }

    virtual int _size_for_replace () const
    {
        return this -> mnM * this -> mnN;
    }

    virtual void _swap_rows(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "scbmatrix");}
    virtual void _swap_cols(int, int)   throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "scbmatrix");}
    virtual void _check_geru()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "geru", "scbmatrix");}
    virtual void _check_gerc()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gerc", "scbmatrix");}
    virtual void _check_rank1update_u() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update_u", "scbmatrix");}
    virtual void _check_rank1update_c() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update_c", "scbmatrix");}
    virtual void _check_gemm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "scbmatrix");}
    virtual void _check_hemm()          throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "hemm", "scbmatrix");}
    virtual void _check_cholesky()      throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "cholesky", "scbmatrix");}
    virtual void _check_bunch_kaufman() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "bunch_kaufman", "scbmatrix");}
};


template <typename TR>
class basic_srsmatrix : public basic_srmatrix<TR>
{
    typedef basic_rvector<TR> RVector;
    typedef Array<TR,TR> BaseArray;
    typedef Matrix<TR,TR> BaseMatrix;
    typedef SqMatrix<TR, TR> BaseSqMatrix;
    typedef basic_rmatrix<TR> BaseRMatrix;
    typedef basic_srmatrix<TR> BaseSRMatrix;

public:
    basic_srsmatrix()
    {
    }

    explicit basic_srsmatrix (int nMN)
        : BaseSRMatrix (nMN)
    {
    }

    basic_srsmatrix (TR* pD, int nMN)
        : BaseSRMatrix (pD, nMN, nMN, nMN * nMN)
    {
        this -> _check_symmetric();
    }

    basic_srsmatrix (const basic_srsmatrix& m)
        : BaseSRMatrix (m.mnM, m.mnM, false)
    {
        this -> _massign(m);
    }

    basic_srsmatrix (const BaseRMatrix& m)
        : BaseSRMatrix (m.msize(), m.nsize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        this -> _check_symmetric();
    }

    // diagonal square symmetric matrix constructor
    explicit basic_srsmatrix (const RVector& v)
        : BaseSRMatrix (v)
    {
    }

    // submatrix constructor
    basic_srsmatrix (basic_srsmatrix& m, int nRowCol, int nSize)        // 1-based
        : BaseSRMatrix (m._sub_pointer (nRowCol, nRowCol, nSize, nSize), nSize, m.ld(), nSize * nSize)
    {
    }

    TR operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const RVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns row which CAN NOT be l-value
    const RVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns diagonal which IS NOT l-value (since it could break symmetricity)
    // 0 - main, negative - low, positive - up
    const RVector diag (int nDiag) const throw (cvmexception)
    {
        RVector vRet (this -> _diag(nDiag));
        return vRet;
    }

    basic_srsmatrix& operator = (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    // assigns a foregn array (nIncr = 1)
    basic_srsmatrix& assign (const TR* pD) throw (cvmexception)
    {
        this -> _assign (pD, 1);
        this -> _check_symmetric();
        return *this;
    }

    basic_srsmatrix& assign (int nRowCol, const basic_srsmatrix& m) throw (cvmexception)         // subvector assignment
    {
        if (nRowCol <= 0 || nRowCol > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRowCol);
        if (m.mnM + nRowCol - 1 > this -> mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRowCol, nRowCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    basic_srsmatrix& set (TR d)
    {
        this -> _set (d);
        return *this;
    }

    // sets both elements to keep matrix symmetric
    basic_srsmatrix& set (int nRow, int nCol, TR d)
    {
        this -> _set_at (nRow - 1, nCol - 1, d);
        return *this;
    }

    // sets both diagonals (keeps matrix symmetric)
    basic_srsmatrix& set_diag (int nDiag, const RVector& vDiag) throw (cvmexception)
    {
        RVector (this -> _diag(nDiag)) = vDiag;
        if (nDiag != 0)
        {
            RVector (this -> _diag(-nDiag)) = vDiag;
        }
        return *this;
    }

    basic_srsmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    bool operator == (const basic_srsmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_srsmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_srsmatrix& operator << (const basic_srsmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_srsmatrix operator + (const basic_srsmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srsmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_srsmatrix operator - (const basic_srsmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_srsmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_srsmatrix& sum (const basic_srsmatrix& m1, const basic_srsmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_srsmatrix& diff (const basic_srsmatrix& m1, const basic_srsmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_srsmatrix& operator += (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_srsmatrix& operator -= (const basic_srsmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_srsmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_srsmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_srsmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_srsmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_srsmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_srsmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_srsmatrix operator * (TR dMult) const
    {
        basic_srsmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_srsmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_srsmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    basic_srsmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_srsmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_srsmatrix& normalize()
    {
        this -> _normalize();
        return *this;
    }

    // transposed Matrix - does nothing
    basic_srsmatrix operator ~ () const
    {
        return *this;
    }

    basic_srsmatrix& transpose (const basic_srsmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return *this;
    }

    basic_srsmatrix& transpose()
    {
        return *this;
    }

    RVector operator * (const RVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseRMatrix operator * (const BaseRMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseRMatrix mRes(m.msize(), m.nsize());
        mRes._symm (true, *this, m, one, zero);
        return mRes;
    }

    BaseSRMatrix operator * (const BaseSRMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseSRMatrix mRes(this -> mnM);
        mRes._symm (true, *this, m, one, zero);
        return mRes;
    }

    // this = alpha*v*v' + beta*this
    basic_srsmatrix& syrk (TR alpha, const RVector& v, TR beta) throw (cvmexception)
    {
        if (this -> mnM != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __syrk<TR, basic_srsmatrix> (false, alpha, 1, v, v.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*a' + beta*this or this = alpha*a'*a + beta*this
    basic_srsmatrix& syrk (bool bTransp, TR alpha, const BaseRMatrix& mA, TR beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __syrk<TR, basic_srsmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*v1*v2' + alpha*v2*v1' + beta*this
    basic_srsmatrix& syr2k (TR alpha, const RVector& v1, const RVector& v2, TR beta) throw (cvmexception)
    {
        if (this -> mnM != v1.size() || this -> mnM != v2.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __syr2k<TR, basic_srsmatrix> (false, alpha, 1, v1, v1.size(), v2, v2.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*b' + alpha*b*a' + beta*this or this = alpha*a'*b + alpha*b'*a + beta*this
    basic_srsmatrix& syr2k (bool bTransp, TR alpha, const BaseRMatrix& mA, const BaseRMatrix& mB, TR beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize()) ||
            this -> mnM != (bTransp ? mB.nsize() : mB.msize()) ||
            bTransp ? mA.msize() != mB.msize() : mA.nsize() != mB.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __syr2k<TR, basic_srsmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), mB, mB.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // matrix inversion
    basic_srsmatrix& inv (const basic_srsmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_srsmatrix>(*this, mArg);
        return *this;
    }

    basic_srsmatrix inv() const throw (cvmexception)
    {
        basic_srsmatrix mRes (this -> mnM);
        __inv<basic_srsmatrix>(mRes, *this);
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_srsmatrix& exp (const basic_srsmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp_symm<basic_srsmatrix, TR> (*this, mArg, tol);
        return *this;
    }

    basic_srsmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_srsmatrix msRes (this -> mnM);
        __exp_symm<basic_srsmatrix, TR> (msRes, *this, tol);
        return msRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_srsmatrix& polynom (const basic_srsmatrix& m, const RVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), v.incr() > 1 ? v1 : v);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_srsmatrix polynom (const RVector& v) const
    {
        basic_srsmatrix msRes (this -> mnM);
        RVector v1;
        if (v.incr() > 1) v1 << v;   // to make sure incr = 1
        __polynom<TR, RVector> (msRes, msRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, v.incr() > 1 ? v1 : v);
        return msRes;
    }

    // eigenvalues
    // we don't use _eig here since this is the special case - symmetric matrix
    RVector eig (BaseSRMatrix& mEigVect) const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_srsmatrix, BaseSRMatrix> (vEig, *this, &mEigVect, true);
        return vEig;
    }

    RVector eig() const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_srsmatrix, BaseSRMatrix> (vEig, *this, NULL, true);
        return vEig;
    }

    // Cholesky factorization
    BaseSRMatrix cholesky () const throw (cvmexception)
    {
        BaseSRMatrix mRes (*this);
        int nOutInfo = __cholesky<BaseSRMatrix> (mRes);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        mRes._clean_low_triangle ();
        return mRes;
    }

    // Bunch-Kaufman factorization
    BaseSRMatrix bunch_kaufman (int* nPivots) const throw (cvmexception)
    {
        BaseSRMatrix mRes (*this);
        __bunch_kaufman<BaseSRMatrix> (mRes, nPivots);
        return mRes;
    }

    basic_srsmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_srsmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_srsmatrix& randomize (TR dFrom, TR dTo)
    {
        this -> _randomize (dFrom, dTo);
        return *this;
    }

    // matrix equilibration (useful for further solve and solve_lu calling)
    // returns true if equilibration was needed and performed
    bool equilibrate (RVector& vScalings, RVector& vB) throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vB[i] *= vScalings[i];
            }
        }
        return bRes;
    }

    bool equilibrate (RVector& vScalings, BaseRMatrix& mB) throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int j = 1; j <= mB.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mB(i,j) *= vScalings[i];
                }
            }
        }
        return bRes;
    }

    // special care for symmetric matrices
    basic_srsmatrix& _factorize (const basic_srsmatrix& m, int* nPivots, bool& bPositiveDefinite) throw (cvmexception)
    {
        (*this) = m;
        int nOutInfo = __cholesky<BaseSRMatrix>(*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0)
        {
            (*this) = m;
            __bunch_kaufman<BaseSRMatrix>(*this, nPivots);
            bPositiveDefinite = false;
        }
        else {
            bPositiveDefinite = true;
        }
        return *this;
    }

    virtual TR* _pd()                                           // redefinition of Array's function
    {
        return this -> mpD;
    }

    virtual const TR* _pd() const                               // redefinition of Array's function
    {
        return this -> mpD;
    }

    virtual TR norminf() const                                  // infinity-norm - the same as 1-norm for symmetric matrices
    {
        return this -> norm1();
    }

    // makes lower triangle to be equal to upper one
    virtual void _flip()
    {
        if (this -> mnM > 1)
        {
            const int nM1 = this -> mnLD + 1, nM1m = this -> mnLD - 1, nM2m = this -> mnM - 1;
            int i = 1, j = 1, m;
            for (;;)
            {
                m = this -> mnM - i;
                __copy<TR> (m, this -> mpD + j + nM1m, this -> mnLD, this -> mpD + j, 1);
                if (i >= nM2m)
                {
                    break;
                }
                i++;
                j += nM1;
            }
        }
    }

    virtual void _check_submatrix ()  throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "srsmatrix");}

protected:
    virtual const TR* _p (const BaseMatrix& m) const  // for _msum _mdiff etc.
    {
        return m.get();
    }

    virtual void _check_ger()         throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "ger", "srsmatrix");}
    virtual void _check_rank1update() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update", "srsmatrix");}
    virtual void _check_gemm()        throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "srsmatrix");}
    virtual void _swap_rows(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "srsmatrix");}
    virtual void _swap_cols(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "srsmatrix");}

    // we do nothing here - it's symmetric 
    virtual void _transp()
    {
    }

    // returns main diagonal of low_up factorization
    // this call is useless for symmetric matrices. this call would mean serious CVM internal error
    virtual RVector _low_up_diag (basic_array<int>&) const throw (cvmexception)
    {
        throw cvmexception (CVM_NOTIMPLEMENTED);
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TR> (this -> mpD, this -> mnSize, this -> mnIncr, d);     // zero tails are supposed here
    }

    virtual void _multiply (RVector& vRes, const RVector& v, bool) const
    {
        RVector vTmp;
        basic_srsmatrix mTmp;
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        const TR* pDm = this -> mpD;
        const TR* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == pDm) mTmp << *this;
        __symv<TR, basic_srsmatrix, RVector> (vRes.get() == pDm ? mTmp : *this, one, 
                                              vRes.get() == pDv ? vTmp : v, zero, vRes);
    }

    virtual void _randomize (TR dFrom, TR dTo)
    {
        this -> BaseSRMatrix::_randomize (dFrom, dTo);
        this -> _flip();
    }

    virtual void _solve (const RVector& vB, RVector& vX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        RVector vB1 (vB);    // to make sure incr = 1 and equilibrate
        RVector vScalings (this -> mnM);
        basic_srsmatrix m (*this);
        const bool bEquilibrated = m.equilibrate (vScalings, vB1);
        RVector vX1 (vB1);   // to make sure incr = 1
        __solve<TR, TR, basic_srsmatrix> (m, 1, vB1, vB1.size(), vX1, vX1.size(), dErr, pLU, pPivots);
        if (bEquilibrated)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vX[i] = vX1[i] * vScalings[i];
            }
        }
        else
        {
            vX = vX1;
        }
    }

    virtual void _solve (const BaseRMatrix& mB, BaseRMatrix& mX, TR& dErr, const TR* pLU, const int* pPivots) const throw (cvmexception)
    {
        BaseRMatrix mB1 (mB);    // to equilibrate
        RVector vScalings (this -> mnM);
        basic_srsmatrix m (*this);
        const bool bEquilibrated = m.equilibrate (vScalings, mB1);
        mX = mB1;
        __solve<TR, TR, basic_srsmatrix> (m, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
        if (bEquilibrated)
        {
            for (int j = 1; j <= mX.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mX(i,j) *= vScalings[i];
                }
            }
        }
    }

    // matrix determinant
    virtual TR _det() const throw (cvmexception)
    {
        TR dDet = TR(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                dDet = this -> _ij_val (0, 0);
                break;
            case 2:
                dDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TR one(1.);
                    bool bPositiveDefinite = false;
                    basic_srsmatrix m(this -> mnM);
                    basic_array<int> nPivots (this -> mnM);
                    m._factorize (*this, nPivots, bPositiveDefinite);

                    int i;
                    dDet = one;
                    if (bPositiveDefinite)
                    {
                        const RVector vUpDiag = m.diag(0);
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vUpDiag[i] * vUpDiag[i];    //here we use Cholesky factorization
                        }
                    }
                    else
                    {
                        const RVector vEig = this -> eig();
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vEig[i];                    //here we use eigenvalues. probably there is a better way.
                        }
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_WRONGBUNCHKAUFMANFACTOR) throw e;
                }
                break;
        }
        return dDet;
    }

    // matrix equilibration helper
    bool _equilibrate (RVector& vScalings) throw (cvmexception)
    {
        if (this -> mnM != vScalings.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = false;
        TR dCond(0.);
        TR dMax  = TR();
        static const TR sp = basic_cvmMachSp<TR>();
        static const TR sp_inv = TR(1.) / sp;

        __poequ<TR, basic_srsmatrix, RVector> (*this, vScalings, dCond, dMax);

        if (dCond < TR(0.1) || _abs(dMax) <= sp || _abs(dMax) >= sp_inv)
        {
            bRes = true;
            for (int i = 0; i < this -> mnM; ++i)
            {
                for (int j = i; j < this -> mnM; ++j)
                {
                    this -> mpD[this -> mnLD * j + i] *= vScalings[i+1] * vScalings[j+1];
                }
            }
        }
        return bRes;
    }

    // sets both elements to keep matrix symmetric, checks ranges
    void _set_at (int nRow, int nCol, TR val) throw (cvmexception)      // zero based
    {
        if (nRow < 0 || nRow >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nRow);
        if (nCol < 0 || nCol >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nCol);
        this -> mpD[this -> mnLD * nCol + nRow] = val;
        if (nRow != nCol)
        {
            this -> mpD[this -> mnLD * nRow + nCol] = val;
        }
    }

    void _check_symmetric () const throw (cvmexception)
    {
        for (int j = 0; j < this -> mnN; ++j)
        {
            for (int i = 0; i < this -> mnM; ++i)
            {
                if (i == j) continue;
                if (_abs (this -> mpD[this -> mnLD * j + i] - this -> mpD[this -> mnLD * i + j]) > basic_cvmMachMin<TR>()) throw cvmexception (CVM_MATRIXNOTSYMMETRIC);
            }
        }
    }
};

template <typename TR, typename TC>
class basic_schmatrix : public basic_scmatrix<TR,TC>
{
    typedef basic_rvector<TR> RVector;
    typedef basic_cvector<TR,TC> CVector;
    typedef Matrix<TR,TC> BaseMatrix;
    typedef basic_cmatrix<TR,TC> BaseCMatrix;
    typedef basic_scmatrix<TR,TC> BaseSCMatrix;

public:
    basic_schmatrix()
    {
    }

    explicit basic_schmatrix (int nMN)
        : BaseSCMatrix (nMN)
    {
    }

    basic_schmatrix (TC* pD, int nMN)
        : BaseSCMatrix (pD, nMN, nMN, nMN * nMN)
    {
        this -> _check_hermitian();
    }

    basic_schmatrix (const basic_schmatrix& m)
        : BaseSCMatrix (m.mnM, m.mnM, false)
    {
        this -> _massign(m);
    }

    basic_schmatrix (const BaseCMatrix& m)
        : BaseSCMatrix (m.msize(), m.nsize(), false)
    {
        if (this -> mnM != this -> mnN) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        this -> _check_hermitian();
    }

    // diagonal square symmetric matrix constructor
    explicit basic_schmatrix (const RVector& v)
        : BaseSCMatrix (v.size(), v.size(), true)
    {
        __copy2<TR,TC> (this -> mpD, v.size(), this -> mnLD + 1, v._pd(), NULL, v.incr());
    }

    explicit basic_schmatrix (const basic_srsmatrix<TR>& m)
        : BaseSCMatrix (m.msize(), m.msize(), true)
    {
        __copy2<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, m._pd(), NULL);
    }

    basic_schmatrix (const TR* pRe, const TR* pIm, int nMN)
        : BaseSCMatrix (pRe, pIm, nMN)
    {
        this -> _check_hermitian();
    }

    basic_schmatrix (const basic_srmatrix<TR>& mRe, const basic_srmatrix<TR>& mIm)
        : BaseSCMatrix (mRe, mIm)
    {
        this -> _check_hermitian();
    }

    // submatrix constructor
    basic_schmatrix (basic_schmatrix& m, int nRowCol, int nSize)        // 1-based
        : BaseSCMatrix (m._sub_pointer (nRowCol, nRowCol, nSize, nSize), nSize, m.ld(), nSize * nSize)
    {
    }

    TC operator () (int nIm, int nIn) const throw (cvmexception)
    {
        if (nIm <= 0 || nIm > this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nIm);
        if (nIn <= 0 || nIn > this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nIn);
        return this -> _ij_val (nIm - 1, nIn - 1);
    }

    // returns column which CAN NOT be l-value
    const CVector operator () (int nFI) const throw (cvmexception)
    {
        return this -> _col (nFI - 1);
    }

    // returns row which CAN NOT be l-value
    const CVector operator [] (int nFI) const throw (cvmexception)
    {
        return this -> _row (nFI - 1);
    }

    // returns diagonal which IS NOT l-value
    // 0 - main, negative - low, positive - up
    const CVector diag (int nDiag) const throw (cvmexception)
    {
        CVector vRet (this -> _diag(nDiag));
        return vRet;
    }

    // real part (symmetric)
    const basic_srsmatrix<TR> real() const
    {
        basic_srsmatrix<TR> mRet (this -> mnM);
        __copy<TR> (this -> mnSize, __get_real_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    // imaginary part (NOT symmetric)
    const basic_srmatrix<TR> imag() const
    {
        basic_srmatrix<TR> mRet (this -> mnM);
        __copy<TR> (this -> mnSize, __get_imag_p<TR>(this -> mpD), this -> mnIncr * 2, mRet, mRet.incr());
        return mRet;
    }

    basic_schmatrix& operator = (const basic_schmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _massign(m);
        return *this;
    }

    // assigns foregn array (nIncr = 1)
    basic_schmatrix& assign (const TC* pD) throw (cvmexception)
    {
        this -> _assign (pD, 1);
        this -> _check_hermitian();
        return *this;
    }

    basic_schmatrix& assign (int nRowCol, const basic_schmatrix& m) throw (cvmexception)         // submatrix assignment
    {
        if (nRowCol <= 0 || nRowCol > this -> mnM) throw cvmexception (CVM_OUTOFRANGE, nRowCol);
        if (m.mnM + nRowCol - 1 > this -> mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _assign_shifted (this -> _sub_pointer_nocheck (nRowCol, nRowCol), m._pd(), m.mnM, m.mnN, m.mnLD);
        return *this;
    }

    // sets both elements to keep matrix hermitian
    basic_schmatrix& set (int nRow, int nCol, TC d) throw (cvmexception)
    {
        this -> _set_at (nRow - 1, nCol - 1, d);
        return *this;
    }

    // sets both diagonals (keeps matrix to be hermitian)
    basic_schmatrix& set_diag (int nDiag, const CVector& vDiag) throw (cvmexception)
    {
        if (nDiag == 0) throw cvmexception (CVM_BREAKS_HERMITIANITY, "set_main_diag");
//        CVector (this -> _diag(nDiag)) = vDiag;               Borland's "optimizer" fails here
//        (CVector (this -> _diag(-nDiag)) = vDiag).conj();
        const int nD = abs (nDiag);
        const int nSize = vDiag.size();
        if (this -> mnM - nD != nSize) throw cvmexception (CVM_SIZESMISMATCH);
        const int nShift = nDiag * this -> mnLD;
        const int nIncr = this -> mnLD + 1;
        TC* pD1 = this -> mpD + (nDiag > 0 ? nShift : nD);
        TC* pD2 = this -> mpD + (nDiag > 0 ? nD : nShift);
        __copy<TC> (nSize, vDiag, vDiag.incr(), pD1, nIncr);
        __copy<TC> (nSize, vDiag, vDiag.incr(), pD2, nIncr);
        __conj<TC> (pD2, nSize, nIncr);
        return *this;
    }

    // sets main diagonal (keeps matrix to be hermitian)
    basic_schmatrix& set_main_diag (const RVector& vDiag) throw (cvmexception)
    {
//        CVector vRet (this -> _diag(0));                      Borland's "optimizer" fails here
//        vRet.real() = vDiag;
        if (this -> mnM != vDiag.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR,TC> (this -> mpD, this -> mnM, this -> mnLD + 1, vDiag, vDiag.incr());
        return *this;
    }

    // assigns real matrix to real part
    basic_schmatrix& assign_real (const basic_srsmatrix<TR>& mRe) throw (cvmexception)
    {
        if (this -> mnM != mRe.msize() || this -> mnN != mRe.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __copy_real<TR,TC> (this -> mpD, this -> mnSize, this -> mnIncr, mRe._pd(), mRe.incr());
        return *this;
    }

    // fills real part
    basic_schmatrix& set_real (TR d)
    {
        this -> _set_real_number (d);
        return *this;
    }

    basic_schmatrix& resize (int nNewMN) throw (cvmexception)
    {
        this -> _resize (nNewMN, nNewMN);
        return *this;
    }

    bool operator == (const basic_schmatrix& m) const
    {
        return this -> mnM == m.mnM && this -> mnN == m.mnN && this -> _mequals (m);
    }

    bool operator != (const basic_schmatrix& m) const
    {
        return !(this -> operator == (m));
    }

    basic_schmatrix& operator << (const basic_schmatrix& m) throw (cvmexception)
    {
        this -> _replace (m);
        this -> _massign (m);
        return *this;
    }

    basic_schmatrix operator + (const basic_schmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_schmatrix mSum (*this);
        mSum._mincr (m);
        return mSum;
    }

    basic_schmatrix operator - (const basic_schmatrix& m) const throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        basic_schmatrix mDiff (*this);
        mDiff._mdecr (m);
        return mDiff;
    }

    basic_schmatrix& sum (const basic_schmatrix& m1, const basic_schmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _msum (m1, m2);
        return *this;
    }

    basic_schmatrix& diff (const basic_schmatrix& m1, const basic_schmatrix& m2) throw (cvmexception)
    {
        if (this -> mnM != m1.mnM || this -> mnM != m2.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdiff (m1, m2);
        return *this;
    }

    basic_schmatrix& operator += (const basic_schmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mincr (m);
        return *this;
    }

    basic_schmatrix& operator -= (const basic_schmatrix& m) throw (cvmexception)
    {
        if (this -> mnM != m.mnM) throw cvmexception (CVM_SIZESMISMATCH);
        this -> _mdecr (m);
        return *this;
    }

    basic_schmatrix operator - () const
    {
        static const TR mone(-1.);
        basic_schmatrix mRes (*this);
        mRes._scal (mone);
        return mRes;
    }

    // plus identity, prefix
    basic_schmatrix& operator ++ ()
    {
        this -> _plus_plus();
        return *this;
    }

    // plus identity, postfix
    basic_schmatrix& operator ++ (int)
    {
        this -> _plus_plus();
        return *this;
    }

    // minus identity, prefix
    basic_schmatrix& operator -- ()
    {
        this -> _minus_minus();
        return *this;
    }

    // minus identity, postfix
    basic_schmatrix& operator -- (int)
    {
        this -> _minus_minus();
        return *this;
    }

    basic_schmatrix operator * (TR dMult) const
    {
        basic_schmatrix mRes (*this);
        mRes._scal (dMult);
        return mRes;
    }

    basic_schmatrix operator / (TR dDiv) const throw (cvmexception)
    {
        basic_schmatrix mRes (*this);
        mRes._div (dDiv);
        return mRes;
    }

    BaseSCMatrix operator * (TC cMult) const
    {
        BaseSCMatrix mRes (*this);
        mRes._scal (cMult);
        return mRes;
    }

    BaseSCMatrix operator / (TC cDiv) const throw (cvmexception)
    {
        BaseSCMatrix mRes (*this);
        mRes._div (cDiv);
        return mRes;
    }

    basic_schmatrix& operator *= (TR dMult)
    {
        this -> _scal (dMult);
        return *this;
    }

    basic_schmatrix& operator /= (TR dDiv) throw (cvmexception)
    {
        this -> _div (dDiv);
        return *this;
    }

    basic_schmatrix& normalize()                                           // array normalizing
    {
        this -> _normalize();
        return *this;
    }

    // hermitian conjugated Matrix
    basic_schmatrix operator ~ () const
    {
        return *this;
    }

    // to override scmatrix
    // this = hermitian conjugate (m)
    basic_schmatrix& conj (const basic_schmatrix& m) throw (cvmexception)
    {
        (*this) = m;
        return *this;
    }

    // this = hermitian conjugate (this)
    basic_schmatrix& conj()
    {
        return *this;
    }

    CVector operator * (const CVector& v) const throw (cvmexception)
    {
        if (this -> mnN != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vRes (this -> mnM);
        this -> _multiply (vRes, v, false);
        return vRes;
    }

    // special exclusion since matrix product is not commutative
    BaseCMatrix operator * (const BaseCMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseCMatrix mRes(m.msize(), m.nsize());
        mRes._hemm (true, *this, m, one, zero);
        return mRes;
    }

    BaseSCMatrix operator * (const BaseSCMatrix& m) const throw (cvmexception)
    {
        if (this -> mnN != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        BaseSCMatrix mRes(this -> mnM);
        mRes._hemm (true, *this, m, one, zero);
        return mRes;
    }

    // this = alpha*v*v' + beta*this
    basic_schmatrix& herk (TC alpha, const CVector& v, TC beta) throw (cvmexception)
    {
        if (this -> mnM != v.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __herk<TC, basic_schmatrix> (false, alpha, 1, v, v.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*a' + beta*this or this = alpha*a'*a + beta*this
    basic_schmatrix& herk (bool bTransp, TC alpha, const BaseCMatrix& mA, TC beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize())) throw cvmexception (CVM_SIZESMISMATCH);
        __herk<TC, basic_schmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*v1*v2' + alpha*v2*v1' + beta*this
    basic_schmatrix& her2k (TC alpha, const CVector& v1, const CVector& v2, TC beta) throw (cvmexception)
    {
        if (this -> mnM != v1.size() || this -> mnM != v2.size()) throw cvmexception (CVM_SIZESMISMATCH);
        __her2k<TC, basic_schmatrix> (false, alpha, 1, v1, v1.size(), v2, v2.size(), beta, *this);
        this -> _flip();
        return *this;
    }

    // this = alpha*a*b' + alpha*b*a' + beta*this or this = alpha*a'*b + alpha*b'*a + beta*this
    basic_schmatrix& her2k (bool bTransp, TC alpha, const BaseCMatrix& mA, const BaseCMatrix& mB, TC beta) throw (cvmexception)
    {
        if (this -> mnM != (bTransp ? mA.nsize() : mA.msize()) ||
            this -> mnM != (bTransp ? mB.nsize() : mB.msize()) ||
            bTransp ? mA.msize() != mB.msize() : mA.nsize() != mB.nsize()) throw cvmexception (CVM_SIZESMISMATCH);
        __her2k<TC, basic_schmatrix> (bTransp, alpha, bTransp ? mA.msize() : mA.nsize(), mA, mA.ld(), mB, mB.ld(), beta, *this);
        this -> _flip();
        return *this;
    }

    // matrix inversion
    basic_schmatrix& inv (const basic_schmatrix& mArg) throw (cvmexception)
    {
        __inv<basic_schmatrix> (*this, mArg);
        return *this;
    }

    basic_schmatrix inv() const throw (cvmexception)
    {
        basic_schmatrix mRes (this -> mnM);
        __inv<basic_schmatrix> (mRes, *this);
        return mRes;
    }

    // matrix exponent with given tolerance
    basic_schmatrix& exp (const basic_schmatrix& mArg, TR tol = basic_cvmMachSp<TR>()) throw (cvmexception)
    {
        __exp_symm<basic_schmatrix, TR> (*this, mArg, tol);
        return *this;
    }

    basic_schmatrix exp (TR tol = basic_cvmMachSp<TR>()) const throw (cvmexception)
    {
        basic_schmatrix msRes (this -> mnM);
        __exp_symm<basic_schmatrix, TR> (msRes, *this, tol);
        return msRes;
    }

    // this = v(1)*I + v(2)*m + v(3)*m^2 + ... + v(N)*m^(N-1)
    basic_schmatrix& polynom (const basic_schmatrix& m, const RVector& v) throw (cvmexception)
    {
        if (this -> mnM != m.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        CVector vc(v);
        __polynom<TC, CVector> (this -> mpD, this -> mnLD, this -> mnM, m._pd(), m._ldm(), vc);
        return *this;
    }

    // returns v(1)*I + v(2)*this + v(3)*this^2 + ... + v(N)*this^(N-1)
    basic_schmatrix polynom (const RVector& v) const
    {
        basic_schmatrix msRes (this -> mnM);
        CVector vc(v);
        __polynom<TC, CVector> (msRes, msRes.mnLD, this -> mnM, this -> mpD, this -> mnLD, vc);
        return msRes;
    }

    // eigenvalues
    // we don't use _eig here since this is the special case - hermitian matrix
    RVector eig (BaseSCMatrix& mEigVect) const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_schmatrix, BaseSCMatrix> (vEig, *this, &mEigVect, true);
        return vEig;
    }

    RVector eig() const throw (cvmexception)
    {
        RVector vEig(this -> mnM);
        __eig<RVector, basic_schmatrix, BaseSCMatrix> (vEig, *this, NULL, true);
        return vEig;
    }

    // Cholesky factorization
    BaseSCMatrix cholesky () const throw (cvmexception)
    {
        BaseSCMatrix mRes (*this);
        int nOutInfo = __cholesky<BaseSCMatrix> (mRes);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0) throw cvmexception (CVM_NOTPOSITIVEDEFINITE, nOutInfo);
        mRes._clean_low_triangle ();
        return mRes;
    }

    // Bunch-Kaufman factorization
    BaseSCMatrix bunch_kaufman (int* nPivots) const throw (cvmexception)
    {
        BaseSCMatrix mRes (*this);
        __bunch_kaufman<BaseSCMatrix> (mRes, nPivots);
        return mRes;
    }

    basic_schmatrix& identity()
    {
        this -> _vanish();
        this -> _plus_plus();
        return *this;
    }

    basic_schmatrix& vanish()
    {
        this -> _vanish();
        return *this;
    }

    basic_schmatrix& randomize_real (TR dFrom, TR dTo)
    {
        this -> _randomize_real (dFrom, dTo);
        return *this;
    }

    basic_schmatrix& randomize_imag (TR dFrom, TR dTo)
    {
        this -> _randomize_imag (dFrom, dTo);
        return *this;
    }

    // matrix equilibration (useful for further solve and solve_lu calling)
    // returns true if equilibration was needed and performed
    bool equilibrate (basic_rvector<TR>& vScalings, CVector& vB) throw (cvmexception)
    {
        if (this -> mnM != vB.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vB[i] *= vScalings[i];
            }
        }
        return bRes;
    }

    bool equilibrate (basic_rvector<TR>& vScalings, BaseCMatrix& mB) throw (cvmexception)
    {
        if (this -> mnM != mB.msize()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = this -> _equilibrate (vScalings);
        if (bRes)
        {
            for (int j = 1; j <= mB.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mB(i,j) *= vScalings[i];
                }
            }
        }
        return bRes;
    }

    // special care for symmetric matrices
    basic_schmatrix& _factorize (const basic_schmatrix& m, int* nPivots, bool& bPositiveDefinite) throw (cvmexception)
    {
        (*this) = m;
        int nOutInfo = __cholesky<BaseSCMatrix>(*this);
        if (nOutInfo < 0) throw cvmexception (CVM_WRONGMKLARG);
        if (nOutInfo > 0)
        {
            (*this) = m;
            __bunch_kaufman<BaseSCMatrix>(*this, nPivots);
            bPositiveDefinite = false;
        }
        else {
            bPositiveDefinite = true;
        }
        return *this;
    }

    // redefinition of Array's function
    virtual TC* _pd()
    {
        return this -> mpD;
    }

    // redefinition of Array's function
    virtual const TC* _pd() const
    {
        return this -> mpD;
    }

    // infinity-norm - the same as 1-norm for symmetric matrices
    virtual TR norminf() const
    {
        return this -> norm1();
    }

    // makes lower triangle to be equal to conjugated upper one
    virtual void _flip()
    {
        if (this -> mnM > 1)
        {
            const int nM1 = this -> mnLD + 1, nM1m = this -> mnLD - 1, nM2m = this -> mnM - 1;
            int i = 1, j = 1, m;
            for (;;)
            {
                m = this -> mnM - i;
                __copy<TC> (m, this -> mpD + j + nM1m, this -> mnLD, this -> mpD + j, 1);
                __conj<TC> (this -> mpD + j, m, 1);
                if (i >= nM2m) break;
                ++i;
                j += nM1;
            }
        }
    }

    virtual void _check_submatrix ()  throw (cvmexception) {throw cvmexception (CVM_SUBMATRIXNOTAVAILABLE, "schmatrix");}

protected:
    // not applicable - see below
    basic_schmatrix& set (TC z) throw (cvmexception)
    {
        this -> _set (z);
        return *this;
    }

    virtual void _set (TC) throw (cvmexception)
    {
        throw cvmexception (CVM_BREAKS_HERMITIANITY, "set_real");
    }

    virtual const TC* _p (const BaseMatrix& m) const  // for _msum _mdiff etc.
    {
        return m.get();
    }

    virtual void _randomize_real (TR dFrom, TR dTo)
    {
        this -> BaseSCMatrix::_randomize_real (dFrom, dTo);
        this -> _flip();
        this -> _make_main_diag_real ();
    }

    virtual void _randomize_imag (TR dFrom, TR dTo)
    {
        this -> BaseSCMatrix::_randomize_imag (dFrom, dTo);
        this -> _flip();
        this -> _make_main_diag_real ();
    }

    virtual void _transp()                                      // __conj execution assumed after that
    {
        this -> _sq_transp();
    }

    // returns main diagonal of low_up factorization
    virtual CVector _low_up_diag (basic_array<int>&) const throw (cvmexception)
    {
        throw cvmexception (CVM_NOTIMPLEMENTED);    // well, this stuff is useless for symmetric matrices. this call would mean serious CVM internal error
    }

    virtual void _scal (TR d)
    {
        __scal<TR, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    virtual void _scal (TC d)
    {
        __scal<TC, TC> (this -> mpD, this -> mnSize, this -> mnIncr, d);
    }

    virtual void _multiply (CVector& vRes, const CVector& v, bool) const
    {
        CVector vTmp;
        basic_schmatrix mTmp;
        static const TR zero = TR(0.);
        static const TR one = TR(1.);
        const TC* pDm = this -> mpD;
        const TC* pDv = v;
        if (vRes.get() == pDv) vTmp << v;
        if (vRes.get() == pDm) mTmp << *this;
        __shmv<TC, basic_schmatrix, CVector> (vRes.get() == pDm ? mTmp : *this, one, 
                                              vRes.get() == pDv ? vTmp : v, zero, vRes);
    }

    // sets both elements to keep matrix hermitian, checks ranges
    void _set_at (int nRow, int nCol, TC val) throw (cvmexception)      // zero based
    {
        if (nRow < 0 || nRow >= this -> mnM) throw cvmexception (CVM_OUTOFRANGE1, nRow);
        if (nCol < 0 || nCol >= this -> mnN) throw cvmexception (CVM_OUTOFRANGE2, nCol);
        if (nRow == nCol && _abs (val.imag()) > basic_cvmMachMin<TR>())     // only reals on main diagonal
        {
            throw cvmexception (CVM_BREAKS_HERMITIANITY, "real number");
        }
        this -> mpD[this -> mnLD * nCol + nRow] = val;
        if (nRow != nCol)
        {
            this -> mpD[this -> mnLD * nRow + nCol] = _conjugate(val);
        }
    }

    void _check_hermitian () const throw (cvmexception)
    {
        for (int j = 0; j < this -> mnN; ++j)
        {
            for (int i = 0; i < this -> mnM; ++i)
            {
                if (i == j) {
                    if (_abs (this -> mpD[this -> mnLD * j + i].imag()) > basic_cvmMachMin<TR>())    // real numbers on main diagonal
                    {
                        throw cvmexception (CVM_MATRIXNOTHERMITIAN);
                    }
                    continue;
                }
                if (!_conjugated<TR> (this -> mpD[this -> mnLD * j + i], this -> mpD[this -> mnLD * i + j])) throw cvmexception (CVM_MATRIXNOTHERMITIAN);
            }
        }
    }

    virtual void _solve (const CVector& vB, CVector& vX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        CVector vB1 (vB);    // to make sure incr = 1 and equilibrate
        basic_rvector<TR> vScalings (this -> mnM);
        basic_schmatrix m(*this);
        const bool bEquilibrated = m.equilibrate (vScalings, vB1);
        CVector vX1 (vB1);   // to make sure incr = 1
        __solve<TR, TC, basic_schmatrix> (m, 1, vB1, vB1.size(), vX1, vX1.size(), dErr, pLU, pPivots);

        if (bEquilibrated)
        {
            for (int i = 1; i <= this -> mnM; ++i)
            {
                vX[i] = vX1[i] * vScalings[i];
            }
        }
        else
        {
            vX = vX1;
        }
    }

    virtual void _solve (const BaseCMatrix& mB, BaseCMatrix& mX, TR& dErr, const TC* pLU, const int* pPivots) const throw (cvmexception)
    {
        BaseCMatrix mB1 (mB);    // to equilibrate
        basic_rvector<TR> vScalings (this -> mnM);
        basic_schmatrix m(*this);
        const bool bEquilibrated = m.equilibrate (vScalings, mB1);
        mX = mB1;
        __solve<TR, TC, basic_schmatrix> (m, mB.nsize(), mB, mB.ld(), mX, mX.ld(), dErr, pLU, pPivots);
        if (bEquilibrated)
        {
            for (int j = 1; j <= mX.nsize(); ++j)
            {
                for (int i = 1; i <= this -> mnM; ++i)
                {
                    mX(i,j) *= vScalings[i];
                }
            }
        }
    }

    // matrix determinant
    virtual TC _det() const throw (cvmexception)
    {
        TC dDet = TC(0.);
        switch (this -> mnM)
        {
            case 0:
                break;
            case 1:
                dDet = this -> _ij_val (0, 0);
                break;
            case 2:
                dDet = this -> _ij_val (0, 0) * this -> _ij_val (1, 1) - 
                       this -> _ij_val (1, 0) * this -> _ij_val (0, 1);
                break;
            default:
                try
                {
                    static const TC one(1.);
                    bool bPositiveDefinite = false;
                    basic_schmatrix m(this -> mnM);
                    basic_array<int> nPivots (this -> mnM);
                    m._factorize (*this, nPivots, bPositiveDefinite);

                    int i;
                    dDet = one;
                    if (bPositiveDefinite)
                    {
                        const CVector vUpDiag = m.diag(0);
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vUpDiag[i] * vUpDiag[i];    //here we use Cholesky factorization
                        }
                    }
                    else
                    {
                        const basic_rvector<TR> vEig = this -> eig();
                        for (i = 1; i <= this -> mnM; ++i)
                        {
                            dDet *= vEig[i];                    //here we use eigenvalues. probably there is a better way.
                        }
                    }
                }
                catch (const cvmexception& e)
                {
                    if (e.cause() != CVM_WRONGBUNCHKAUFMANFACTOR) throw e;
                }
                break;
        }
        return dDet;
    }

    // matrix equilibration helper
    bool _equilibrate (basic_rvector<TR>& vScalings) throw (cvmexception)
    {
        if (this -> mnM != vScalings.size()) throw cvmexception (CVM_SIZESMISMATCH);
        bool bRes = false;
        TR dCond = TR();
        TR dMax  = TR();
        static const TR sp = basic_cvmMachSp<TR>();
        static const TR sp_inv = TR(1.) / sp;

        __poequ<TR, basic_schmatrix, basic_rvector<TR> > (*this, vScalings, dCond, dMax);

        if (dCond < TR(0.1) || _abs(dMax) <= sp || _abs(dMax) >= sp_inv)
        {
            bRes = true;
            for (int i = 0; i < this -> mnM; ++i)
            {
                for (int j = i; j < this -> mnM; ++j)
                {
                    this -> mpD[this -> mnLD * j + i] *= vScalings[i+1] * vScalings[j+1];
                }
            }
        }
        return bRes;
    }

    void _make_main_diag_real ()
    {
        static const TR zero = TR(0.);
        __scal<TR,TR> (__get_imag_p<TR> (this -> mpD), this -> mnM, (this -> mnLD + 1) * 2, zero);
    }

    virtual void _check_ger()         throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "ger", "schmatrix");}
    virtual void _check_rank1update() throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "rank1update", "schmatrix");}
    virtual void _check_gemm()        throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "gemm", "schmatrix");}
    virtual void _swap_rows(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_rows", "schmatrix");}
    virtual void _swap_cols(int, int) throw (cvmexception) {throw cvmexception (CVM_METHODNOTAVAILABLE, "swap_cols", "schmatrix");}
};


template <typename TR>
inline basic_rvector<TR> operator * (TR d, const basic_rvector<TR>& v)
{
    return v * d;
}
template <typename TR>
inline basic_rmatrix<TR> operator * (TR d, const basic_rmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srmatrix<TR> operator * (TR d, const basic_srmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srbmatrix<TR> operator * (TR d, const basic_srbmatrix<TR>& m)
{
    return m * d;
}
template <typename TR>
inline basic_srsmatrix<TR> operator * (TR d, const basic_srsmatrix<TR>& m)
{
    return m * d;
}

template <typename TR, typename TC>
inline basic_cvector<TR,TC> operator * (TR d, const basic_cvector<TR,TC>& v)
{
    return v * d;
}
template <typename TR, typename TC>
inline basic_cmatrix<TR,TC> operator * (TR d, const basic_cmatrix<TR,TC>& m)
{
    return m * d;
}
template <typename TR, typename TC>
inline basic_scmatrix<TR,TC> operator * (TR d, const basic_scmatrix<TR,TC>& m)
{
    return m * d;
}
template <typename TR, typename TC>
inline basic_scbmatrix<TR,TC> operator * (TR d, const basic_scbmatrix<TR,TC>& m)
{
    return m * d;
}
template <typename TR, typename TC>
inline basic_schmatrix<TR,TC> operator * (TR d, const basic_schmatrix<TR,TC>& m)
{
    return m * d;
}

template <typename TR, typename TC>
inline basic_cvector<TR,TC> operator * (std::complex<TR> c, const basic_cvector<TR,TC>& v)
{
    return v * c;
}
template <typename TR, typename TC>
inline basic_cmatrix<TR,TC> operator * (std::complex<TR> c, const basic_cmatrix<TR,TC>& m)
{
    return m * c;
}
template <typename TR, typename TC>
inline basic_scmatrix<TR,TC> operator * (std::complex<TR> c, const basic_scmatrix<TR,TC>& m)
{
    return m * c;
}
template <typename TR, typename TC>
inline basic_scbmatrix<TR,TC> operator * (std::complex<TR> c, const basic_scbmatrix<TR,TC>& m)
{
    return m * c;
}
template <typename TR, typename TC>
inline basic_schmatrix<TR,TC> operator * (std::complex<TR> c, const basic_schmatrix<TR,TC>& m)
{
    return m * c;
}

template <typename TR>
inline basic_rvector<TR>  operator * (CVM_LONGEST_INT d, const basic_rvector<TR>& v)
{
    return v * static_cast<TR>(d);
}
template <typename TR>
inline basic_rmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_rmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srbmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srbmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR>
inline basic_srsmatrix<TR>  operator * (CVM_LONGEST_INT d, const basic_srsmatrix<TR>& m)
{
    return m * static_cast<TR>(d);
}

template <typename TR, typename TC>
inline basic_cvector<TR,TC>  operator * (CVM_LONGEST_INT d, const basic_cvector<TR,TC>& v)
{
    return v * static_cast<TR>(d);
}
template <typename TR, typename TC>
inline basic_cmatrix<TR,TC>  operator * (CVM_LONGEST_INT d, const basic_cmatrix<TR,TC>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR, typename TC>
inline basic_scmatrix<TR,TC>  operator * (CVM_LONGEST_INT d, const basic_scmatrix<TR,TC>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR, typename TC>
inline basic_scbmatrix<TR,TC>  operator * (CVM_LONGEST_INT d, const basic_scbmatrix<TR,TC>& m)
{
    return m * static_cast<TR>(d);
}
template <typename TR, typename TC>
inline basic_schmatrix<TR,TC>  operator * (CVM_LONGEST_INT d, const basic_schmatrix<TR,TC>& m)
{
    return m * static_cast<TR>(d);
}



#if defined (CVM_FLOAT)
typedef float  treal;
#else
typedef double treal;
#endif

typedef std::complex<treal> tcomplex;

typedef basic_array    <int>             iarray;
typedef basic_rvector  <treal>           rvector;
typedef basic_rmatrix  <treal>           rmatrix;
typedef basic_srmatrix <treal>           srmatrix;
typedef basic_cvector  <treal, tcomplex> cvector;
typedef basic_cmatrix  <treal, tcomplex> cmatrix;
typedef basic_scmatrix <treal, tcomplex> scmatrix;
typedef basic_srbmatrix<treal>           srbmatrix;
typedef basic_scbmatrix<treal, tcomplex> scbmatrix;
typedef basic_srsmatrix<treal>           srsmatrix;
typedef basic_schmatrix<treal, tcomplex> schmatrix;

// identity matrices creation
template <typename TR>
inline const basic_srmatrix<TR> basic_eye_real (int nM)
{
    basic_srmatrix<TR> mI(nM);
    ++mI;
    return mI;
}
template <typename TR, typename TC>
inline const basic_scmatrix<TR,TC> basic_eye_complex (int nM)
{
    basic_scmatrix<TR,TC> mI(nM);
    ++mI;
    return mI;
}

inline const srmatrix eye_real (int nM)
{
    return basic_eye_real<treal>(nM);
}
inline const scmatrix eye_complex (int nM)
{
    return basic_eye_complex<treal, tcomplex>(nM);
}

inline treal cvmMachMin()
{
    return basic_cvmMachMin<treal>();
}
inline treal cvmMachSp()
{
    return basic_cvmMachSp<treal>();
}

CVM_NAMESPACE_END


// BLAS callback error handler
#if !defined (_MSC_VER)
#define XERBLA xerbla_
#endif

extern "C" {
    void   __stdcall XERBLA    (const char* szSubName,
#ifdef CVM_PASS_STRING_LENGTH_TO_FTN_SUBROUTINES
                                const unsigned int nLen,
#endif
                                const int* pnParam) throw (cvm::cvmexception);
}

#if defined (_MSC_VER)
#   pragma warning(pop)
#endif

#endif                  // _CVM_H

