// KL1p - A portable C++ compressed sensing library.
// Copyright (c) 2011-2012 Ren? Gebel
// 
// This file is part of the KL1p C++ library.
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY of fitness for any purpose. 
//
// This library is free software; You can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License (LGPL) 
// as published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// See http://www.opensource.org/licenses for more info.

#ifndef KL1P_WAVELET2DOPERATOR_H
#define KL1P_WAVELET2DOPERATOR_H

#include "Operator.h"




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter=klab::TDelegateWaveletFilter<T> >
class TWavelet2DOperator : public TOperator<T, T>
{
public:

    TWavelet2DOperator(klab::UInt32 height, klab::UInt32 width);
    TWavelet2DOperator(klab::UInt32 height, klab::UInt32 width, const TFilter& filter);
    TWavelet2DOperator(const TWavelet2DOperator<T, TFilter>& op);
    virtual ~TWavelet2DOperator();    

    klab::UInt32     width() const;
    klab::UInt32     height() const;
    const TFilter&          filter() const;
    TFilter&                filter();

    virtual void            apply(const arma::Col<T>& in, arma::Col<T>& out);
    virtual void            applyAdjoint(const arma::Col<T>& in, arma::Col<T>& out);


private:

    TWavelet2DOperator();
    TWavelet2DOperator<T, TFilter>& operator=(const TWavelet2DOperator<T, TFilter>& op);


private:

    klab::UInt32                _height;
    klab::UInt32                _width;

    klab::TFWT2D<T, TFilter>    _fwt;
};

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
class TWavelet2DOperatorSpecialisation
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, TFilter>& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, TFilter>& spec);
    TWavelet2DOperatorSpecialisation<T, TFilter>&   operator=(const TWavelet2DOperatorSpecialisation<T, TFilter>& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::THaarWaveletFilter<T> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::THaarWaveletFilter<T> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::THaarWaveletFilter<T> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::THaarWaveletFilter<T> >&  operator=(const TWavelet2DOperatorSpecialisation<T, klab::THaarWaveletFilter<T> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 2> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 2> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 2> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 2> >&     operator=(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 2> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 4> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 4> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 4> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 4> >&     operator=(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 4> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 6> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 6> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 6> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 6> >&     operator=(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 6> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 8> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 8> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 8> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 8> >&     operator=(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 8> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T>
class TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 10> >
{
public:

    static void     ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 10> >& fwt); 


private:

    TWavelet2DOperatorSpecialisation();
    TWavelet2DOperatorSpecialisation(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 10> >& spec);
    TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 10> >&    operator=(const TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 10> >& spec);
};

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline TWavelet2DOperator<T, TFilter>::TWavelet2DOperator(klab::UInt32 height, klab::UInt32 width) : TOperator<T, T>(width*height),
_height(height), _width(width), 
_fwt()
{}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline TWavelet2DOperator<T, TFilter>::TWavelet2DOperator(klab::UInt32 height, klab::UInt32 width, const TFilter& filter) : TOperator<T, T>(width*height),
_height(height), _width(width), 
_fwt(filter)
{}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline TWavelet2DOperator<T, TFilter>::TWavelet2DOperator(const TWavelet2DOperator<T, TFilter>& op) : TOperator<T, T>(op),
_height(op._height), _width(op._width),
_fwt(op._fwt)
{}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline TWavelet2DOperator<T, TFilter>::~TWavelet2DOperator()
{}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline klab::UInt32 TWavelet2DOperator<T, TFilter>::width() const
{
    return _width;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline klab::UInt32 TWavelet2DOperator<T, TFilter>::height() const
{
    return _height;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline const TFilter& TWavelet2DOperator<T, TFilter>::filter() const
{
    return _fwt.filter();
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline TFilter& TWavelet2DOperator<T, TFilter>::filter()
{
    return _fwt.filter();
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline void TWavelet2DOperator<T, TFilter>::apply(const arma::Col<T>& in, arma::Col<T>& out)
{
    ThrowTraceExceptionIf(KIncompatibleSizeOperatorException, in.n_rows!=this->n());

    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(_height, _width);

    _fwt.transform(inMatrix, outMatrix);

    outMatrix.reshape(this->m(), 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline void TWavelet2DOperator<T, TFilter>::applyAdjoint(const arma::Col<T>& in, arma::Col<T>& out)
{
    ThrowTraceExceptionIf(KIncompatibleSizeOperatorException, in.n_rows!=this->m());

    TWavelet2DOperatorSpecialisation<T, TFilter>::ApplyAdjoint(in, out, _height, _width, _fwt);
}

// ---------------------------------------------------------------------------------------------------- //

template<class T, class TFilter>
inline void TWavelet2DOperatorSpecialisation<T, TFilter>::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, TFilter>& fwt)
{
    arma::Mat<T> fwtMatrix(in.n_rows, in.n_rows);
    for(klab::UInt32 j=0; j<in.n_rows; ++j)
    {
        arma::Mat<T> base(in.n_rows, 1);
        arma::Mat<T> col;
        base.fill(0.0);
        base[j] = 1.0;
        base.reshape(height, width);

        fwt.transform(base, col);

        col.reshape(width*height, 1);
        for(klab::UInt32 i=0; i<col.n_rows; ++i)
            fwtMatrix(i, j) = col[i];
    } 

    out = arma::trans(fwtMatrix) * in;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::THaarWaveletFilter<T> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::THaarWaveletFilter<T> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 2> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 2> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 4> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 4> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 6> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 6> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 8> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 8> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

template<class T>
inline void TWavelet2DOperatorSpecialisation<T, klab::TDaubechiesWaveletFilter<T, 10> >::ApplyAdjoint(const arma::Col<T>& in, arma::Col<T>& out, klab::UInt32 height, klab::UInt32 width, klab::TFWT2D<T, klab::TDaubechiesWaveletFilter<T, 10> >& fwt)
{
    arma::Mat<T> inMatrix(in);
    arma::Mat<T> outMatrix;
    inMatrix.reshape(height, width);

    fwt.invert(inMatrix, outMatrix);

    outMatrix.reshape(width*height, 1);
    out = outMatrix;
}

// ---------------------------------------------------------------------------------------------------- //

}


#endif
