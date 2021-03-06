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

#include "InverseFourier1DOperatorUnitTest.h"
#include "../include/InverseFourier1DOperator.h"

using namespace klab;
using namespace kl1p;




// ---------------------------------------------------------------------------------------------------- //

KInverseFourier1DOperatorUnitTest::KInverseFourier1DOperatorUnitTest() : Test::Suite()
{ 
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testNConstructor)  
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testNShiftConstructor)  
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testCopyConstructor)

    TEST_ADD(KInverseFourier1DOperatorUnitTest::testIsShift)
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testApply)
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testApplyAdjoint)
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testApplyAdjoint_Real)
    TEST_ADD(KInverseFourier1DOperatorUnitTest::testApplyAdjoint_Real_02)
}

// ---------------------------------------------------------------------------------------------------- //

KInverseFourier1DOperatorUnitTest::~KInverseFourier1DOperatorUnitTest()
{}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testNConstructor()
{
    try
    {
        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(128);
        TEST_ASSERT(op.m()==128 && op.n()==128)
        TEST_ASSERT(op.isShift()==false)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testNShiftConstructor()
{
    try
    {
        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(128, true);
        TEST_ASSERT(op.m()==128 && op.n()==128)
        TEST_ASSERT(op.isShift()==true)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testCopyConstructor()
{
    try
    {
        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(128, true);
        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op2(op);
        TEST_ASSERT(op2.m()==128 && op2.n()==128)
        TEST_ASSERT(op2.isShift() == true)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testIsShift()
{
    try
    {
        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(0);
        TEST_ASSERT(op.isShift() == false)

        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op2(128, true);
        TEST_ASSERT(op2.isShift() == true)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testApply()
{
    try
    {
        klab::DoubleReal precision = 1e-12;

        arma::Col<std::complex<klab::DoubleReal> > x(8);
        arma::Col<std::complex<klab::DoubleReal> > y;

        x[0] = std::complex<klab::DoubleReal>(36.0, 0.0);
        x[1] = std::complex<klab::DoubleReal>(-4.0, 9.656854249492380);
        x[2] = std::complex<klab::DoubleReal>(-4.0, 4.0);
        x[3] = std::complex<klab::DoubleReal>(-4.0, 1.656854249492381);
        x[4] = std::complex<klab::DoubleReal>(-4.0, 0.0);
        x[5] = std::complex<klab::DoubleReal>(-4.0, -1.656854249492381);
        x[6] = std::complex<klab::DoubleReal>(-4.0, -4.0);
        x[7] = std::complex<klab::DoubleReal>(-4.0, -9.656854249492380);


        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(8);
        op.apply(x, y);

        TEST_ASSERT(y.n_rows == 8)
        TEST_ASSERT(klab::Equals(y[0], std::complex<klab::DoubleReal>(1.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[1], std::complex<klab::DoubleReal>(2.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[2], std::complex<klab::DoubleReal>(3.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[3], std::complex<klab::DoubleReal>(4.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[4], std::complex<klab::DoubleReal>(5.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[5], std::complex<klab::DoubleReal>(6.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[6], std::complex<klab::DoubleReal>(7.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[7], std::complex<klab::DoubleReal>(8.0, 0.0), precision))
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testApplyAdjoint()
{
    try
    {
        klab::DoubleReal precision = 1e-12;

        arma::Col<std::complex<klab::DoubleReal> > x(8);
        arma::Col<std::complex<klab::DoubleReal> > y;
        arma::Col<std::complex<klab::DoubleReal> > xr;

        x[0] = std::complex<klab::DoubleReal>(36.0, 0.0);
        x[1] = std::complex<klab::DoubleReal>(-4.0, 9.656854249492380);
        x[2] = std::complex<klab::DoubleReal>(-4.0, 4.0);
        x[3] = std::complex<klab::DoubleReal>(-4.0, 1.656854249492381);
        x[4] = std::complex<klab::DoubleReal>(-4.0, 0.0);
        x[5] = std::complex<klab::DoubleReal>(-4.0, -1.656854249492381);
        x[6] = std::complex<klab::DoubleReal>(-4.0, -4.0);
        x[7] = std::complex<klab::DoubleReal>(-4.0, -9.656854249492380);
        

        klab::TFFT1D<std::complex<klab::DoubleReal> > fft;
        arma::Mat<std::complex<klab::DoubleReal> > fftMatrix(8, 8);
        for(klab::UInt32 j=0; j<8; ++j)
        {
            arma::Col<std::complex<klab::DoubleReal> > base(8);
            arma::Col<std::complex<klab::DoubleReal> > col;
            base.fill(0.0);
            base[j] = 1.0;

            fft.invert(base, col);

            for(klab::UInt32 i=0; i<col.n_rows; ++i)
                fftMatrix(i, j) = col[i];
        } 

        TInverseFourier1DOperator<std::complex<klab::DoubleReal> > op(8);

        op.apply(x, y);
        TEST_ASSERT(y.n_rows == 8)
        TEST_ASSERT(klab::Equals(y[0], std::complex<klab::DoubleReal>(1.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[1], std::complex<klab::DoubleReal>(2.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[2], std::complex<klab::DoubleReal>(3.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[3], std::complex<klab::DoubleReal>(4.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[4], std::complex<klab::DoubleReal>(5.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[5], std::complex<klab::DoubleReal>(6.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[6], std::complex<klab::DoubleReal>(7.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(y[7], std::complex<klab::DoubleReal>(8.0, 0.0), precision))

        arma::Col<std::complex<klab::DoubleReal> > xt = arma::trans(fftMatrix) * y;

        op.applyAdjoint(y, xr);
        TEST_ASSERT(xr.n_rows == 8)
        TEST_ASSERT(klab::Equals(xr[0], xt[0], precision))
        TEST_ASSERT(klab::Equals(xr[1], xt[1], precision))
        TEST_ASSERT(klab::Equals(xr[2], xt[2], precision))
        TEST_ASSERT(klab::Equals(xr[3], xt[3], precision))
        TEST_ASSERT(klab::Equals(xr[4], xt[4], precision))
        TEST_ASSERT(klab::Equals(xr[5], xt[5], precision))
        TEST_ASSERT(klab::Equals(xr[6], xt[6], precision))
        TEST_ASSERT(klab::Equals(xr[7], xt[7], precision))
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testApplyAdjoint_Real()
{
    try
    {
        klab::DoubleReal precision = 1e-12;

        arma::Col<klab::DoubleReal> x(8);
        arma::Col<klab::DoubleReal> y;
        arma::Col<klab::DoubleReal> xr;

        x[0] =  36.0;
        x[1] = -4.0;
        x[2] = -4.0;
        x[3] = -4.0;
        x[4] = -4.0;
        x[5] =  9.65685424949238;
        x[6] =  4.0;
        x[7] =  1.65685424949238;        

        klab::TFFT1D<klab::DoubleReal> fft;
        arma::Mat<klab::DoubleReal> fftMatrix(8, 8);
        for(klab::UInt32 j=0; j<8; ++j)
        {
            arma::Col<klab::DoubleReal> base(8);
            arma::Col<klab::DoubleReal> col;
            base.fill(0.0);
            base[j] = 1.0;

            fft.invert(base, col);

            for(klab::UInt32 i=0; i<col.n_rows; ++i)
                fftMatrix(i, j) = col[i];
        } 

        arma::Col<klab::DoubleReal> xt = arma::trans(fftMatrix) * x;

        TInverseFourier1DOperator<klab::DoubleReal> op(8);

        op.applyAdjoint(x, xr);
        TEST_ASSERT(xr.n_rows == 8)
        TEST_ASSERT(klab::Equals(xr[0], xt[0], precision))
        TEST_ASSERT(klab::Equals(xr[1], xt[1], precision))
        TEST_ASSERT(klab::Equals(xr[2], xt[2], precision))
        TEST_ASSERT(klab::Equals(xr[3], xt[3], precision))
        TEST_ASSERT(klab::Equals(xr[4], xt[4], precision))
        TEST_ASSERT(klab::Equals(xr[5], xt[5], precision))
        TEST_ASSERT(klab::Equals(xr[6], xt[6], precision))
        TEST_ASSERT(klab::Equals(xr[7], xt[7], precision))
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KInverseFourier1DOperatorUnitTest::testApplyAdjoint_Real_02()
{
    try
    {
        klab::DoubleReal precision = 1e-12;

        arma::Col<klab::DoubleReal> x(9);
        arma::Col<klab::DoubleReal> y;
        arma::Col<klab::DoubleReal> xr;

        x[0] = 1.0;
        x[1] = 2.0;
        x[2] = 3.0;
        x[3] = 4.0;
        x[4] = 5.0;
        x[5] = 6.0;
        x[6] = 7.0;
        x[7] = 8.0;   
        x[8] = 9.0;   

        klab::TFFT1D<klab::DoubleReal> fft;
        arma::Mat<klab::DoubleReal> fftMatrix(9, 9);
        for(klab::UInt32 j=0; j<9; ++j)
        {
            arma::Col<klab::DoubleReal> base(9);
            arma::Col<klab::DoubleReal> col;
            base.fill(0.0);
            base[j] = 1.0;

            fft.invert(base, col);

            for(klab::UInt32 i=0; i<col.n_rows; ++i)
                fftMatrix(i, j) = col[i];
        } 

        arma::Col<klab::DoubleReal> xt = arma::trans(fftMatrix) * x;

        TInverseFourier1DOperator<klab::DoubleReal> op(9);

        op.applyAdjoint(x, xr);
        TEST_ASSERT(xr.n_rows == 9)
        TEST_ASSERT(klab::Equals(xr[0], xt[0], precision))
        TEST_ASSERT(klab::Equals(xr[1], xt[1], precision))
        TEST_ASSERT(klab::Equals(xr[2], xt[2], precision))
        TEST_ASSERT(klab::Equals(xr[3], xt[3], precision))
        TEST_ASSERT(klab::Equals(xr[4], xt[4], precision))
        TEST_ASSERT(klab::Equals(xr[5], xt[5], precision))
        TEST_ASSERT(klab::Equals(xr[6], xt[6], precision))
        TEST_ASSERT(klab::Equals(xr[7], xt[7], precision))
        TEST_ASSERT(klab::Equals(xr[8], xt[8], precision))
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

