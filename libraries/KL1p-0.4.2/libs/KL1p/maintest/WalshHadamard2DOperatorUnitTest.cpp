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

#include "WalshHadamard2DOperatorUnitTest.h"
#include "../include/WalshHadamard2DOperator.h"

using namespace klab;
using namespace kl1p;




// ---------------------------------------------------------------------------------------------------- //

KWalshHadamard2DOperatorUnitTest::KWalshHadamard2DOperatorUnitTest() : Test::Suite()
{
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testHeightWidthConstructor) 
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testCopyConstructor)

    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testWidth) 
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testHeight) 
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testApply)
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testApply_Complex)
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testApplyAdjoint)
    TEST_ADD(KWalshHadamard2DOperatorUnitTest::testApplyAdjoint_Complex)
}

// ---------------------------------------------------------------------------------------------------- //

KWalshHadamard2DOperatorUnitTest::~KWalshHadamard2DOperatorUnitTest()
{}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testHeightWidthConstructor()
{
    try
    {
        TWalshHadamard2DOperator<klab::DoubleReal> op1(32, 64);
        TEST_ASSERT(op1.m()==2048 && op1.n()==2048)
        TEST_ASSERT(op1.height()==32 && op1.width()==64)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op2(64, 32);
        TEST_ASSERT(op2.m()==2048 && op2.n()==2048)
        TEST_ASSERT(op2.height()==64 && op2.width()==32)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testCopyConstructor()
{
    try
    {
        TWalshHadamard2DOperator<klab::DoubleReal> op1(32, 64);
        TWalshHadamard2DOperator<klab::DoubleReal> op1Copy(op1);
        TEST_ASSERT(op1Copy.m()==2048 && op1Copy.n()==2048)
        TEST_ASSERT(op1Copy.height()==32 && op1Copy.width()==64)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op2(64, 32);
        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op2Copy(op2);
        TEST_ASSERT(op2Copy.m()==2048 && op2Copy.n()==2048)
        TEST_ASSERT(op2Copy.height()==64 && op2Copy.width()==32)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testWidth()
{
    try
    {
        TWalshHadamard2DOperator<klab::DoubleReal> op1(0, 0);
        TEST_ASSERT(op1.height()==0 && op1.width()==0)

        TWalshHadamard2DOperator<klab::DoubleReal> op2(32, 64);
        TEST_ASSERT(op2.height()==32 && op2.width()==64)

        TWalshHadamard2DOperator<klab::DoubleReal> op3(64, 32);
        TEST_ASSERT(op3.height()==64 && op3.width()==32)


        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op4(0, 0);
        TEST_ASSERT(op4.height()==0 && op4.width()==0)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op5(32, 64);
        TEST_ASSERT(op5.height()==32 && op5.width()==64)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op6(64, 32);
        TEST_ASSERT(op6.height()==64 && op6.width()==32)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testHeight()
{
    try
    {
        TWalshHadamard2DOperator<klab::DoubleReal> op1(0, 0);
        TEST_ASSERT(op1.height()==0 && op1.width()==0)

        TWalshHadamard2DOperator<klab::DoubleReal> op2(32, 64);
        TEST_ASSERT(op2.height()==32 && op2.width()==64)

        TWalshHadamard2DOperator<klab::DoubleReal> op3(64, 32);
        TEST_ASSERT(op3.height()==64 && op3.width()==32)


        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op4(0, 0);
        TEST_ASSERT(op4.height()==0 && op4.width()==0)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op5(32, 64);
        TEST_ASSERT(op5.height()==32 && op5.width()==64)

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op6(64, 32);
        TEST_ASSERT(op6.height()==64 && op6.width()==32)
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testApply()
{
    try
    {
        klab::UInt32 m = 8;
        klab::UInt32 n = 8;
        klab::DoubleReal precision = 1e-12;

        arma::Mat<klab::DoubleReal> f(m, n);
        arma::Mat<klab::DoubleReal> F;

        f(0, 0) = 150.0;
        f(0, 1) = 111.0;
        f(0, 2) = 129.0;
        f(0, 3) = 139.0;
        f(0, 4) = 127.0;
        f(0, 5) = 143.0;
        f(0, 6) = 160.0;
        f(0, 7) = 90.0;
        f(1, 0) = 127.0;
        f(1, 1) = 108.0;
        f(1, 2) = 120.0;
        f(1, 3) = 167.0;
        f(1, 4) = 182.0;
        f(1, 5) = 136.0;
        f(1, 6) = 118.0;
        f(1, 7) = 85.0;
        f(2, 0) = 124.0;
        f(2, 1) = 114.0;
        f(2, 2) = 132.0;
        f(2, 3) = 141.0;
        f(2, 4) = 171.0;
        f(2, 5) = 184.0;
        f(2, 6) = 110.0;
        f(2, 7) = 141.0;
        f(3, 0) = 127.0;
        f(3, 1) = 115.0;
        f(3, 2) = 107.0;
        f(3, 3) = 85.0;
        f(3, 4) = 161.0;
        f(3, 5) = 122.0;
        f(3, 6) = 97.0;
        f(3, 7) = 157.0;
        f(4, 0) = 121.0;
        f(4, 1) = 107.0;
        f(4, 2) = 73.0;
        f(4, 3) = 118.0;
        f(4, 4) = 153.0;
        f(4, 5) = 97.0;
        f(4, 6) = 126.0;
        f(4, 7) = 170.0;
        f(5, 0) = 115.0;
        f(5, 1) = 98.0;
        f(5, 2) = 90.0;
        f(5, 3) = 91.0;
        f(5, 4) = 140.0;
        f(5, 5) = 87.0;
        f(5, 6) = 144.0;
        f(5, 7) = 195.0;
        f(6, 0) = 120.0;
        f(6, 1) = 81.0;
        f(6, 2) = 84.0;
        f(6, 3) = 77.0;
        f(6, 4) = 151.0;
        f(6, 5) = 138.0;
        f(6, 6) = 144.0;
        f(6, 7) = 135.0;
        f(7, 0) = 126.0;
        f(7, 1) = 71.0;
        f(7, 2) = 70.0;
        f(7, 3) = 108.0;
        f(7, 4) = 149.0;
        f(7, 5) = 175.0;
        f(7, 6) = 125.0;
        f(7, 7) = 93.0;

        f.reshape(64, 1);
        arma::Col<klab::DoubleReal> x(f);
        arma::Col<klab::DoubleReal> y;
        f.reshape(8, 8);

        TWalshHadamard2DOperator<klab::DoubleReal> op(8, 8);
        op.apply(x, y);
        F = y;
        F.reshape(8, 8);

        TEST_ASSERT(F.n_rows==m && F.n_cols==n)
        TEST_ASSERT(klab::Equals(F(0, 0), 7952.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 1), 194.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 2), 310.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 3), 520.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 4), -860.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 5), -26.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 6), -142.0, precision))
        TEST_ASSERT(klab::Equals(F(0, 7), 132.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 0), 170.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 1), -16.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 2), -64.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 3), -130.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 4), 22.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 5), 28.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 6), 36.0, precision))
        TEST_ASSERT(klab::Equals(F(1, 7), 114.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 0), 82.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 1), 72.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 2), -336.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 3), 126.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 4), 282.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 5), -296.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 6), 208.0, precision))
        TEST_ASSERT(klab::Equals(F(2, 7), -10.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 0), -148.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 1), 6.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 2), 106.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 3), -28.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 4), 84.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 5), 22.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 6), 122.0, precision))
        TEST_ASSERT(klab::Equals(F(3, 7), 140.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 0), 408.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 1), 14.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 2), 138.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 3), -184.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 4), 484.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 5), -38.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 6), -482.0, precision))
        TEST_ASSERT(klab::Equals(F(4, 7), 28.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 0), 134.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 1), -32.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 2), -56.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 3), -206.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 4), 10.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 5), 116.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 6), -76.0, precision))
        TEST_ASSERT(klab::Equals(F(5, 7), 270.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 0), -74.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 1), 256.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 2), 192.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 3), -294.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 4), 134.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 5), -168.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 6), -360.0, precision))
        TEST_ASSERT(klab::Equals(F(6, 7), 602.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 0), -132.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 1), 170.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 2), -186.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 3), -100.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 4), -84.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 5), 162.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 6), 478.0, precision))
        TEST_ASSERT(klab::Equals(F(7, 7), -196.0, precision))


        // Limit cases.
        try                                         {op.apply(arma::Col<klab::DoubleReal>(), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<klab::DoubleReal>(56), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<klab::DoubleReal>(49), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<klab::DoubleReal>(63), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testApply_Complex()
{
    try
    {
        klab::UInt32 m = 8;
        klab::UInt32 n = 8;
        klab::DoubleReal precision = 1e-12;

        arma::Mat<std::complex<klab::DoubleReal> > f(m, n);
        arma::Mat<std::complex<klab::DoubleReal> > F;

        f(0, 0) = std::complex<klab::DoubleReal>(150.0, 0.0);
        f(0, 1) = std::complex<klab::DoubleReal>(111.0, 0.0);
        f(0, 2) = std::complex<klab::DoubleReal>(129.0, 0.0);
        f(0, 3) = std::complex<klab::DoubleReal>(139.0, 0.0);
        f(0, 4) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(0, 5) = std::complex<klab::DoubleReal>(143.0, 0.0);
        f(0, 6) = std::complex<klab::DoubleReal>(160.0, 0.0);
        f(0, 7) = std::complex<klab::DoubleReal>(90.0, 0.0);
        f(1, 0) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(1, 1) = std::complex<klab::DoubleReal>(108.0, 0.0);
        f(1, 2) = std::complex<klab::DoubleReal>(120.0, 0.0);
        f(1, 3) = std::complex<klab::DoubleReal>(167.0, 0.0);
        f(1, 4) = std::complex<klab::DoubleReal>(182.0, 0.0);
        f(1, 5) = std::complex<klab::DoubleReal>(136.0, 0.0);
        f(1, 6) = std::complex<klab::DoubleReal>(118.0, 0.0);
        f(1, 7) = std::complex<klab::DoubleReal>(85.0, 0.0);
        f(2, 0) = std::complex<klab::DoubleReal>(124.0, 0.0);
        f(2, 1) = std::complex<klab::DoubleReal>(114.0, 0.0);
        f(2, 2) = std::complex<klab::DoubleReal>(132.0, 0.0);
        f(2, 3) = std::complex<klab::DoubleReal>(141.0, 0.0);
        f(2, 4) = std::complex<klab::DoubleReal>(171.0, 0.0);
        f(2, 5) = std::complex<klab::DoubleReal>(184.0, 0.0);
        f(2, 6) = std::complex<klab::DoubleReal>(110.0, 0.0);
        f(2, 7) = std::complex<klab::DoubleReal>(141.0, 0.0);
        f(3, 0) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(3, 1) = std::complex<klab::DoubleReal>(115.0, 0.0);
        f(3, 2) = std::complex<klab::DoubleReal>(107.0, 0.0);
        f(3, 3) = std::complex<klab::DoubleReal>(85.0, 0.0);
        f(3, 4) = std::complex<klab::DoubleReal>(161.0, 0.0);
        f(3, 5) = std::complex<klab::DoubleReal>(122.0, 0.0);
        f(3, 6) = std::complex<klab::DoubleReal>(97.0, 0.0);
        f(3, 7) = std::complex<klab::DoubleReal>(157.0, 0.0);
        f(4, 0) = std::complex<klab::DoubleReal>(121.0, 0.0);
        f(4, 1) = std::complex<klab::DoubleReal>(107.0, 0.0);
        f(4, 2) = std::complex<klab::DoubleReal>(73.0, 0.0);
        f(4, 3) = std::complex<klab::DoubleReal>(118.0, 0.0);
        f(4, 4) = std::complex<klab::DoubleReal>(153.0, 0.0);
        f(4, 5) = std::complex<klab::DoubleReal>(97.0, 0.0);
        f(4, 6) = std::complex<klab::DoubleReal>(126.0, 0.0);
        f(4, 7) = std::complex<klab::DoubleReal>(170.0, 0.0);
        f(5, 0) = std::complex<klab::DoubleReal>(115.0, 0.0);
        f(5, 1) = std::complex<klab::DoubleReal>(98.0, 0.0);
        f(5, 2) = std::complex<klab::DoubleReal>(90.0, 0.0);
        f(5, 3) = std::complex<klab::DoubleReal>(91.0, 0.0);
        f(5, 4) = std::complex<klab::DoubleReal>(140.0, 0.0);
        f(5, 5) = std::complex<klab::DoubleReal>(87.0, 0.0);
        f(5, 6) = std::complex<klab::DoubleReal>(144.0, 0.0);
        f(5, 7) = std::complex<klab::DoubleReal>(195.0, 0.0);
        f(6, 0) = std::complex<klab::DoubleReal>(120.0, 0.0);
        f(6, 1) = std::complex<klab::DoubleReal>(81.0, 0.0);
        f(6, 2) = std::complex<klab::DoubleReal>(84.0, 0.0);
        f(6, 3) = std::complex<klab::DoubleReal>(77.0, 0.0);
        f(6, 4) = std::complex<klab::DoubleReal>(151.0, 0.0);
        f(6, 5) = std::complex<klab::DoubleReal>(138.0, 0.0);
        f(6, 6) = std::complex<klab::DoubleReal>(144.0, 0.0);
        f(6, 7) = std::complex<klab::DoubleReal>(135.0, 0.0);
        f(7, 0) = std::complex<klab::DoubleReal>(126.0, 0.0);
        f(7, 1) = std::complex<klab::DoubleReal>(71.0, 0.0);
        f(7, 2) = std::complex<klab::DoubleReal>(70.0, 0.0);
        f(7, 3) = std::complex<klab::DoubleReal>(108.0, 0.0);
        f(7, 4) = std::complex<klab::DoubleReal>(149.0, 0.0);
        f(7, 5) = std::complex<klab::DoubleReal>(175.0, 0.0);
        f(7, 6) = std::complex<klab::DoubleReal>(125.0, 0.0);
        f(7, 7) = std::complex<klab::DoubleReal>(93.0, 0.0);

        f.reshape(64, 1);
        arma::Col<std::complex<klab::DoubleReal> > x(f);
        arma::Col<std::complex<klab::DoubleReal> > y;
        f.reshape(8, 8);

        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op(8, 8);
        op.apply(x, y);
        F = y;
        F.reshape(8, 8);

        TEST_ASSERT(F.n_rows==m && F.n_cols==n)
        TEST_ASSERT(klab::Equals(F(0, 0), std::complex<klab::DoubleReal>(7952.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 1), std::complex<klab::DoubleReal>(194.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 2), std::complex<klab::DoubleReal>(310.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 3), std::complex<klab::DoubleReal>(520.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 4), std::complex<klab::DoubleReal>(-860.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 5), std::complex<klab::DoubleReal>(-26.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 6), std::complex<klab::DoubleReal>(-142.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(0, 7), std::complex<klab::DoubleReal>(132.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 0), std::complex<klab::DoubleReal>(170.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 1), std::complex<klab::DoubleReal>(-16.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 2), std::complex<klab::DoubleReal>(-64.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 3), std::complex<klab::DoubleReal>(-130.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 4), std::complex<klab::DoubleReal>(22.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 5), std::complex<klab::DoubleReal>(28.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 6), std::complex<klab::DoubleReal>(36.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(1, 7), std::complex<klab::DoubleReal>(114.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 0), std::complex<klab::DoubleReal>(82.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 1), std::complex<klab::DoubleReal>(72.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 2), std::complex<klab::DoubleReal>(-336.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 3), std::complex<klab::DoubleReal>(126.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 4), std::complex<klab::DoubleReal>(282.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 5), std::complex<klab::DoubleReal>(-296.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 6), std::complex<klab::DoubleReal>(208.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(2, 7), std::complex<klab::DoubleReal>(-10.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 0), std::complex<klab::DoubleReal>(-148.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 1), std::complex<klab::DoubleReal>(6.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 2), std::complex<klab::DoubleReal>(106.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 3), std::complex<klab::DoubleReal>(-28.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 4), std::complex<klab::DoubleReal>(84.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 5), std::complex<klab::DoubleReal>(22.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 6), std::complex<klab::DoubleReal>(122.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(3, 7), std::complex<klab::DoubleReal>(140.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 0), std::complex<klab::DoubleReal>(408.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 1), std::complex<klab::DoubleReal>(14.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 2), std::complex<klab::DoubleReal>(138.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 3), std::complex<klab::DoubleReal>(-184.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 4), std::complex<klab::DoubleReal>(484.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 5), std::complex<klab::DoubleReal>(-38.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 6), std::complex<klab::DoubleReal>(-482.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(4, 7), std::complex<klab::DoubleReal>(28.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 0), std::complex<klab::DoubleReal>(134.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 1), std::complex<klab::DoubleReal>(-32.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 2), std::complex<klab::DoubleReal>(-56.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 3), std::complex<klab::DoubleReal>(-206.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 4), std::complex<klab::DoubleReal>(10.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 5), std::complex<klab::DoubleReal>(116.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 6), std::complex<klab::DoubleReal>(-76.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(5, 7), std::complex<klab::DoubleReal>(270.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 0), std::complex<klab::DoubleReal>(-74.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 1), std::complex<klab::DoubleReal>(256.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 2), std::complex<klab::DoubleReal>(192.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 3), std::complex<klab::DoubleReal>(-294.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 4), std::complex<klab::DoubleReal>(134.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 5), std::complex<klab::DoubleReal>(-168.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 6), std::complex<klab::DoubleReal>(-360.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(6, 7), std::complex<klab::DoubleReal>(602.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 0), std::complex<klab::DoubleReal>(-132.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 1), std::complex<klab::DoubleReal>(170.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 2), std::complex<klab::DoubleReal>(-186.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 3), std::complex<klab::DoubleReal>(-100.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 4), std::complex<klab::DoubleReal>(-84.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 5), std::complex<klab::DoubleReal>(162.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 6), std::complex<klab::DoubleReal>(478.0, 0.0), precision))
        TEST_ASSERT(klab::Equals(F(7, 7), std::complex<klab::DoubleReal>(-196.0, 0.0), precision))


        // Limit cases.
        try                                         {op.apply(arma::Col<std::complex<klab::DoubleReal> >(), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<std::complex<klab::DoubleReal> >(56), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<std::complex<klab::DoubleReal> >(49), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.apply(arma::Col<std::complex<klab::DoubleReal> >(63), y); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testApplyAdjoint()
{
    try
    {
        klab::UInt32 m = 8;
        klab::UInt32 n = 8;
        klab::DoubleReal precision = 1e-12;

        arma::Mat<klab::DoubleReal> f(m, n);
        arma::Mat<klab::DoubleReal> F;

        f(0, 0) = 150.0;
        f(0, 1) = 111.0;
        f(0, 2) = 129.0;
        f(0, 3) = 139.0;
        f(0, 4) = 127.0;
        f(0, 5) = 143.0;
        f(0, 6) = 160.0;
        f(0, 7) = 90.0;
        f(1, 0) = 127.0;
        f(1, 1) = 108.0;
        f(1, 2) = 120.0;
        f(1, 3) = 167.0;
        f(1, 4) = 182.0;
        f(1, 5) = 136.0;
        f(1, 6) = 118.0;
        f(1, 7) = 85.0;
        f(2, 0) = 124.0;
        f(2, 1) = 114.0;
        f(2, 2) = 132.0;
        f(2, 3) = 141.0;
        f(2, 4) = 171.0;
        f(2, 5) = 184.0;
        f(2, 6) = 110.0;
        f(2, 7) = 141.0;
        f(3, 0) = 127.0;
        f(3, 1) = 115.0;
        f(3, 2) = 107.0;
        f(3, 3) = 85.0;
        f(3, 4) = 161.0;
        f(3, 5) = 122.0;
        f(3, 6) = 97.0;
        f(3, 7) = 157.0;
        f(4, 0) = 121.0;
        f(4, 1) = 107.0;
        f(4, 2) = 73.0;
        f(4, 3) = 118.0;
        f(4, 4) = 153.0;
        f(4, 5) = 97.0;
        f(4, 6) = 126.0;
        f(4, 7) = 170.0;
        f(5, 0) = 115.0;
        f(5, 1) = 98.0;
        f(5, 2) = 90.0;
        f(5, 3) = 91.0;
        f(5, 4) = 140.0;
        f(5, 5) = 87.0;
        f(5, 6) = 144.0;
        f(5, 7) = 195.0;
        f(6, 0) = 120.0;
        f(6, 1) = 81.0;
        f(6, 2) = 84.0;
        f(6, 3) = 77.0;
        f(6, 4) = 151.0;
        f(6, 5) = 138.0;
        f(6, 6) = 144.0;
        f(6, 7) = 135.0;
        f(7, 0) = 126.0;
        f(7, 1) = 71.0;
        f(7, 2) = 70.0;
        f(7, 3) = 108.0;
        f(7, 4) = 149.0;
        f(7, 5) = 175.0;
        f(7, 6) = 125.0;
        f(7, 7) = 93.0;

        f.reshape(64, 1);
        arma::Col<klab::DoubleReal> x(f);
        arma::Col<klab::DoubleReal> y;
        arma::Col<klab::DoubleReal> xr;
        f.reshape(8, 8);

        klab::TFWHT2D<klab::DoubleReal> fwht;
        arma::Mat<klab::DoubleReal> fwhtMatrix(64, 64);
        for(klab::UInt32 j=0; j<64; ++j)
        {
            arma::Col<klab::DoubleReal> base(64);
            arma::Col<klab::DoubleReal> col;
            base.fill(0.0);
            base[j] = 1.0;

            arma::Mat<klab::DoubleReal> inMatrix(base);
            arma::Mat<klab::DoubleReal> outMatrix;
            inMatrix.reshape(8, 8);
            fwht.transform(inMatrix, outMatrix);
            outMatrix.reshape(64, 1);

            for(klab::UInt32 i=0; i<outMatrix.n_rows; ++i)
                fwhtMatrix(i, j) = outMatrix[i];
        } 


        TWalshHadamard2DOperator<klab::DoubleReal> op(8, 8);
        op.apply(x, y);

        arma::Col<klab::DoubleReal> xt = arma::trans(fwhtMatrix) * y;

        op.applyAdjoint(y, xr);
        F = xr;
        F.reshape(8, 8);

        TEST_ASSERT(xr.n_rows==64)
        TEST_ASSERT(F.n_rows==8 && F.n_cols==8)
        TEST_ASSERT(klab::Equals(xr[0], xt[0], precision))
        TEST_ASSERT(klab::Equals(xr[1], xt[1], precision))
        TEST_ASSERT(klab::Equals(xr[2], xt[2], precision))
        TEST_ASSERT(klab::Equals(xr[3], xt[3], precision))
        TEST_ASSERT(klab::Equals(xr[4], xt[4], precision))
        TEST_ASSERT(klab::Equals(xr[5], xt[5], precision))
        TEST_ASSERT(klab::Equals(xr[6], xt[6], precision))
        TEST_ASSERT(klab::Equals(xr[7], xt[7], precision))
        TEST_ASSERT(klab::Equals(xr[8], xt[8], precision))
        TEST_ASSERT(klab::Equals(xr[9], xt[9], precision))
        TEST_ASSERT(klab::Equals(xr[10], xt[10], precision))
        TEST_ASSERT(klab::Equals(xr[11], xt[11], precision))
        TEST_ASSERT(klab::Equals(xr[12], xt[12], precision))
        TEST_ASSERT(klab::Equals(xr[13], xt[13], precision))
        TEST_ASSERT(klab::Equals(xr[14], xt[14], precision))
        TEST_ASSERT(klab::Equals(xr[15], xt[15], precision))
        TEST_ASSERT(klab::Equals(xr[16], xt[16], precision))
        TEST_ASSERT(klab::Equals(xr[17], xt[17], precision))
        TEST_ASSERT(klab::Equals(xr[18], xt[18], precision))
        TEST_ASSERT(klab::Equals(xr[19], xt[19], precision))
        TEST_ASSERT(klab::Equals(xr[20], xt[20], precision))
        TEST_ASSERT(klab::Equals(xr[21], xt[21], precision))
        TEST_ASSERT(klab::Equals(xr[22], xt[22], precision))
        TEST_ASSERT(klab::Equals(xr[23], xt[23], precision))
        TEST_ASSERT(klab::Equals(xr[24], xt[24], precision))
        TEST_ASSERT(klab::Equals(xr[25], xt[25], precision))
        TEST_ASSERT(klab::Equals(xr[26], xt[26], precision))
        TEST_ASSERT(klab::Equals(xr[27], xt[27], precision))
        TEST_ASSERT(klab::Equals(xr[28], xt[28], precision))
        TEST_ASSERT(klab::Equals(xr[29], xt[29], precision))
        TEST_ASSERT(klab::Equals(xr[30], xt[30], precision))
        TEST_ASSERT(klab::Equals(xr[31], xt[31], precision))
        TEST_ASSERT(klab::Equals(xr[32], xt[32], precision))
        TEST_ASSERT(klab::Equals(xr[33], xt[33], precision))
        TEST_ASSERT(klab::Equals(xr[34], xt[34], precision))
        TEST_ASSERT(klab::Equals(xr[35], xt[35], precision))
        TEST_ASSERT(klab::Equals(xr[36], xt[36], precision))
        TEST_ASSERT(klab::Equals(xr[37], xt[37], precision))
        TEST_ASSERT(klab::Equals(xr[38], xt[38], precision))
        TEST_ASSERT(klab::Equals(xr[39], xt[39], precision))
        TEST_ASSERT(klab::Equals(xr[40], xt[40], precision))
        TEST_ASSERT(klab::Equals(xr[41], xt[41], precision))
        TEST_ASSERT(klab::Equals(xr[42], xt[42], precision))
        TEST_ASSERT(klab::Equals(xr[43], xt[43], precision))
        TEST_ASSERT(klab::Equals(xr[44], xt[44], precision))
        TEST_ASSERT(klab::Equals(xr[45], xt[45], precision))
        TEST_ASSERT(klab::Equals(xr[46], xt[46], precision))
        TEST_ASSERT(klab::Equals(xr[47], xt[47], precision))
        TEST_ASSERT(klab::Equals(xr[48], xt[48], precision))
        TEST_ASSERT(klab::Equals(xr[49], xt[49], precision))
        TEST_ASSERT(klab::Equals(xr[50], xt[50], precision))
        TEST_ASSERT(klab::Equals(xr[51], xt[51], precision))
        TEST_ASSERT(klab::Equals(xr[52], xt[52], precision))
        TEST_ASSERT(klab::Equals(xr[53], xt[53], precision))
        TEST_ASSERT(klab::Equals(xr[54], xt[54], precision))
        TEST_ASSERT(klab::Equals(xr[55], xt[55], precision))
        TEST_ASSERT(klab::Equals(xr[56], xt[56], precision))
        TEST_ASSERT(klab::Equals(xr[57], xt[57], precision))
        TEST_ASSERT(klab::Equals(xr[58], xt[58], precision))
        TEST_ASSERT(klab::Equals(xr[59], xt[59], precision))
        TEST_ASSERT(klab::Equals(xr[60], xt[60], precision))
        TEST_ASSERT(klab::Equals(xr[61], xt[61], precision))
        TEST_ASSERT(klab::Equals(xr[62], xt[62], precision))
        TEST_ASSERT(klab::Equals(xr[63], xt[63], precision))


        // Limit cases.
        try                                         {op.applyAdjoint(arma::Col<klab::DoubleReal>(), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<klab::DoubleReal>(56), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<klab::DoubleReal>(49), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<klab::DoubleReal>(63), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //

void KWalshHadamard2DOperatorUnitTest::testApplyAdjoint_Complex()
{
    try
    {
        klab::UInt32 m = 8;
        klab::UInt32 n = 8;
        klab::DoubleReal precision = 1e-12;

        arma::Mat<std::complex<klab::DoubleReal> > f(m, n);
        arma::Mat<std::complex<klab::DoubleReal> > F;

        f(0, 0) = std::complex<klab::DoubleReal>(150.0, 0.0);
        f(0, 1) = std::complex<klab::DoubleReal>(111.0, 0.0);
        f(0, 2) = std::complex<klab::DoubleReal>(129.0, 0.0);
        f(0, 3) = std::complex<klab::DoubleReal>(139.0, 0.0);
        f(0, 4) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(0, 5) = std::complex<klab::DoubleReal>(143.0, 0.0);
        f(0, 6) = std::complex<klab::DoubleReal>(160.0, 0.0);
        f(0, 7) = std::complex<klab::DoubleReal>(90.0, 0.0);
        f(1, 0) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(1, 1) = std::complex<klab::DoubleReal>(108.0, 0.0);
        f(1, 2) = std::complex<klab::DoubleReal>(120.0, 0.0);
        f(1, 3) = std::complex<klab::DoubleReal>(167.0, 0.0);
        f(1, 4) = std::complex<klab::DoubleReal>(182.0, 0.0);
        f(1, 5) = std::complex<klab::DoubleReal>(136.0, 0.0);
        f(1, 6) = std::complex<klab::DoubleReal>(118.0, 0.0);
        f(1, 7) = std::complex<klab::DoubleReal>(85.0, 0.0);
        f(2, 0) = std::complex<klab::DoubleReal>(124.0, 0.0);
        f(2, 1) = std::complex<klab::DoubleReal>(114.0, 0.0);
        f(2, 2) = std::complex<klab::DoubleReal>(132.0, 0.0);
        f(2, 3) = std::complex<klab::DoubleReal>(141.0, 0.0);
        f(2, 4) = std::complex<klab::DoubleReal>(171.0, 0.0);
        f(2, 5) = std::complex<klab::DoubleReal>(184.0, 0.0);
        f(2, 6) = std::complex<klab::DoubleReal>(110.0, 0.0);
        f(2, 7) = std::complex<klab::DoubleReal>(141.0, 0.0);
        f(3, 0) = std::complex<klab::DoubleReal>(127.0, 0.0);
        f(3, 1) = std::complex<klab::DoubleReal>(115.0, 0.0);
        f(3, 2) = std::complex<klab::DoubleReal>(107.0, 0.0);
        f(3, 3) = std::complex<klab::DoubleReal>(85.0, 0.0);
        f(3, 4) = std::complex<klab::DoubleReal>(161.0, 0.0);
        f(3, 5) = std::complex<klab::DoubleReal>(122.0, 0.0);
        f(3, 6) = std::complex<klab::DoubleReal>(97.0, 0.0);
        f(3, 7) = std::complex<klab::DoubleReal>(157.0, 0.0);
        f(4, 0) = std::complex<klab::DoubleReal>(121.0, 0.0);
        f(4, 1) = std::complex<klab::DoubleReal>(107.0, 0.0);
        f(4, 2) = std::complex<klab::DoubleReal>(73.0, 0.0);
        f(4, 3) = std::complex<klab::DoubleReal>(118.0, 0.0);
        f(4, 4) = std::complex<klab::DoubleReal>(153.0, 0.0);
        f(4, 5) = std::complex<klab::DoubleReal>(97.0, 0.0);
        f(4, 6) = std::complex<klab::DoubleReal>(126.0, 0.0);
        f(4, 7) = std::complex<klab::DoubleReal>(170.0, 0.0);
        f(5, 0) = std::complex<klab::DoubleReal>(115.0, 0.0);
        f(5, 1) = std::complex<klab::DoubleReal>(98.0, 0.0);
        f(5, 2) = std::complex<klab::DoubleReal>(90.0, 0.0);
        f(5, 3) = std::complex<klab::DoubleReal>(91.0, 0.0);
        f(5, 4) = std::complex<klab::DoubleReal>(140.0, 0.0);
        f(5, 5) = std::complex<klab::DoubleReal>(87.0, 0.0);
        f(5, 6) = std::complex<klab::DoubleReal>(144.0, 0.0);
        f(5, 7) = std::complex<klab::DoubleReal>(195.0, 0.0);
        f(6, 0) = std::complex<klab::DoubleReal>(120.0, 0.0);
        f(6, 1) = std::complex<klab::DoubleReal>(81.0, 0.0);
        f(6, 2) = std::complex<klab::DoubleReal>(84.0, 0.0);
        f(6, 3) = std::complex<klab::DoubleReal>(77.0, 0.0);
        f(6, 4) = std::complex<klab::DoubleReal>(151.0, 0.0);
        f(6, 5) = std::complex<klab::DoubleReal>(138.0, 0.0);
        f(6, 6) = std::complex<klab::DoubleReal>(144.0, 0.0);
        f(6, 7) = std::complex<klab::DoubleReal>(135.0, 0.0);
        f(7, 0) = std::complex<klab::DoubleReal>(126.0, 0.0);
        f(7, 1) = std::complex<klab::DoubleReal>(71.0, 0.0);
        f(7, 2) = std::complex<klab::DoubleReal>(70.0, 0.0);
        f(7, 3) = std::complex<klab::DoubleReal>(108.0, 0.0);
        f(7, 4) = std::complex<klab::DoubleReal>(149.0, 0.0);
        f(7, 5) = std::complex<klab::DoubleReal>(175.0, 0.0);
        f(7, 6) = std::complex<klab::DoubleReal>(125.0, 0.0);
        f(7, 7) = std::complex<klab::DoubleReal>(93.0, 0.0);

        f.reshape(64, 1);
        arma::Col<std::complex<klab::DoubleReal> > x(f);
        arma::Col<std::complex<klab::DoubleReal> > y;
        arma::Col<std::complex<klab::DoubleReal> > xr;
        f.reshape(8, 8);

        klab::TFWHT2D<std::complex<klab::DoubleReal> > fwht;
        arma::Mat<std::complex<klab::DoubleReal> > fwhtMatrix(64, 64);
        for(klab::UInt32 j=0; j<64; ++j)
        {
            arma::Col<std::complex<klab::DoubleReal> > base(64);
            arma::Col<std::complex<klab::DoubleReal> > col;
            base.fill(0.0);
            base[j] = 1.0;

            arma::Mat<std::complex<klab::DoubleReal> > inMatrix(base);
            arma::Mat<std::complex<klab::DoubleReal> > outMatrix;
            inMatrix.reshape(8, 8);
            fwht.transform(inMatrix, outMatrix);
            outMatrix.reshape(64, 1);

            for(klab::UInt32 i=0; i<outMatrix.n_rows; ++i)
                fwhtMatrix(i, j) = outMatrix[i];
        } 


        TWalshHadamard2DOperator<std::complex<klab::DoubleReal> > op(8, 8);
        op.apply(x, y);

        arma::Col<std::complex<klab::DoubleReal> > xt = arma::trans(fwhtMatrix) * y;

        op.applyAdjoint(y, xr);
        F = xr;
        F.reshape(8, 8);

        TEST_ASSERT(xr.n_rows==64)
        TEST_ASSERT(F.n_rows==8 && F.n_cols==8)
        TEST_ASSERT(klab::Equals(xr[0], xt[0], precision))
        TEST_ASSERT(klab::Equals(xr[1], xt[1], precision))
        TEST_ASSERT(klab::Equals(xr[2], xt[2], precision))
        TEST_ASSERT(klab::Equals(xr[3], xt[3], precision))
        TEST_ASSERT(klab::Equals(xr[4], xt[4], precision))
        TEST_ASSERT(klab::Equals(xr[5], xt[5], precision))
        TEST_ASSERT(klab::Equals(xr[6], xt[6], precision))
        TEST_ASSERT(klab::Equals(xr[7], xt[7], precision))
        TEST_ASSERT(klab::Equals(xr[8], xt[8], precision))
        TEST_ASSERT(klab::Equals(xr[9], xt[9], precision))
        TEST_ASSERT(klab::Equals(xr[10], xt[10], precision))
        TEST_ASSERT(klab::Equals(xr[11], xt[11], precision))
        TEST_ASSERT(klab::Equals(xr[12], xt[12], precision))
        TEST_ASSERT(klab::Equals(xr[13], xt[13], precision))
        TEST_ASSERT(klab::Equals(xr[14], xt[14], precision))
        TEST_ASSERT(klab::Equals(xr[15], xt[15], precision))
        TEST_ASSERT(klab::Equals(xr[16], xt[16], precision))
        TEST_ASSERT(klab::Equals(xr[17], xt[17], precision))
        TEST_ASSERT(klab::Equals(xr[18], xt[18], precision))
        TEST_ASSERT(klab::Equals(xr[19], xt[19], precision))
        TEST_ASSERT(klab::Equals(xr[20], xt[20], precision))
        TEST_ASSERT(klab::Equals(xr[21], xt[21], precision))
        TEST_ASSERT(klab::Equals(xr[22], xt[22], precision))
        TEST_ASSERT(klab::Equals(xr[23], xt[23], precision))
        TEST_ASSERT(klab::Equals(xr[24], xt[24], precision))
        TEST_ASSERT(klab::Equals(xr[25], xt[25], precision))
        TEST_ASSERT(klab::Equals(xr[26], xt[26], precision))
        TEST_ASSERT(klab::Equals(xr[27], xt[27], precision))
        TEST_ASSERT(klab::Equals(xr[28], xt[28], precision))
        TEST_ASSERT(klab::Equals(xr[29], xt[29], precision))
        TEST_ASSERT(klab::Equals(xr[30], xt[30], precision))
        TEST_ASSERT(klab::Equals(xr[31], xt[31], precision))
        TEST_ASSERT(klab::Equals(xr[32], xt[32], precision))
        TEST_ASSERT(klab::Equals(xr[33], xt[33], precision))
        TEST_ASSERT(klab::Equals(xr[34], xt[34], precision))
        TEST_ASSERT(klab::Equals(xr[35], xt[35], precision))
        TEST_ASSERT(klab::Equals(xr[36], xt[36], precision))
        TEST_ASSERT(klab::Equals(xr[37], xt[37], precision))
        TEST_ASSERT(klab::Equals(xr[38], xt[38], precision))
        TEST_ASSERT(klab::Equals(xr[39], xt[39], precision))
        TEST_ASSERT(klab::Equals(xr[40], xt[40], precision))
        TEST_ASSERT(klab::Equals(xr[41], xt[41], precision))
        TEST_ASSERT(klab::Equals(xr[42], xt[42], precision))
        TEST_ASSERT(klab::Equals(xr[43], xt[43], precision))
        TEST_ASSERT(klab::Equals(xr[44], xt[44], precision))
        TEST_ASSERT(klab::Equals(xr[45], xt[45], precision))
        TEST_ASSERT(klab::Equals(xr[46], xt[46], precision))
        TEST_ASSERT(klab::Equals(xr[47], xt[47], precision))
        TEST_ASSERT(klab::Equals(xr[48], xt[48], precision))
        TEST_ASSERT(klab::Equals(xr[49], xt[49], precision))
        TEST_ASSERT(klab::Equals(xr[50], xt[50], precision))
        TEST_ASSERT(klab::Equals(xr[51], xt[51], precision))
        TEST_ASSERT(klab::Equals(xr[52], xt[52], precision))
        TEST_ASSERT(klab::Equals(xr[53], xt[53], precision))
        TEST_ASSERT(klab::Equals(xr[54], xt[54], precision))
        TEST_ASSERT(klab::Equals(xr[55], xt[55], precision))
        TEST_ASSERT(klab::Equals(xr[56], xt[56], precision))
        TEST_ASSERT(klab::Equals(xr[57], xt[57], precision))
        TEST_ASSERT(klab::Equals(xr[58], xt[58], precision))
        TEST_ASSERT(klab::Equals(xr[59], xt[59], precision))
        TEST_ASSERT(klab::Equals(xr[60], xt[60], precision))
        TEST_ASSERT(klab::Equals(xr[61], xt[61], precision))
        TEST_ASSERT(klab::Equals(xr[62], xt[62], precision))
        TEST_ASSERT(klab::Equals(xr[63], xt[63], precision))


        // Limit cases.
        try                                         {op.applyAdjoint(arma::Col<std::complex<klab::DoubleReal> >(), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<std::complex<klab::DoubleReal> >(56), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<std::complex<klab::DoubleReal> >(49), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}

        try                                         {op.applyAdjoint(arma::Col<std::complex<klab::DoubleReal> >(63), xr); TEST_ASSERT(false)}
        catch(KIncompatibleSizeOperatorException&)  {TEST_ASSERT(true)}
        catch(...)                                  {TEST_ASSERT(false)}
    }
    catch(...)
    {
        TEST_ASSERT(false)
    }
}

// ---------------------------------------------------------------------------------------------------- //
