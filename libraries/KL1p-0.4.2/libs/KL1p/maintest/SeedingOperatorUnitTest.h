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

#ifndef KL1P_SEEDINGOPERATORUNITTEST_H
#define KL1P_SEEDINGOPERATORUNITTEST_H

#include <cpptest.h>




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

class KSeedingOperatorUnitTest : public Test::Suite
{
public:
    
    KSeedingOperatorUnitTest();
    virtual ~KSeedingOperatorUnitTest();

    
protected:

    // Constructor/Affectation tests.
    void    testConstructor(); 
    void    testConstructor2(); 
    void    testCopyConstructor();

    // Methods tests.
    void    testCountSeededBlocks();
    void    testM0();
    void    testN0();
    void    testMb();
    void    testNb();
    void    testDiagonalVariance();
    void    testLowerVariance();
    void    testUpperVariance();
    void    testCountBlockRows();
    void    testCountBlockColumns();
    void    testBlock();
    void    testIsZeroBlock();
    void    testInBlock();
    void    testApply();
    void    testApplyAdjoint();
    void    testColumn();
    void    testColumnAdjoint();


private:

    KSeedingOperatorUnitTest(const KSeedingOperatorUnitTest& test);
    KSeedingOperatorUnitTest&	operator=(const KSeedingOperatorUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}


#endif
