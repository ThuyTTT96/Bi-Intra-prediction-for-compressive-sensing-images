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

#ifndef KL1P_COMPLEXPROXYOPERATORUNITTEST_H
#define KL1P_COMPLEXPROXYOPERATORUNITTEST_H

#include <cpptest.h>




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

class KComplexProxyOperatorUnitTest : public Test::Suite
{
public:
    
    KComplexProxyOperatorUnitTest();
    virtual ~KComplexProxyOperatorUnitTest();
    
    
protected:

    // Constructor/Affectation tests.
    void    testOpConstructor(); 
    void    testCopyConstructor();

    // Methods tests.
    void    testOp();
    void    testIsZero();
    void    testSum();
    void    testNormFrobenius();
    void    testSquaredNormFrobenius();
    void    testMean();
    void    testVariance();
    void    testApply();
    void    testApplyAdjoint();
    void    testColumn();
    void    testColumnAdjoint();


private:

    KComplexProxyOperatorUnitTest(const KComplexProxyOperatorUnitTest& test);
    KComplexProxyOperatorUnitTest&  operator=(const KComplexProxyOperatorUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}

#endif
