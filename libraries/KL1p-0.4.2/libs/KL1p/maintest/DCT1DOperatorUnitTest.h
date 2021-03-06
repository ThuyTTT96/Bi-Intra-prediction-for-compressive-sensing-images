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

#ifndef KL1P_DCT1DOPERATORUNITTEST_H
#define KL1P_DCT1DOPERATORUNITTEST_H

#include <cpptest.h>




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

class KDCT1DOperatorUnitTest : public Test::Suite
{
public:
    
    KDCT1DOperatorUnitTest();
    virtual ~KDCT1DOperatorUnitTest();
    
    
protected:

    // Constructor/Affectation tests.

    void    testNConstructor(); 
    void    testCopyConstructor();

    // Methods tests.
    void    testApply();
    void    testApply_Complex();
    void    testApplyAdjoint();
    void    testApplyAdjoint_Complex();


private:

    KDCT1DOperatorUnitTest(const KDCT1DOperatorUnitTest& test);
    KDCT1DOperatorUnitTest&     operator=(const KDCT1DOperatorUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}

#endif
