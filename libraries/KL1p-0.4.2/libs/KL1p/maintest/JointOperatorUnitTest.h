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

#ifndef KL1P_JOINTOPERATORUNITTEST_H
#define KL1P_JOINTOPERATORUNITTEST_H

#include <cpptest.h>
#include "../include/Operator.h"




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

class KJointOperatorUnitTest : public Test::Suite
{
public:
    
    KJointOperatorUnitTest();
    virtual ~KJointOperatorUnitTest();


public:

    static void     GenerateTestOperators_01(std::vector<klab::TSmartPointer<kl1p::TOperator<klab::DoubleReal> > >& operators, klab::UInt32& m, klab::UInt32& n, klab::UInt32& blockRows, klab::UInt32& blockCols);
    static void     GenerateTestOperators_02(std::vector<klab::TSmartPointer<kl1p::TOperator<klab::DoubleReal> > >& operators, klab::UInt32& m, klab::UInt32& n, klab::UInt32& blockRows, klab::UInt32& blockCols);
    
    
protected:

    // Constructor/Affectation tests.
    void    testOpConstructor(); 
    void    testOpArrayConstructor(); 
    void    testCopyConstructor();

    // Methods tests.
    void    testCountBlockRows();
    void    testCountBlockColumns();
    void    testBlock();
    void    testIsZeroBlock();
    void    testInBlock();
    void    testApply();
    void    testApplyAdjoint();
    void    testColumn();
    void    testColumnAdjoint();
    void    testApplyBlockVariance();
    void    testApplyBlockVarianceAdjoint();


private:

    KJointOperatorUnitTest(const KJointOperatorUnitTest& test);
    KJointOperatorUnitTest& operator=(const KJointOperatorUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}


#endif
