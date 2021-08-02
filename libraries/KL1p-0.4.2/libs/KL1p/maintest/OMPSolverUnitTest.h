// KL1p - A portable C++ compressed sensing library.
// Copyright (c) 2011-2012 Ren� Gebel
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

#ifndef KL1P_OMPSOLVERUNITTEST_H
#define KL1P_OMPSOLVERUNITTEST_H

#include <cpptest.h>




namespace kl1p
{

// ---------------------------------------------------------------------------------------------------- //

class KOMPSolverUnitTest : public Test::Suite
{
public:
    
    KOMPSolverUnitTest();
    virtual ~KOMPSolverUnitTest();
    
    
protected:

    // Constructor/Affectation tests.
    void    testDefaultConstructor();
    void    testToleranceConstructor(); 
    void    testToleranceIterationLimitConstructor(); 
    void    testCopyConstructor();
    void    testAffectationOperator();

    // Methods tests.
    void    testSetTolerance();
    void    testSetIterationLimit();
    void    testTolerance();
    void    testIterationLimit();
    void    testResidualNorm();
    void    testRelativeResidualNorm();
    void    testSparsity();
    void    testIterations();
    void    testLeastSquareSolver();
    void    testEnableHistory();
    void    testIsHistoryEnabled();
    void    testHistory();
    void    testSolve();
    void    testSolve_02();


private:

    KOMPSolverUnitTest(const KOMPSolverUnitTest& test);
    KOMPSolverUnitTest&     operator=(const KOMPSolverUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}

#endif
