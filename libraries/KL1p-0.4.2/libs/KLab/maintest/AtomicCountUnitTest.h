// KLab - A portable C++ collection of classes for general purpose.
// Copyright (c) 2011-2012 Ren? Gebel
// 
// This file is part of the KLab C++ library.
// This library is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY of fitness for any purpose. 
//
// This library is free software; You can redistribute it and/or modify it 
// under the terms of the GNU Lesser General Public License (LGPL) 
// as published by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// See http://www.opensource.org/licenses for more info.

#ifndef KLAB_ATOMICCOUNTUNITTEST_H
#define KLAB_ATOMICCOUNTUNITTEST_H

#include <cpptest.h>




namespace klab
{

// ---------------------------------------------------------------------------------------------------- //

class KAtomicCountUnitTest : public Test::Suite
{
public:
    
    KAtomicCountUnitTest();
    virtual ~KAtomicCountUnitTest();
    
    
protected:

    // Construction/Affectation tests.
    void    testDefaultConstructor();
    void    testValueConstructor();

    // Operators tests.
    void    testIncrementOperator();
    void    testDecrementOperator();
    void    testCastOperator();


private:

    KAtomicCountUnitTest(const KAtomicCountUnitTest& test);
    KAtomicCountUnitTest&   operator=(const KAtomicCountUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}

#endif
