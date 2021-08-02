// KLab - A portable C++ collection of classes for general purpose.
// Copyright (c) 2011-2012 Ren� Gebel
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

#ifndef KLAB_DATETIMEUNITTEST_H
#define KLAB_DATETIMEUNITTEST_H

#include <cpptest.h>




namespace klab
{

// ---------------------------------------------------------------------------------------------------- //

class KDateTimeUnitTest : public Test::Suite
{
public:
    
    KDateTimeUnitTest();
    virtual ~KDateTimeUnitTest();
    
    
protected:

    // Construction/Affectation tests.
    void    testDefaultConstructor();   
    void    testLocalConstructor();   
    void    testTimestampConstructor();   
    void    testDateConstructor();  
    void    testCopyConstructor();
    void    testAffectationOperator();

    // Operators tests.
    void    testEqualOperator();
    void    testNotEqualOperator();
    void    testLessOperator();
    void    testLessOrEqualOperator();
    void    testGreaterOrEqualOperator();
    void    testGreaterOperator();
    void    testSubtractionOperator();

    // Methods tests.
    void    testClear();
    void    testNow();
    void    testYear();
    void    testMonth();
    void    testDay();
    void    testWeekday();
    void    testHour();
    void    testMinute();
    void    testSecond();
    void    testMilliseconds();
    void    testMicroseconds();
    void    testIsValid();
    void    testEquals();
    void    testCompare();
    void    testToString();


private:

    KDateTimeUnitTest(const KDateTimeUnitTest& test);
    KDateTimeUnitTest&  operator=(const KDateTimeUnitTest& test);
};

// ---------------------------------------------------------------------------------------------------- //

}

#endif
