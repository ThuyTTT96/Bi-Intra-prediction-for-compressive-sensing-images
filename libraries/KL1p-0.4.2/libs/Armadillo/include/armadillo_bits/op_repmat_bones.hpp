// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_repmat
//! @{



class op_repmat
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_repmat>& in);
  };



//! @}