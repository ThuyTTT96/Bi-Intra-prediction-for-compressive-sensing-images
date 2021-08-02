KL1p v0.4.2 - A portable C++ compressed sensing library.
http://kl1p.sourceforge.net




Table of contents
=================

1: General Information
   1.1: Introduction
   1.2: Features

2: Installation
   2.1: Preliminaries
   2.2: Installation on Linux
   2.3: Installation on Mac OS X
   2.4: Installation on MS Windows
   2.5: Manual installation
   
3: Compiling and Linking
   3.1: Compiling and linking on Linux
   3.2: Compiling and linking on Mac OS X
   3.3: Compiling and linking on MS Windows
   3.4: Examples
   3.5: Test cases
   
4: Documentation

5: Bug Reports

6: Credits

7: License

8: References and Resources
   8.1: References
   8.2: Resources




1: General Information
   ===================
   
1.1: Introduction
     ------------
KL1p is a portable C++ framework for handling sparse recovery of inverse problems
of underdetermined linear systems, like in compressed sensing technique (CS). 
It is distributed under a license that is useful in both open-source 
and proprietary contexts.

The focus of the KL1p library lies on usability and extendability. Several of the 
most common CS algorithms are implemented and problem definitions are made through
combination of operators. These operators played the same role as matrices 
in linear algebra but are implemented with their equivalent efficient functions 
when possible (e.g Fourier matrix implemented with FFT function) 
and new ones can be easily added. This method combines the flexibility of 
the matricial method with the efficiency of the functional method.

Please note that KL1p is still in development and is not recommended for 
a production environment.


1.2: Features
     --------	 
Currently implemented compressed sensing algorithms :
  - Basis Pursuit [1][2]
  - Orthogonal Matching Pursuit (OMP) [3]
  - Regularized Orthogonal Matching Pursuit (ROMP) [4]
  - Compressive Sampling Matching Pursuit (CoSaMP) [5]
  - Subspace Pursuit [6]
  - Smoothed L0 (SL0) [7]
  - Approximate Message Passing (AMP) [8]
  - Expectation Maximization Belief Propagation (EMBP) [9][10]
  
Currently implemented operators :
  - Random matrices (normal and uniform random).
  - Fourier transformation.
  - Discrete Cosine transformation.
  - Walsh-Hadamard transformation.
  - Daubechies Wavelets transformation.  
  - Combination of operators through addition, multiplication, concatenation, ...
  - Seeding matrices [9] (some special block tridiagonal sensing matrices used in conjunction 
    with EMBP algorithm to outperform the Donoho-Tanner phase transition). 
  - Some common matrices used in various CS reconstruction techniques 
    like downsampling or permutation matrices.
  - And many others for various purposes ...
  
For a list of changes please read the file "VERSIONS.txt"


2: Installation 
   ============

2.1: Preliminaries
     -------------
The KL1p library is dependent on the following tools and libraries :     

  * Standard compliant C++ compiler :
      KL1p makes use of template meta-programming, recursive templates 
      and template based function overloading. As such, C++ compilers 
      which do not fully implement the C++ standard may not work correctly.
	  For information, KL1p was successfully built with the following compilers :
	    - Microsoft Visual C++ 2008 Express (32/64 bits compilers on Windows XP)
	    - Microsoft Visual C++ 2010 Express (32/64 bits compilers on Windows Seven)
		- GCC 4.2.1 (64 bits compiler on Mac OS X Snow Leopard)
		- GCC 4.3.2 (64 bits compiler on Linux Debian)
	    - GCC 4.4.3 (32 bits compiler on Linux Ubuntu)	   

  * The cross platform build system CMake v2.8.1 : 
      KL1p use the open source software CMake as build system. The minimum 
	  required version is 2.8.1. 
	  For installation and information about CMake, please refer to 
	  http://www.cmake.org
	   
  * The C++ linear algebra library Armadillo v3.1.92 [11] : 
      KL1p use the Armadillo C++ linear algebra library for matrix operations.
	  The version 3.1.92 of the library is already provided with the KL1p package 
	  under the folder "libs/Armadillo". 
	  For information about Armadillo or installation of recent version, please refer to 
	  http://arma.sourceforge.net
	   	  
		  
The following modules are facultative for standard installation of KL1p : 

  * The C Fast Fourier Transform library FFTW v3.3.1 :
      By default, KL1p already provides basic Fast Fourier Transform operations. 
	  But for performance increase during the computation of FFTs or derived 
	  operations (DCT, convolution, correlation, ...) the installation of the 
	  FFTW GPL library is recommended. The minimum required version is 3.3.1.
	  For installation and information about FFTW, please refer to 
	  http://www.fftw.org
	  See section 3.3 for example of configuration of KL1p with FFTW.

  * The linear algebra packages BLAS or LAPACK : 
      Linear algebra operations may be sped up with the activation
	  of the BLAS or LAPACK packages through the Armadillo library. 
	  For information about this configuration of the Armadillo library, please refer to 
	  http://arma.sourceforge.net/download.html
	  or see the file "libs/Armadillo/README.txt", section 2 Installation.

  * The C++ test framework CppTest v1.1.0 :
      The KL1p package provides the unit test cases used during the development process (see section 3.4). 
	  These test cases are essentially provided for documentation purpose. But if 
	  for any reason you need to build these tests cases, you must install 
	  the C++ test framework CppTest. The minimum required version is 1.1.0.
	  For installation and information about CppTest, please refer to 
	  http://cpptest.sourceforge.net


2.2: Installation on Linux
     ---------------------
You can use the manual installation process as described in section 2.5, 
or the following automatic installation for GCC compilers.

If CMake is not already be present on your system, download it from http://www.cmake.org
and install it on your system. On major Linux systems (such as Fedora, Ubuntu, Debian, etc), 
cmake is available as a pre-built package, though it may need to be explicitly installed 
(using a tool such as PackageKit, yum, rpm, apt, aptitude).	 

KL1p provides scripts to automatically build GCC versions of the library. These files 
are located in the "build/unix" folder and are named :
  - "UnixMakefile_x86_Make.sh" for 32 bits versions (debug and release).
  - "UnixMakefile_x64_Make.sh" for 64 bits versions (debug and release).
  
Open a shell (command line), change into the directory "build/unix" from the directory 
that was created by unpacking the KL1p archive, and type for example the following command for 
building 32 bits versions of the library :
  ./UnixMakefile_x86_Make.sh
  
For compilers other than GCC, please use the manual installation process as described in section 2.5.

When the build is complete, the binaries will be found into the directory "bin/unix". Depending of the build
chosen, the following files will be present :
  - "libKLab_d.a" for 32 bits debug version.
  - "libKLab.a" for 32 bits release version.
  - "libKLab64_d.a" for 64 bits debug version.
  - "libKLab64.a" for 64 bits release version.  

You must configure your compiler to link to the corresponding KL1p files : 
  - Configure your compiler to use the location for header files (in addition to the locations 
    it uses already), the folder "include" under the directory that was created by unpacking the KL1p archive.
  - Configure your compiler to link with the previously builded libraries under directory "bin/unix". 

  
2.3: Installation on Mac OS X
     ------------------------
You can use the manual installation process as described in section 2.5, 
or the following automatic installation for GCC compilers.

If CMake is not already be present on your system, download it from http://www.cmake.org
and install it on your system. 	 

KL1p provides scripts to automatically build GCC versions of the library. These files 
are located in the "build/macos" folder and are named :
  - "UnixMakefile_x86_Make.sh" for 32 bits versions (debug and release).
  - "UnixMakefile_x64_Make.sh" for 64 bits versions (debug and release).
  
Open a shell (command line), change into the directory "build/macos" from the directory 
that was created by unpacking the KL1p archive, and type for example the following command for 
building 32 bits versions of the library :
  ./UnixMakefile_x86_Make.sh
  
For compilers other than GCC, please use the manual installation process as described in section 2.5.

When the build is complete, the binaries will be found into the directory "bin/macos". Depending of the build
chosen, the following files will be present :
  - "libKLab_d.a" for 32 bits debug version.
  - "libKLab.a" for 32 bits release version.
  - "libKLab64_d.a" for 64 bits debug version.
  - "libKLab64.a" for 64 bits release version.  

You must configure your compiler to link to the corresponding KL1p files : 
  - Configure your compiler to use the location for header files (in addition to the locations 
    it uses already), the folder "include" under the directory that was created by unpacking the KL1p archive.
  - Configure your compiler to link with the previously builded libraries under directory "bin/macos". 
  

2.4: Installation on MS Windows
     --------------------------
You can use the manual installation process as described in section 2.5, 
or the following automatic installation for Visual C++ compilers.
	 
If CMake is not already be present on your system, download it from http://www.cmake.org
and install it on your system.
	 
KL1p provides scripts to automatically build Visual C++ versions of the library. These files 
are located in the "build/win" folder and are named :
  - "Visual2005_x86_Make.bat" for 32 bits versions built on Visual C++ 2005 and Visual C++ 2005 Express (debug and release).
  - "Visual2005_x64_Make.bat" for 64 bits versions built on Visual C++ 2005 and Visual C++ 2005 Express (debug and release).
  - "Visual2008_x86_Make.bat" for 32 bits versions built on Visual C++ 2008 and Visual C++ 2008 Express (debug and release).
  - "Visual2008_x64_Make.bat" for 64 bits versions built on Visual C++ 2008 and Visual C++ 2008 Express (debug and release).
  - "Visual2010_x86_Make.bat" for 32 bits versions built on Visual C++ 2010 (debug and release).
  - "Visual2010_x64_Make.bat" for 64 bits versions built on Visual C++ 2010 (debug and release).
  - "Visual2010Express_x86_Make.bat" for 32 bits versions built on Visual C++ 2010 Express (debug and release).
  - "Visual2010Express_x64_Make.bat" for 64 bits versions built on Visual C++ 2010 Express (debug and release).
  
Open a file explorer, change into the directory "build/win" from the directory 
that was created by unpacking the KL1p archive, and run for example the following script for
building 32 bits versions of the library with Visual C++ 2010 :
  Visual2010_x86_Make.bat	  
  
As a courtesy, the KL1p package already contains pre-generated versions of Visual C++ projects. 
These projects are located in folder "build/win" and are named : 
  - "KLab.vcproj" for Visual C++ 2008 and Visual C++ 2008 Express environments.
  - "KLab.vcxproj" for Visual C++ 2010 and Visual C++ 2010 Express environments.

For compilers other than Visual C++, please use the manual installation process as described in section 2.5.

When the build is complete, the binaries will be found into the directory "bin/win". Depending of the build
chosen, the following files will be present :
  - "KLab_d.lib" for 32 bits debug version.
  - "KLab.lib" for 32 bits release version.
  - "KLab64_d.lib" for 64 bits debug version.
  - "KLab64.lib" for 64 bits release version. 
  
You must configure your compiler to link to the corresponding KL1p files : 
  - Configure your compiler to use the location for header files (in addition to the locations 
    it uses already), the folder "include" under the directory that was created by unpacking the KL1p archive.
  - Configure your compiler to link with the previously builded libraries under directory "bin/win". 


2.5: Manual installation
     -------------------
If CMake is not already be present on your system, download it from http://www.cmake.org
and install it on your system.
	 
Open a shell (command line), change into the directory "libs/KLab/main" from the directory 
that was created by unpacking the KL1p archive. This directory contains the file "CMakeLists.txt"
required by cmake to build the binaries. The "CMakeLists.txt" file may need to be configured for your system.
Please refer to cmake documentation at http://www.cmake.org/cmake/help/documentation.html
to see how to create projects from a "CMakeLists.txt" file. The standard cmake command looks like :
  cmake -DCMAKE_BUILD_TYPE=[BUILD_TYPE] -G[COMPILER_NAME] [SOURCE_PATH]
where :
 - [BUILD_TYPE] corresponds to the type of the build ("Release" or "Debug" for example). 
 - [COMPILER_NAME] corresponds to the name of the desired compiler (see cmake documentation for a complete list of compiler names). 
 - [SOURCE_PATH] corresponds to the location of the "CMakeLists.txt" file.
  
As the project is created by cmake command, build it with the chosen compiler. 
When the build is complete, the binaries will be found into the directory "bin/unix" for Linux systems
or into the directory "bin/macos" for Mac OS X systems or into the directory "bin/win" for Windows systems.

You must configure your compiler to link to the corresponding KL1p files : 
  - Configure your compiler to use the location for header files (in addition to the locations 
    it uses already), the folder "include" under the directory that was created by unpacking the KL1p archive.
  - Configure your compiler to link with the previously builded libraries under directory "bin/unix", "bin/macos" or "bin/win"
    depending on your system. 


3: Compiling and Linking 
   =====================

3.1: Compiling and linking on Linux
     ------------------------------ 
In order to use KL1p library, your program need to include the header file "KL1pInclude.h" 
located in the "include" folder :
  #include <KL1pInclude.h>
  
Depending on your build type, your program also need to link to one of KL1p library previously generated
in the "build/unix" folder (see section 2) :
  - "libKLab_d.a" for 32 bits debug version.
  - "libKLab.a" for 32 bits release version.
  - "libKLab64_d.a" for 64 bits debug version.
  - "libKLab64.a" for 64 bits release version.  
  
See section 3.4 to see examples using KL1p.


3.2: Compiling and linking on Mac OS X
     ---------------------------------
In order to use KL1p library, your program need to include the header file "KL1pInclude.h" 
located in the "include" folder :
  #include <KL1pInclude.h>
  
Depending on your build type, your program also need to link to one of KL1p library previously generated
in the "build/macos" folder (see section 2) :
  - "libKLab_d.a" for 32 bits debug version.
  - "libKLab.a" for 32 bits release version.
  - "libKLab64_d.a" for 64 bits debug version.
  - "libKLab64.a" for 64 bits release version.  
  
See section 3.4 to see examples using KL1p.


3.3: Compiling and linking on MS Windows
     ----------------------------------- 
In order to use KL1p library, your program need to include the header file "KL1pInclude.h" 
located in the "include" folder :
  #include <KL1pInclude.h>
  
Depending on your build type, your program also need to link to one of KL1p library previously generated
in the "build/win" folder (see section 2) :	 
  - "KLab_d.lib" for 32 bits debug version.
  - "KLab.lib" for 32 bits release version.
  - "KLab64_d.lib" for 64 bits debug version.
  - "KLab64.lib" for 64 bits release version. 

See section 3.4 to see examples using KL1p.


3.4: Examples
     -------- 
The "examples" directory contains several quick examples that use the KL1p library :

  * CompressedSensingExample
      A basic reconstruction example. A random sparse signal is created, downsampled and reconstructed with
	  each of the algorithms implemented in KL1p.
  
  * FourierCompressedSensingExample
      A reconstruction example based on Fourier transformation. A random sparse signal is created, 
	  transformed with Fourier transformation, downsampled and reconstructed with each of the algorithms 
	  implemented in KL1p.
	  To speed up the reconstruction process, FFTW may be enabled with one of the following preprocessor 
	  command (before inclusion of "KL1pInclude.h" file) :
	
	    // Use dynamic library version of FFTW.
	    // The FFTW dynamic library must be accessible by the example executable (e.g. must be located
	    // in the same folder or in the system folder) and must be named "libfftw3-3.so" for linux system 
	    // or "libfftw3-3.dll" for windows system.
	    #define KSCI_ENABLE_DYNAMIC_FFTW	
	    #include <KL1pInclude.h>
	
        or
	
	    // Use static library version of FFTW.
	    #define KSCI_ENABLE_STATIC_FFTW
	    #include <KL1pInclude.h>	   
  
  * SeededEMBPCompressedSensingExample
      A reconstruction example with a special sensing matrix called "seeding matrix". The combination of
	  a seeding matrix and the EMBP algorithm may provide better reconstruction performance [9].
	  A random sparse signal is created, measured with the seeding matrix and reconstructed thanks to 
	  the EMBP algorithm.
  
  * ComplexSeededEMBPCompressedSensingExample
      Same as SeededEMBPCompressedSensingExample above but with the recontruction of a complex (imaginary) signal.


3.4: Test cases
     ---------- 
The KL1p package provides the unit test cases used during the development process. 
These test cases are essentially provided for documentation purpose. But if for any reason
you need to build these tests cases, you must install the C++ test framework CppTest. 
The minimum required version is 1.1.0. For installation and information about CppTest, please refer to 
http://cpptest.sourceforge.net/

The unit test cases are located for each sub-libraries under the sub-folder "maintest" : 
  - folder "libs/KLab/maintest" for unit test cases of KLab library.
  - folder "libs/KSci/maintest" for unit test cases of KSci library.
  - folder "libs/KL1p/maintest" for unit test cases of KL1p library.
  
Each "maintest" subfolder contains a file named "CMakeLists.txt" required by cmake to build 
the test executable. The "CMakeLists.txt" file may need to be configured for your system. 
Please refer to cmake documentation at http://www.cmake.org/cmake/help/documentation.html
to see how to create projects from a "CMakeLists.txt" file. The standard cmake command looks like :
  cmake -DCMAKE_BUILD_TYPE=[BUILD_TYPE] -G[COMPILER_NAME] [SOURCE_PATH]
where :
 - [BUILD_TYPE] corresponds to the type of the build ("Release" or "Debug" for example). 
 - [COMPILER_NAME] corresponds to the name of the desired compiler (see cmake documentation for a complete list of compiler names). 
 - [SOURCE_PATH] corresponds to the location of the "CMakeLists.txt" file.
  
As the project is created by cmake command, build it with the chosen compiler. 
When the build is complete, the corresponding test executable will be found into the directory "bin/unix" for Linux systems 
or into the directory "bin/macos" for Mac OS X systems or into the directory "bin/win" for Windows systems.


4: Documentation 
   =============
For information about KL1p, please visit http://kl1p.sourceforge.net


5: Bug Reports
   ===========
If you find a bug in the library, please make a small self-contained program 
which exposes the bug and send the program source as well as the bug description 
to kl1p-contact@lists.sourceforge.net

In the bug description please include:
  - Information about your system and KL1p version.
  - If the bug was a crash, supply the exact text that was printed out 
    when the exception or error occured.
  - Any other relevant information.


6: Credits
   =======
Developer: 
  Rene Gebel


7: License
   =======
This library is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY of fitness for any purpose. 
This library is free software; You can redistribute it and/or modify it 
under the terms of the GNU Lesser General Public License (LGPL) as published 
by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.
See the "LICENSE.txt" file and http://www.opensource.org/licenses for more info.


8: References and Resources
   ========================	
   
8.1: References
     ----------
[1]  Emmanuel Candes and Justin Romberg.
     l1-MAGIC : Recovery of Sparse Signals via Convex Programming.
	 Caltech, 2005
	 
[2]  Justin Romberg.
     l1-MAGIC
     http://users.ece.gatech.edu/~justin/l1magic, 2005
	 
[3]  Joel A. Tropp and Anna C. Gilbert.
     Signal Recovery From Random Measurements Via Orthogonal Matching Pursuit.
	 IEEE Transactions on Information Theory, Vol.53, NO.12, 2007
	 
[4]  D. Needell and R. Vershynin.
     Signal recovery from incomplete and inaccurate measurements via regularized orthogonal matching pursuit.
     IEEE J. Selected Topics in Signal Process., vol.4, pp.310-316, 2010. 	 

[5]  D. Needell and J.A. Tropp.
     CoSaMP: Iterative signal recovery from incomplete and inaccurate samples.
	 Applied and Computational Harmonic Analysis, 2009

[6]  Wei Dai and Olgica Milenkovic.
     Subspace Pursuit for Compressive Sensing Signal Reconstruction.
	 arXiv:0803.0811v3, 2009

[7]  G. Hosein Mohimani, Massoud Babaie-Zadeh and Christian Jutten.
     Smoothed L0 (SL0) Algorithm for Sparse Decomposition
	 http://ee.sharif.edu/~SLzero, 2008

[8]  D. L. Donoho, A. Maleki, and A. Montanari.
     Message-passing algorithms for compressed sensing.
	 Proc. National Academy of the Sciences, vol.106, no.45, pp.18914-18919, 2009

[9]  F. Krzakala, M. Mezard, F. Sausset, Y. F. Sun and L. Zdeborova.
     Statistical-physics-based reconstruction in compressed sensing.
	 arXiv:1109.4424v3, 2012
	 
[10] J. Barbier and F. Krzakala. 
     ASPICS: Applying Statistical Physics to Inference in Compressed Sensing.
     http://aspics.krzakala.org, 2012
	 
[11] Conrad Sanderson.
     Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments.
     Technical Report, NICTA, 2010.
	 
	 
8.2: Resources
     ---------
- Nuit-Blanche, blog of Igor Carron.
  http://nuit-blanche.blogspot.fr
    Daily news and a huge list of resources from compressed sensing and matrix factorization communities. 
	
- Pursuits in the Null Space, blog of Bob L. Sturm.
  http://media.aau.dk/null_space_pursuits 
    Many Matlab codes and comparative studies of various compressed sensing algorithms.
