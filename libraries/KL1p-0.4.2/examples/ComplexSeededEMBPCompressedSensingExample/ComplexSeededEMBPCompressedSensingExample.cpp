// KL1p - A portable C++ compressed sensing library.
// Copyright (c) 2011-2012 René Gebel
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

#include "ComplexSeededEMBPCompressedSensingExample.h"

using namespace kl1p;




// ---------------------------------------------------------------------------------------------------- //

void	kl1p::CreateGaussianSignal(klab::UInt32 size, klab::UInt32 sparsity, klab::DoubleReal mean, klab::DoubleReal sigma, arma::Col<klab::DoubleReal>& out)
{
	out.set_size(size);
	out.fill(0.0);

	std::vector<klab::TArrayElement<klab::DoubleReal> > indices;
    for(klab::UInt32 i=0; i<size; ++i)
        indices.push_back(klab::TArrayElement<klab::DoubleReal>(i, klab::KRandom::Instance().generateDoubleReal(0.0, 1.0)));

    std::partial_sort(indices.begin(), indices.begin()+klab::Min(size, sparsity), indices.end(), std::greater<klab::TArrayElement<klab::DoubleReal> >());  

	for(klab::UInt32 i=0; i<sparsity; ++i)
	{
		klab::DoubleReal u1 = klab::KRandom::Instance().generateDoubleReal(0.0, 1.0);
		klab::DoubleReal u2 = klab::KRandom::Instance().generateDoubleReal(0.0, 1.0);

		klab::DoubleReal sign = klab::KRandom::Instance().generateBool() ? -1.0 : 1.0;
		out[indices[i].i()] = sign * ((klab::Sqrt(-2.0*klab::Log(u1)) * klab::Cos(2.0*klab::PI*u2))*sigma + mean);
	}
}

// ---------------------------------------------------------------------------------------------------- //

void	kl1p::WriteToCSVFile(const arma::Col<klab::DoubleReal>& signal, const std::string& filePath)
{
	std::ofstream of(filePath.c_str());
	if(of.is_open())
	{
		for(klab::UInt32 i=0; i<signal.n_rows; ++i)
			of<<i<<";"<<signal[i]<<std::endl;

		of.close();
	}
	else
	{
		std::cout<<"ERROR! Unable to open file \""<<filePath<<"\" !"<<std::endl;
	}
}

// ---------------------------------------------------------------------------------------------------- //

void	kl1p::RunExample()
{
	try
	{
		std::cout<<"Start of KL1p complex seeded-EMBP compressed-sensing example."<<std::endl;
		std::cout<<"Try to determine a sparse complex vector x "<<std::endl;
		std::cout<<"from an underdetermined set of linear measurements y=A*x, "<<std::endl;
		std::cout<<"where A is a complex random block sensing matrix."<<std::endl;

		klab::UInt32 n = 256;					// Size of the original signal x0.
		klab::DoubleReal alpha = 0.20;			// Ratio of the cs-measurements.
		klab::DoubleReal rho = 0.1;				// Ratio of the sparsity of the signal x0.
		klab::UInt32 m = klab::UInt32(alpha*n);	// Number of cs-measurements.
		klab::UInt32 k = klab::UInt32(rho*n);	// Sparsity of the signal x0 (number of non-zero elements).
		klab::UInt64 seed = 0;					// Seed used for random number generation (0 if regenerate random numbers on each launch).
		klab::UInt32 blocks = 8;				// Count number of blocks of sensing matrix.
		klab::UInt32 nb = n / blocks;			// Si ze of the blocks of the sensing matrix.
		klab::DoubleReal alpha0 = 0.5;			// Ratio of the measurements in the first block of the sensing matrix.
		klab::DoubleReal alphaB = 0.0;			// Ratio of the measurements in the others blocks of the sensing matrix.
		klab::DoubleReal dvariance = 1.0;		// Variance of the main diagonal of the sensing matrix.		
		klab::DoubleReal lvariance = 20.0;		// Variance of the lower diagonal of the sensing matrix.
		klab::DoubleReal uvariance = 0.1;		// Variance of the upper diagonal of the sensing matrix.
		bool bWrite = false;					// Write signals to files ?		

		// Initialize random seed if needed.
		if(seed > 0)
			klab::KRandom::Instance().setSeed(seed);

		if(blocks > 1)
			alphaB = klab::DoubleReal(m-klab::UInt32(alpha0*nb)) / klab::DoubleReal((blocks-1)*nb);

		// Display signal informations.
		std::cout<<"=============================="<<std::endl;
		std::cout<<"N="<<n<<" (signal size)"<<std::endl;
		std::cout<<"M="<<m<<"="<<std::setprecision(5)<<(alpha*100.0)<<"% (number of measurements)"<<std::endl;
		std::cout<<"K="<<k<<"="<<std::setprecision(5)<<(rho*100.0)<<"% (signal sparsity)"<<std::endl;
		std::cout<<"B="<<blocks<<" (number of blocks)"<<std::endl;
		std::cout<<"Random Seed="<<klab::KRandom::Instance().seed()<<std::endl;				
				
		arma::Col<klab::DoubleReal> x0;					// Original signal x0 of size n.
		kl1p::CreateGaussianSignal(n, k, 0.0, 1.0, x0);	// Create randomly the original signal x0.

		arma::Col<std::complex<klab::DoubleReal> > cx0(x0.n_rows);	// Complex version of the original signal x0 of size n.
		for(klab::UInt32 i=0; i<x0.n_rows; ++i)
			cx0[i] = std::complex<klab::DoubleReal>(x0[i], 0.0);

		if(bWrite)
			kl1p::WriteToCSVFile(x0, "OriginalSignal.csv");	// Write x0 to a file.

		// Create random gaussian i.i.d matrix A of size (m,n).
		klab::TSmartPointer<kl1p::TOperator<std::complex<klab::DoubleReal> > > A0 = new kl1p::TNormalRandomMatrixOperator<std::complex<klab::DoubleReal> >(m, n, std::complex<klab::DoubleReal>(0.0, 0.0), std::complex<klab::DoubleReal>(1.0, 1.0));
		A0 = new kl1p::TScalingOperator<std::complex<klab::DoubleReal> >(A0, 1.0/klab::Sqrt(klab::DoubleReal(A0->m())));	// Pseudo-normalization of the matrix (required for homogeneous EMBP solver).		

		// Create random seeding matrix A of size(m, n).
		klab::TSmartPointer<kl1p::TBlockOperator<std::complex<klab::DoubleReal> > > A = new kl1p::TSeedingOperator<std::complex<klab::DoubleReal> >(n, blocks, alpha0, alphaB, std::complex<klab::DoubleReal>(dvariance, dvariance), std::complex<klab::DoubleReal>(lvariance, lvariance), std::complex<klab::DoubleReal>(uvariance, uvariance));

		// Use seeding matrix as standard matrix for homogeneous version of EMBP solver.
		klab::TSmartPointer<kl1p::TOperator<std::complex<klab::DoubleReal> > > Ah = klab::static_pointer_cast<kl1p::TOperator<std::complex<klab::DoubleReal> > >(A);
		Ah = new kl1p::TScalingOperator<std::complex<klab::DoubleReal> >(Ah, 1.0/klab::Sqrt(klab::DoubleReal(Ah->m())));	// Pseudo-normalization of the matrix (required for homogeneous EMBP solver).		
		
		// Create random generic-seeding matrix A of size(m, n).
		// A generic-seeding matrix is a block-matrix as a seeding matrix,
		// but the random blocks are replaced with randomly downsampled and permuted blocks 
		// created from a choosen operator (a Fourier operator in the example below).
		klab::TSmartPointer<TOperator<std::complex<klab::DoubleReal> > > op = new TScalingOperator<std::complex<klab::DoubleReal> >(new TFourier1DOperator<std::complex<klab::DoubleReal> >(nb), klab::Sqrt(1.0/n));
		klab::TSmartPointer<kl1p::TBlockOperator<std::complex<klab::DoubleReal> > > Ag = new kl1p::TGenericSeedingOperator<std::complex<klab::DoubleReal> >(op, blocks, klab::UInt32(alpha0*nb), klab::UInt32(alphaB*nb), std::complex<klab::DoubleReal>(dvariance, dvariance), std::complex<klab::DoubleReal>(lvariance, lvariance), std::complex<klab::DoubleReal>(uvariance, uvariance));

		// Use generic-seeding matrix as standard matrix for homogeneous version of EMBP solver.
		klab::TSmartPointer<kl1p::TOperator<std::complex<klab::DoubleReal> > > Agh = klab::static_pointer_cast<kl1p::TOperator<std::complex<klab::DoubleReal> > >(Ag);
		Agh = new kl1p::TScalingOperator<std::complex<klab::DoubleReal> >(Agh, 1.0/klab::Sqrt(klab::DoubleReal(Agh->m())));	// Pseudo-normalization of the matrix (required for homogeneous EMBP solver).		


		std::cout<<"eM1="<<A->m()<<"="<<std::setprecision(5)<<((klab::DoubleReal(A->m())/klab::DoubleReal(n))*100.0)<<"% (effective number of measurements for seeded matrix)"<<std::endl;
		std::cout<<"eM2="<<Ag->m()<<"="<<std::setprecision(5)<<((klab::DoubleReal(Ag->m())/klab::DoubleReal(n))*100.0)<<"% (effective number of measurements for generic-seeded matrix)"<<std::endl;
		std::cout<<"=============================="<<std::endl;

		// Perform cs-measurements of size m.
		arma::Col<std::complex<klab::DoubleReal> > y;
		arma::Col<std::complex<klab::DoubleReal> > yg;
		arma::Col<std::complex<klab::DoubleReal> > yh;
		A0->apply(cx0, yh);
		A->apply(cx0, y);
		Ag->apply(cx0, yg);
		
		klab::DoubleReal tolerance = 1e-3;	// Tolerance of the solution.
		arma::Col<std::complex<klab::DoubleReal> > cx;
		arma::Col<klab::DoubleReal> x;		// Will contain the solution of the reconstruction.

		klab::KTimer timer;

		
		// Compute homogeneous version of EMBP on standard matrix.
		std::cout<<"------------------------------"<<std::endl;
		std::cout<<"[EMBP-h0] Start."<<std::endl;
		timer.start();
		kl1p::TEMBPSolver<klab::DoubleReal, std::complex<klab::DoubleReal> > embph0(tolerance, 500);
		embph0.enableHomogeneous(true);
		embph0.solve(yh, A0, k, cx);
		timer.stop();
		x.set_size(cx.n_rows);
		for(klab::UInt32 i=0; i<x.n_rows; ++i)
			x[i] = cx[i].real();
		std::cout<<"[EMBP-h0] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
				  <<"Iterations="<<embph0.iterations()<<std::endl;
		if(bWrite)
			kl1p::WriteToCSVFile(x, "EMBPh0-Signal.csv");	// Write solution to a file.		

		// Compute homogeneous version of EMBP on seeding matrix.
		std::cout<<"------------------------------"<<std::endl;
		std::cout<<"[EMBP-h1] Start."<<std::endl;
		timer.start();
		kl1p::TEMBPSolver<klab::DoubleReal, std::complex<klab::DoubleReal> > embph1(tolerance, 500);
		embph1.enableHomogeneous(true);
		embph1.solve(y, Ah, k, cx);
		timer.stop();
		x.set_size(cx.n_rows);
		for(klab::UInt32 i=0; i<x.n_rows; ++i)
			x[i] = cx[i].real();
		std::cout<<"[EMBP-h1] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
				  <<"Iterations="<<embph1.iterations()<<std::endl;
		if(bWrite)
			kl1p::WriteToCSVFile(x, "EMBPh1-Signal.csv");	// Write solution to a file.

		// Compute seeded version of EMBP on seeding matrix.
		std::cout<<"------------------------------"<<std::endl;
		std::cout<<"[EMBP-s] Start."<<std::endl;
		timer.start();
		kl1p::TEMBPSolver<klab::DoubleReal, std::complex<klab::DoubleReal> > embps(tolerance);
		embps.solve(y, A, k, cx);
		timer.stop();
		x.set_size(cx.n_rows);
		for(klab::UInt32 i=0; i<x.n_rows; ++i)
			x[i] = cx[i].real();
		std::cout<<"[EMBP-s] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
				  <<"Iterations="<<embps.iterations()<<std::endl;
		if(bWrite)
			kl1p::WriteToCSVFile(x, "EMBPs-Signal.csv");	// Write solution to a file.

		// Compute homogeneous version of EMBP on generic-seeding matrix.
		std::cout<<"------------------------------"<<std::endl;
		std::cout<<"[EMBP-h2] Start."<<std::endl;
		timer.start();
		kl1p::TEMBPSolver<klab::DoubleReal, std::complex<klab::DoubleReal> > embph2(tolerance, 500);
		embph2.enableHomogeneous(true);
		embph2.solve(yg, Agh, k, cx);
		timer.stop();
		x.set_size(cx.n_rows);
		for(klab::UInt32 i=0; i<x.n_rows; ++i)
			x[i] = cx[i].real();
		std::cout<<"[EMBP-h2] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
				  <<"Iterations="<<embph2.iterations()<<std::endl;
		if(bWrite)
			kl1p::WriteToCSVFile(x, "EMBPh2-Signal.csv");	// Write solution to a file.

		// Compute seeded version of EMBP on generic-seeding matrix.
		std::cout<<"------------------------------"<<std::endl;
		std::cout<<"[EMBP-gs] Start."<<std::endl;
		timer.start();
		kl1p::TEMBPSolver<klab::DoubleReal, std::complex<klab::DoubleReal> > embpgs(tolerance);
		embpgs.solve(yg, Ag, k, cx);
		timer.stop();
		x.set_size(cx.n_rows);
		for(klab::UInt32 i=0; i<x.n_rows; ++i)
			x[i] = cx[i].real();
		std::cout<<"[EMBP-gs] Done - SNR="<<std::setprecision(5)<<klab::SNR(x, x0)<<" - "
			      <<"Time="<<klab::UInt32(timer.durationInMilliseconds())<<"ms"<<" - "
				  <<"Iterations="<<embpgs.iterations()<<std::endl;
		if(bWrite)
			kl1p::WriteToCSVFile(x, "EMBPgs-Signal.csv");	// Write solution to a file.


		std::cout<<"------------------------------"<<std::endl;

		std::cout<<std::endl;
		std::cout<<"End of example."<<std::endl;
	}
	catch(klab::KException& e)
	{
		std::cout<<"ERROR! KLab exception : "<<klab::FormatExceptionToString(e)<<std::endl;
	}
	catch(std::exception& e)
	{
		std::cout<<"ERROR! Standard exception : "<<klab::FormatExceptionToString(e)<<std::endl;
	}
	catch(...)
	{
		std::cout<<"ERROR! Unknown exception !"<<std::endl;
	}
}

// ---------------------------------------------------------------------------------------------------- //
