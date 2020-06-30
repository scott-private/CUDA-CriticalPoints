/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : CriticalPoints.cpp
============================================================================
*/
#include <cassert>
#include <iostream>
#include <vector>
#include <float.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include "CriticalPoints.h"
#include "KernelCriticalPoints.h"


GpuCPFinder::GpuCPFinder(const char* filename) : BaseCPFinder(filename)
{
    initialization();
}

GpuCPFinder::~GpuCPFinder()
{
    checkCudaErrors(cudaFree(d_tets));
    checkCudaErrors(cudaFree(d_outtets));
    checkCudaErrors(cudaFree(d_outpoints));
    checkCudaErrors(cudaFree(d_outtets2));
    checkCudaErrors(cudaFree(d_outpoints2));
    checkCudaErrors(cudaFree(d_tetsize));
    checkCudaErrors(cudaFree(d_pointsize));
    checkCudaErrors(cudaFree(d_queryvec));
}


void GpuCPFinder::initialization()
{
    // Exactinit
    const int PredParaNum = 14;	        // the size of hostConst
    REAL hostConst[14];					// constants used in exact CUDA CCW test

	//exact computation
	REAL *d_constants;
    checkCudaErrors(cudaMalloc((void **)&d_constants, PredParaNum * sizeof(REAL)) ); 
    cudaExactinit<<<1,1>>>(d_constants); 
	checkCudaErrors(cudaMemcpy(hostConst, d_constants, PredParaNum * sizeof(REAL), cudaMemcpyDeviceToHost) ); 
	
    checkCudaErrors(cudaMemcpyToSymbol(deviceConst, hostConst, PredParaNum * sizeof(REAL)) ); 
    checkCudaErrors(cudaFree(d_constants));

    // Tet init
    tetnum = seqtets.size() / 4;

    checkCudaErrors(cudaMalloc((void **)&d_tets, seqtets.size() * sizeof(PointVelocity3D)) ); 
    checkCudaErrors(cudaMemcpy(d_tets, seqtets.data(), seqtets.size() * sizeof(PointVelocity3D), cudaMemcpyHostToDevice) ); 

    checkCudaErrors(cudaMalloc((void **)&d_outtets, tetnum * sizeof(int)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_outpoints, tetnum * 4 * sizeof(PointTet)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_outtets2, tetnum * sizeof(int)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_outpoints2, tetnum * 4 * sizeof(PointTet)) ); 

    checkCudaErrors(cudaMalloc((void **)&d_tetsize, sizeof(int)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_pointsize, sizeof(int)) );

    checkCudaErrors(cudaMemset(d_tetsize,0,sizeof(int))) ; 
    checkCudaErrors(cudaMemset(d_pointsize,0,sizeof(int))) ;
    
    cudaInit<<<100,256>>>(d_outtets, d_outpoints, tetnum);

    Vector3D query = {0.0,0.0,0.0};
    checkCudaErrors(cudaMalloc((void **)&d_queryvec, sizeof(Vector3D)) );
    checkCudaErrors(cudaMemcpy(d_queryvec, &query, sizeof(Vector3D), cudaMemcpyHostToDevice) ); 
}


void GpuCPFinder::search()
{
    cudaSearch<<<100,256>>>(d_tets, tetnum, d_queryvec, d_outtets, d_outpoints);
    checkCudaErrors(cudaDeviceSynchronize());

    cudaCompress<<<100,256>>>(d_outtets, d_outpoints, tetnum, d_outtets2, d_outpoints2, d_tetsize, d_pointsize);

    int pointSize = 0, tetSize = 0;
    checkCudaErrors(cudaMemcpy(&tetSize, d_tetsize, sizeof(int), cudaMemcpyDeviceToHost) ); 
    checkCudaErrors(cudaMemcpy(&pointSize, d_pointsize, sizeof(int), cudaMemcpyDeviceToHost) ); 
    // std::cout << pointSize << " " << tetSize << std::endl;
    criticalTets.resize(tetSize);
    criticalPoints.resize(pointSize);
    checkCudaErrors(cudaMemcpy(criticalTets.data(), d_outtets2, tetSize * sizeof(int), cudaMemcpyDeviceToHost) ); 
    checkCudaErrors(cudaMemcpy(criticalPoints.data(), d_outpoints2, pointSize * sizeof(PointTet), cudaMemcpyDeviceToHost) ); 
}