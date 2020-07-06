/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : GpuCriticalPoints.cu
============================================================================
*/
#include <cassert>
#include <iostream>
#include <vector>
#include <chrono>
#include <float.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <thrust/copy.h>

#include "CriticalPoints.h"
#include "KernelCriticalPoints.h"


GpuCPFinder::GpuCPFinder(const char* filename) : BaseCPFinder(filename)
{
    initialization();
}

GpuCPFinder::~GpuCPFinder()
{
    checkCudaErrors(cudaFree(d_tets1));
    checkCudaErrors(cudaFree(d_outstencil));
    checkCudaErrors(cudaFree(d_outpoints));
    checkCudaErrors(cudaFree(d_outpoints2));
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

    checkCudaErrors(cudaMalloc((void **)&d_tets1, seqtets.size() * sizeof(PointVelocity3D)) ); 
    checkCudaErrors(cudaMemcpy(d_tets1, seqtets.data(), seqtets.size() * sizeof(PointVelocity3D), cudaMemcpyHostToDevice) ); 

    checkCudaErrors(cudaMalloc((void **)&d_outstencil, tetnum * 4 * sizeof(int)) );
    checkCudaErrors(cudaMalloc((void **)&d_outpoints, tetnum * 4 * sizeof(PointTet)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_outpoints2, tetnum * 4 * sizeof(PointTet)) ); 

    checkCudaErrors(cudaMalloc((void **)&d_pointsize, sizeof(int)) );

    checkCudaErrors(cudaMemset(d_pointsize,0,sizeof(int))) ;
    
    cudaInit<<<100,256>>>(d_outstencil, d_outpoints, tetnum);

    Vector3D query = {0.0,0.0,0.0};
    checkCudaErrors(cudaMalloc((void **)&d_queryvec, sizeof(Vector3D)) );
    checkCudaErrors(cudaMemcpy(d_queryvec, &query, sizeof(Vector3D), cudaMemcpyHostToDevice) ); 
}


void GpuCPFinder::search()
{
    cudaSearch<<<100,256>>>(d_tets1, tetnum, d_queryvec, d_outpoints);
    checkCudaErrors(cudaDeviceSynchronize());

    cudaCompress<<<100,256>>>(d_outpoints, tetnum, d_outpoints2, d_pointsize);

    int pointSize = 0;
    checkCudaErrors(cudaMemcpy(&pointSize, d_pointsize, sizeof(int), cudaMemcpyDeviceToHost) ); 
    // std::cout << pointSize << " " << tetSize << std::endl;
    criticalPoints.resize(pointSize);
    checkCudaErrors(cudaMemcpy(criticalPoints.data(), d_outpoints2, pointSize * sizeof(PointTet), cudaMemcpyDeviceToHost) ); 
}


__global__
void cudaRoughJudge(const PointVelocity3D* tets,
                    const int tetSize, 
                    Vector3D* queryvec,
                    int* stencil,
                    PointTet* criticalPoints,
                    int * pointSize)
{
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t<tetSize; t+=dim){
        const PointVelocity3D* pv = &tets[t*4];
        Vector3D point = {NAN,NAN,NAN};

        // 初步判断
        REAL xMax,xMin,yMax,yMin,zMax,zMin;
        xMax = xMin = pv[0].getVelocity().x;
        yMax = yMin = pv[0].getVelocity().y; 
        zMax = zMin = pv[0].getVelocity().z; 
        for(int it=1;it<4;it++){
            xMax = pv[it].getVelocity().x > xMax ? pv[it].getVelocity().x : xMax;
            xMin = pv[it].getVelocity().x < xMin ? pv[it].getVelocity().x : xMin;
            yMax = pv[it].getVelocity().y > yMax ? pv[it].getVelocity().y : yMax;
            yMin = pv[it].getVelocity().y < yMin ? pv[it].getVelocity().y : yMin;
            zMax = pv[it].getVelocity().z > zMax ? pv[it].getVelocity().z : zMax;
            zMin = pv[it].getVelocity().z < zMin ? pv[it].getVelocity().z : zMin;
        } 
        if(xMax < queryvec->x || xMin > queryvec->x){continue;}
        if(yMax < queryvec->y || yMin > queryvec->y){continue;}
        if(zMax < queryvec->z || zMin > queryvec->z){continue;}

        stencil[t*4] = t+1; 
        stencil[t*4+1] = t+1; 
        stencil[t*4+2] = t+1; 
        stencil[t*4+3] = t+1; 
        continue;
    }
}   

__global__
void cudaDegenerateJudge(const PointVelocity3D* tets,
                    const int tetSize, 
                    Vector3D* queryvec,
                    int* stencil,
                    PointTet* criticalPoints,
                    int* pointSize)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t<tetSize; t+=dim){
        Vector3D point = {NAN,NAN,NAN};
        const PointVelocity3D* pv = &tets[t*4];

        // 判断点是否在四面体内部
        REAL a[3] = {pv[0].getVelocity().x, pv[0].getVelocity().y, pv[0].getVelocity().z};
        REAL b[3] = {pv[1].getVelocity().x, pv[1].getVelocity().y, pv[1].getVelocity().z};
        REAL c[3] = {pv[2].getVelocity().x, pv[2].getVelocity().y, pv[2].getVelocity().z};  
        REAL d[3] = {pv[3].getVelocity().x, pv[3].getVelocity().y, pv[3].getVelocity().z};
        REAL p[3] = {queryvec->x, queryvec->y, queryvec->z};

        REAL F1 = orient3d(a,b,c,d);
        REAL F2 = orient3d(a,b,d,c);
        REAL F3 = orient3d(a,c,d,b);
        REAL F4 = orient3d(b,c,d,a);
    
        REAL Fo1 = orient3d(a,b,c,p);
        REAL Fo2 = orient3d(a,b,d,p);
        REAL Fo3 = orient3d(a,c,d,p);
        REAL Fo4 = orient3d(b,c,d,p);
        if(F1>0 && Fo1<0){continue;}
        if(F1<0 && Fo1>0){continue;}
        if(F2>0 && Fo2<0){continue;}
        if(F2<0 && Fo2>0){continue;}
        if(F3>0 && Fo3<0){continue;}
        if(F3<0 && Fo3>0){continue;}
        if(F4>0 && Fo4<0){continue;}
        if(F4<0 && Fo4>0){continue;}
    
        // 四面体退化情况
        if(F1 == 0 && F2 == 0 && F3 == 0 && F4 == 0){
            if(Fo1 == 0 && Fo2 == 0 && Fo3 == 0 && Fo4 == 0){
                // CoPoint abcd
                if(pv[0].getVelocity() == pv[1].getVelocity() &&  pv[0].getVelocity() == pv[2].getVelocity() && pv[0].getVelocity() == pv[3].getVelocity())
                { 
                    if( pv[0].getVelocity().x > -EPS && pv[0].getVelocity().x < EPS && 
                        pv[0].getVelocity().y > -EPS && pv[0].getVelocity().y < EPS && 
                        pv[0].getVelocity().z > -EPS && pv[0].getVelocity().z < EPS)
                    {
                        point = pv[0].getPoint();
                        criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,4);
                    }
                    continue; 
                }
    
                // Collinear abcd
                Vector3D vab = pv[1].getVelocity()-pv[0].getVelocity();
                Vector3D vac = pv[2].getVelocity()-pv[0].getVelocity();
                Vector3D vad = pv[3].getVelocity()-pv[0].getVelocity();
                Vector3D vaq = (*queryvec) - pv[0].getVelocity();
                // if((pv[1].getVelocity()-pv[0].getVelocity() == pv[2].getVelocity()-pv[0].getVelocity() || pv[1].getVelocity()-pv[0].getVelocity() == pv[0].getVelocity()-pv[2].getVelocity()) &&
                //     (pv[1].getVelocity()-pv[0].getVelocity() == pv[3].getVelocity()-pv[0].getVelocity() || pv[1].getVelocity()-pv[0].getVelocity() == pv[0].getVelocity()-pv[3].getVelocity()))
                // { 
                //     if(pv[1].getVelocity()-pv[0].getVelocity() == queryvec-pv[0].getVelocity() || pv[1].getVelocity()-pv[0].getVelocity() == pv[0].getVelocity()-queryvec){
                if((vab == vac || vab == (vac * -1)) && (vab == vad || vab == (vad * -1))){
                    if(vab == vaq || vab == (vaq * -1)){
                        normalize(vaq);
                        Vector3D vqb = pv[1].getVelocity()-(*queryvec);
                        normalize(vqb);
                        Vector3D vqc = pv[2].getVelocity()-(*queryvec);
                        normalize(vqc);
                        Vector3D vqd = pv[3].getVelocity()-(*queryvec);
                        normalize(vqd);
                        if(vaq == vqb || vaq == vqc || vaq == vqd){
                            lineInter(pv[0].getVelocity(), pv[1].getVelocity(), pv[0].getPoint(), pv[1].getPoint(), (*queryvec), point);
                            criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,300);
                        }
                    } 
                    continue; 
                }
    
                // Coplanar abcd 面积公式
                int tp[4];
                Vector3D p[4];
                Vector3D v[4];
                tp[0] = planarInterpolation2(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                             pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                             (*queryvec),p[0]);
                tp[1] = planarInterpolation2(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                             pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                             (*queryvec),p[1]);
                tp[2] = planarInterpolation2(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                             pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                             (*queryvec),p[2]);
                tp[3] = planarInterpolation2(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                             pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                             (*queryvec),p[3]);
    
                REAL minLenth = 1000000000.0;
                int nearId = -1;
                for(int ti=0;ti<4;ti++){
                    if(tp[ti] == 1){
                        cudaInterpolation(pv,p[ti],v[ti]);
                        REAL len = lenth(v[ti]);
                        if(len < minLenth){
                            minLenth = len;
                            nearId = ti;
                        }
                    }
                }
                if(nearId == -1) {continue;}
                criticalPoints[atomicAdd(pointSize,1)] = PointTet(p[nearId],pv[0].tetId,200);
            }
            continue;
        }

        stencil[t*4] = t+1; 
        stencil[t*4+1] = t+1; 
        stencil[t*4+2] = t+1; 
        stencil[t*4+3] = t+1; 
    }
}

__global__
void cudaFacetJudge(const PointVelocity3D* tets,
                    const int tetSize, 
                    Vector3D* queryvec,
                    int* stencil,
                    PointTet* criticalPoints,
                    int* pointSize)
{
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t<tetSize; t+=dim){
        Vector3D point = {NAN,NAN,NAN};
        const PointVelocity3D* pv = &tets[t*4];

        // 顶点查找是否有可用的
        bool isVert = false;
        for(int i=0;i<4;i++){
            if( pv[i].getVelocity().x > -EPS && pv[i].getVelocity().x < EPS && 
                pv[i].getVelocity().y > -EPS && pv[i].getVelocity().y < EPS && 
                pv[i].getVelocity().z > -EPS && pv[i].getVelocity().z < EPS)
            {
                point = pv[i].getPoint();
                criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,4);
                isVert = true;
            }
        }
        if(isVert) {continue;}

        REAL a[3] = {pv[0].getVelocity().x, pv[0].getVelocity().y, pv[0].getVelocity().z};
        REAL b[3] = {pv[1].getVelocity().x, pv[1].getVelocity().y, pv[1].getVelocity().z};
        REAL c[3] = {pv[2].getVelocity().x, pv[2].getVelocity().y, pv[2].getVelocity().z};  
        REAL d[3] = {pv[3].getVelocity().x, pv[3].getVelocity().y, pv[3].getVelocity().z};
        REAL p[3] = {queryvec->x, queryvec->y, queryvec->z};

        REAL Fo1 = orient3d(a,b,c,p);
        REAL Fo2 = orient3d(a,b,d,p);
        REAL Fo3 = orient3d(a,c,d,p);
        REAL Fo4 = orient3d(b,c,d,p);

        // 请求点与四面体是否共面
        if(Fo1 == 0){
            int type = 21;
            type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                    pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                    (*queryvec), point, pv[0].tetId);
            if(type == 1){type = 21;}
            else if(type == 0){
                type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                      pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                      (*queryvec), point, pv[0].tetId);
            }
            criticalPoints[atomicAdd(pointSize,1)] = PointTet(point, pv[0].tetId,type);
            continue;
        }else if(Fo2 == 0){
            int type = 22;
            type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                    pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                    (*queryvec), point, pv[0].tetId);
            if(type == 1){type = 22;}
            else if(type == 0){
                type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                      pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                      (*queryvec), point, pv[0].tetId);
            }
            criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,type);
            continue;
        }else if(Fo3 == 0){
            int type = 23;
            type = planarInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                    pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                    (*queryvec), point, pv[0].tetId);
            if(type == 1){type = 23;}
            else if(type == 0){
                type = lineInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                      pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                      (*queryvec), point, pv[0].tetId);
            }
            criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,type);
            continue;
        }else if(Fo4 == 0){
            int type = 24;
            type = planarInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                    pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                    (*queryvec), point, pv[0].tetId);
            if(type == 1){type = 24;}
            else if(type == 0){
                type = lineInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                      pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                      (*queryvec), point, pv[0].tetId);
            }
            criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,type);
            continue;
        }

        stencil[t*4] = t+1; 
        stencil[t*4+1] = t+1; 
        stencil[t*4+2] = t+1; 
        stencil[t*4+3] = t+1; 
    }
}

__global__
void cudaTetJudge(const PointVelocity3D* tets,
                    const int tetSize, 
                    Vector3D* queryvec,
                    int* stencil,
                    PointTet* criticalPoints,
                    int* pointSize)
{
    // 四面体内部插值
    // 预备知识：
    // Wa  = Vbcdp / Vabcd
    // Wb  = Vacdp / Vabcd
    // Wc  = Vabdp / Vabcd
    // Wd  = Vabcp / Vabcd
    // Vabcd = dot(Vab, cross(Vac, Vad)) / 6
    // Ep = Wa*Ea + Wb*Eb + Wc*Ec + Wd*Ed
    // 四面体的体积 = 任意3调不共面的向量的混合积的1/6
    // 1 = Wa + Wb + Wc + Wd
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t<tetSize; t+=dim){
        Vector3D point = {NAN,NAN,NAN};
        const PointVelocity3D* pv = &tets[t*4];
        double barycoords[4];

        const Vector3D v1v0 = pv[1].getVelocity() - pv[0].getVelocity();
        const Vector3D v2v0 = pv[2].getVelocity() - pv[0].getVelocity();
        const Vector3D v3v2 = pv[3].getVelocity() - pv[2].getVelocity();
        const Vector3D pv0 = (*queryvec) - pv[0].getVelocity();

        double vol = dot(v2v0, cross(v1v0, v3v2));
        double tetvolumeinv = 1.0f / vol; // 1 / 6Vabcd

        // calculate the barycentric coordinates
        barycoords[0] = dot(pv[2].getVelocity() - (*queryvec), cross(pv[1].getVelocity() - (*queryvec), v3v2)) * tetvolumeinv;
        
        barycoords[1] = dot(v2v0, cross(pv0, v3v2)) * tetvolumeinv;

        barycoords[2] = dot(pv0, cross(v1v0, pv[3].getVelocity() - (*queryvec))) * tetvolumeinv;

        barycoords[3] = 1.0f - barycoords[0] - barycoords[1] - barycoords[2];
        if(barycoords[3] < 0) barycoords[3] = 0;
        // barycoords[3] = dot(v1v0, cross(pv[2].getVelocity()-pv[0].getVelocity(), pv[1].getVelocity() - (*queryvec))) * tetvolumeinv;
        
        // compute barycentric interpolation
        point = pv[0].getPoint() * barycoords[0] + pv[1].getPoint() * barycoords[1]
                + pv[2].getPoint() * barycoords[2] + pv[3].getPoint() * barycoords[3];
        

        if(isinf(point.x) || isinf(point.y) || isinf(point.z) || 
            isnan(point.x) || isnan(point.y) || isnan(point.z)){continue;}

        criticalPoints[atomicAdd(pointSize,1)] = PointTet(point,pv[0].tetId,1);

    }
}


struct is_positive{
    __host__ __device__
	bool operator()(const int x)
	{
		return x >=0;
	}
};

int compactNegative(PointVelocity3D* out, PointVelocity3D* in, int* inValid, int size)
{
    // thrust::counting_iterator<int> pos_pointer(0);	
    thrust::device_ptr<PointVelocity3D> in_pointer((PointVelocity3D*)in);
	thrust::device_ptr<PointVelocity3D> out_pointer((PointVelocity3D*)out);	
	thrust::device_ptr<int> stencil_pointer((int*)inValid);
	thrust::device_ptr<PointVelocity3D> final_pointer = 
		thrust::copy_if(in_pointer, in_pointer + size, stencil_pointer, out_pointer, is_positive());
	int num = final_pointer - out_pointer;
	return num;
};


GpuSSCPFinder::GpuSSCPFinder(const char* filename) : BaseCPFinder(filename)
{
    initialization();
}

GpuSSCPFinder::~GpuSSCPFinder()
{
    checkCudaErrors(cudaFree(d_tets1));
    checkCudaErrors(cudaFree(d_tets2));
    checkCudaErrors(cudaFree(d_outstencil));
    checkCudaErrors(cudaFree(d_outpoints));
    checkCudaErrors(cudaFree(d_pointsize));
    checkCudaErrors(cudaFree(d_queryvec));
}


void GpuSSCPFinder::initialization()
{
    // Exactinit
    const int PredParaNum = 14;	        // the size of hostConst
    REAL hostConst[14];					// constants used in exact CUDA CCW test

	// Exact computation
	REAL *d_constants;
    checkCudaErrors(cudaMalloc((void **)&d_constants, PredParaNum * sizeof(REAL)) ); 
    cudaExactinit<<<1,1>>>(d_constants); 
	checkCudaErrors(cudaMemcpy(hostConst, d_constants, PredParaNum * sizeof(REAL), cudaMemcpyDeviceToHost) ); 
	
    checkCudaErrors(cudaMemcpyToSymbol(deviceConst, hostConst, PredParaNum * sizeof(REAL)) ); 
    checkCudaErrors(cudaFree(d_constants));

    // Tet init
    tetnum = seqtets.size() / 4;

    checkCudaErrors(cudaMalloc((void **)&d_tets1, seqtets.size() * sizeof(PointVelocity3D)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_tets2, seqtets.size() * sizeof(PointVelocity3D)) ); 
    checkCudaErrors(cudaMemcpy(d_tets1, seqtets.data(), seqtets.size() * sizeof(PointVelocity3D), cudaMemcpyHostToDevice) ); 

    checkCudaErrors(cudaMalloc((void **)&d_outstencil, tetnum * 4 * sizeof(int)) );
    checkCudaErrors(cudaMalloc((void **)&d_outpoints, tetnum * 4 * sizeof(PointTet)) ); 
    checkCudaErrors(cudaMalloc((void **)&d_pointsize, sizeof(int)) );

    checkCudaErrors(cudaMemset(d_pointsize,0,sizeof(int))) ;
    
    cudaInit<<<100,256>>>(d_outstencil, d_outpoints, tetnum);

    Vector3D query = {0.0,0.0,0.0};
    checkCudaErrors(cudaMalloc((void **)&d_queryvec, sizeof(Vector3D)) );
    checkCudaErrors(cudaMemcpy(d_queryvec, &query, sizeof(Vector3D), cudaMemcpyHostToDevice) ); 
}


void GpuSSCPFinder::search()
{
    dim3 grid(64,1,1);
    dim3 block(256,1,1);
    int tetsVertNum1 = 0;
    int tetsVertNum2 = 0;
    double loadtime = 0;

    while(true){
        // Rough judgment
        cudaRoughJudge<<<grid,block>>>(d_tets1, tetnum, d_queryvec, d_outstencil,d_outpoints, d_pointsize);
        checkCudaErrors(cudaDeviceSynchronize());
        tetsVertNum2 = compactNegative(d_tets2,d_tets1,d_outstencil,tetnum*4);
        checkCudaErrors(cudaMemset(d_outstencil,-1,tetnum * 4 * sizeof(int)));

        // Degradation judgment
        cudaDegenerateJudge<<<grid,block>>>(d_tets2, tetsVertNum2/4, d_queryvec, d_outstencil, d_outpoints, d_pointsize);
        checkCudaErrors(cudaDeviceSynchronize());
        tetsVertNum1 = compactNegative(d_tets1,d_tets2,d_outstencil,tetsVertNum2);
        checkCudaErrors(cudaMemset(d_outstencil,-1,tetsVertNum2 * sizeof(int)));

        // Vertex line and facet judgment 
        cudaFacetJudge<<<grid,block>>>(d_tets1, tetsVertNum1/4, d_queryvec, d_outstencil, d_outpoints, d_pointsize);
        checkCudaErrors(cudaDeviceSynchronize());
        tetsVertNum2 = compactNegative(d_tets2,d_tets1,d_outstencil,tetsVertNum1);
        checkCudaErrors(cudaMemset(d_outstencil,-1,tetsVertNum1 * sizeof(int)));

        // Inside the tetrahedron
        cudaTetJudge<<<grid,block>>>(d_tets2, tetsVertNum2/4, d_queryvec, d_outstencil, d_outpoints, d_pointsize);
        checkCudaErrors(cudaDeviceSynchronize());
        tetsVertNum1 = compactNegative(d_tets1,d_tets2,d_outstencil,tetsVertNum2);
        checkCudaErrors(cudaMemset(d_outstencil,-1,tetsVertNum2 * sizeof(int)));

        if(isBigFile){
            auto start = std::chrono::system_clock::now();
            auto ret = loadDataBlock();
            tetnum = seqtets.size() / 4;
            checkCudaErrors(cudaMemcpy(d_tets1, seqtets.data(), seqtets.size() * sizeof(PointVelocity3D), cudaMemcpyHostToDevice) ); 
            auto end = std::chrono::system_clock::now();
            auto diff = end - start;
            loadtime += std::chrono::duration <double, std::milli> (diff).count() / 1000;
            std::cout << "    - " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::flush;
            if(ret){break;}
        }else{ break;}
    }
    if(loadtime > 0) std::cout << "    total ld time: " << loadtime << " seconds" << std::endl;
    int pointSize = 0;
    checkCudaErrors(cudaMemcpy(&pointSize, d_pointsize, sizeof(int), cudaMemcpyDeviceToHost) ); 
    criticalPoints.resize(pointSize);
    checkCudaErrors(cudaMemcpy(criticalPoints.data(), d_outpoints, pointSize * sizeof(PointTet), cudaMemcpyDeviceToHost) );
    
}