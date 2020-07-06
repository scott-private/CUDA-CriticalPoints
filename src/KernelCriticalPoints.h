/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : KernelCriticalPoints.h
============================================================================
*/

#ifndef KERNEL_CRITICAL_POINTS_H
#define KERNEL_CRITICAL_POINTS_H

#include "GpuPredicates.h"

__device__
void cudaInterpolation(const PointVelocity3D pv[4], Vector3D point, Vector3D& velocity)
{
    // 第一个点表示这个区域角的数量，---》四面体
    double barycoords[4];
    Vector3D cellpoints[4] = { // 读取四面体的四个顶点
        pv[0].getPoint(),
        pv[1].getPoint(),
        pv[2].getPoint(),
        pv[3].getPoint() 
    };

    // 预备知识：
    // Wa  = Vbcdp / Vabcd
    // Wb  = Vacdp / Vabcd
    // Wc  = Vabdp / Vabcd
    // Wd  = Vabcp / Vabcd
    // Vabcd = dot(Vab, cross(Vac, Vad)) / 6
    // Ep = Wa*Ea + Wb*Eb + Wc*Ec + Wd*Ed
    // 四面体的体积 = 任意3调不共面的向量的混合积的1/6
    // 1 = Wa + Wb + Wc + Wd

    const Vector3D v1v0 = cellpoints[1] - cellpoints[0];
    const Vector3D v2v0 = cellpoints[2] - cellpoints[0];
    const Vector3D v3v2 = cellpoints[3] - cellpoints[2];
    const Vector3D pv0 = point - cellpoints[0];

    double test = dot(v2v0, cross(v1v0, v3v2));
    if(test == 0) {velocity={-10000,-10000,-10000}; return ;}
    double tetvolumeinv = 1.0f / test;// 1 / 6Vabcd

    // calculate the barycentric coordinates
    barycoords[0] = dot(cellpoints[2] - point, cross(cellpoints[1] - point, v3v2)) * tetvolumeinv;
    barycoords[1] = dot(v2v0, cross(pv0, v3v2)) * tetvolumeinv;
    barycoords[2] = dot(pv0, cross(v1v0, cellpoints[3] - point)) * tetvolumeinv;
    barycoords[3] = 1.0 - barycoords[0] - barycoords[1] - barycoords[2];
    
    // compute barycentric interpolation
    velocity = pv[0].getVelocity() * barycoords[0] +
                pv[1].getVelocity() * barycoords[1] +
                pv[2].getVelocity() * barycoords[2] +
                pv[3].getVelocity() * barycoords[3];
}

__device__ __forceinline__
void lineInter(Vector3D va, Vector3D vb,
               Vector3D pa, Vector3D pb, 
               Vector3D queryvec, Vector3D& point)
{
    double lenap = lenth(queryvec-va); 
    double lenbp = lenth(queryvec-vb);
    double maxlen = lenap + lenbp;
    if(pb.x > pa.x){ 
        point.x = (pb.x - pa.x) * (lenap / maxlen) + pa.x;
        point.y = (pb.y - pa.y) * (lenap / maxlen) + pa.y;
        point.z = (pb.z - pa.z) * (lenap / maxlen) + pa.z;
    }else{
        point.x = (pa.x - pb.x) * (lenbp / maxlen) + pb.x;
        point.y = (pa.y - pb.y) * (lenbp / maxlen) + pb.y; 
        point.z = (pa.z - pb.z) * (lenbp / maxlen) + pb.z; 
    }
}

__device__
int lineInterpolation(Vector3D va, Vector3D vb, Vector3D vc, 
                      Vector3D pa, Vector3D pb, Vector3D pc, 
                      Vector3D queryvec, Vector3D& point , int tetId)
{
    int collinearType = 30;
  
    if(point_on_seg(queryvec,Line3D(vb,vc))){
        // if(tetId ==  1377834) std::cout << "--3" << std::endl;
        collinearType = 33;
        lineInter(vb,vc,pb,pc,queryvec,point);
    }else if(point_on_seg(queryvec,Line3D(va,vb))){
        // if(tetId ==  1377834) std::cout << "--1" << std::endl;
        collinearType = 31;
        lineInter(va,vb,pa,pb,queryvec,point);
    }else if(point_on_seg(queryvec,Line3D(va,vc))){
        // if(tetId ==  1377834) std::cout << "--2" << std::endl;
        collinearType = 32;
        lineInter(va,vc,pa,pc,queryvec,point);
    }else{
        Vector3D p1 = {NAN,NAN,NAN};
        Vector3D p2 = {NAN,NAN,NAN};
        Vector3D p3 = {NAN,NAN,NAN};
        lineInter(va,vb,pa,pb,queryvec,p1);
        lineInter(va,vc,pa,pc,queryvec,p2);
        lineInter(vb,vc,pb,pc,queryvec,p3);
        // collinearType = 0;
    }
    // if(tetId ==  1377834) std::cout << "--4: " << point.x << "," << point.y <<"," << point.z << std::endl;
    return collinearType;
}

__device__ 
double count_triangle_area(Vector3D a,Vector3D b,Vector3D c)
{
	double area = -1;
	 
	double side[3];//存储三条边的长度;
 
	side[0] = sqrt(pow(a.x - b.x,2) + pow(a.y - b.y,2) + pow(a.z - b.z,2)); 
	side[1] = sqrt(pow(a.x - c.x,2) + pow(a.y - c.y,2) + pow(a.z - c.z,2));
	side[2] = sqrt(pow(c.x - b.x,2) + pow(c.y - b.y,2) + pow(c.z - b.z,2)); 
 
	//不能构成三角形;
	if(side[0]+side[1]<=side[2] || side[0]+side[2]<=side[1] || side[1]+side[2]<=side[0]){return area;} 
 
	//利用海伦公式。s=sqr(p*(p-a)(p-b)(p-c)); 
	double p = (side[0]+side[1]+side[2])/2; //半周长;
	area = sqrt(p*(p-side[0])*(p-side[1])*(p-side[2])); 
	
	return area;
}


// 预备知识：
// Wa  = Sbcp / Sabc
// Wb  = Sacp / Sabc
// Wc  = Sabp / Sabc
// Ep = Wa*Ea + Wb*Eb + Wc*Ec
// 1 = Wa + Wb + Wc
__device__
int planarInterpolation(Vector3D va, Vector3D vb, Vector3D vc, 
                        Vector3D pa, Vector3D pb, Vector3D pc, 
                        Vector3D queryvec, Vector3D& point, int tetId)
{
    double Sabc = count_triangle_area(va,vb,vc);
    double Sbcv = count_triangle_area(vb,vc,queryvec);
    double Sacv = count_triangle_area(va,vc,queryvec);
    double Sabv = count_triangle_area(va,vb,queryvec);
    // if(tetId ==  1377834) std::cout  << "--1-0 " << Sabc << " " << Sbcv << " " << Sacv << " " << Sabv  << std::endl;

    if(Sabc == -1){return 0;}
    if(Sbcv == -1){lineInter(vb,vc,pb,pc,queryvec,point); return 31; }
    if(Sacv == -1){lineInter(va,vc,pa,pc,queryvec,point); return 32; }
    if(Sabv == -1){lineInter(va,vb,pa,pb,queryvec,point); return 33; }

    // if(tetId ==  1377834) std::cout << "--1*1 " << std::endl;
    double Wa = (Sbcv / Sabc);
    // if(tetId ==  1377834) std::cout << "--1*2:"<< Wa << std::endl;
    if(Wa < 0 && Wa > 1) {return 0;}
    
    double Wb = (Sacv / Sabc);
    // if(tetId ==  1377834) std::cout << "--1*3:"<< Wb << std::endl;
    if(Wb < 0 && Wb > 1) {return 0;}
    
    double Wc = 1 - Wa - Wb;
    // double Wc = (Sabv / Sabc);
    // if(tetId ==  1377834) std::cout << "--1*4:"<< Wc << std::endl;
    if(Wc < 0 && Wc > 1) {return 0;}

    point = pa*Wa + pb*Wb + pc*Wc;

    // if(tetId ==  1377834) std::cout << "--1*5:(" << point.x << " " << point.y << " " << point.z<< ")" << std::endl;
    return 1;
}

__device__ 
int planarInterpolation2(Vector3D va, Vector3D vb, Vector3D vc, 
                         Vector3D pa, Vector3D pb, Vector3D pc, 
                         Vector3D queryvec, Vector3D& point)
{
    double Sabc = count_triangle_area(va,vb,vc);
    double Sbcv = count_triangle_area(vb,vc,queryvec);
    double Sacv = count_triangle_area(va,vc,queryvec);
    double Sabv = count_triangle_area(va,vb,queryvec);

    if(Sabc < EPS && Sabc > -EPS ){ return 0; }
    else if(Sabc == -1 ||Sbcv == -1 || Sacv == -1 || Sabv == -1){ return 0; }
    else if(Sabc - (Sbcv + Sacv + Sabv) < -EPS && Sabc - (Sbcv + Sacv + Sabv) > EPS){ return 0; }

    double Wa = (Sbcv / Sabc);
    if(Wa < 0 && Wa > 1) {return 0;}
    
    double Wb = (Sacv / Sabc);
    if(Wb < 0 && Wb > 1) {return 0;}
    
    double Wc = 1 - Wa - Wb;
    // double Wc = (Sabv / Sabc);
    if(Wc < 0 && Wc > 1) {return 0;}

    point = pa*Wa + pb*Wb + pc*Wc;
    return 1;
}

__device__
bool findCriticalPoint2(const PointVelocity3D pv[4], 
                        Vector3D queryvec,
                        const int  tetId,
                        PointTet* criticalPoints)
{
    int cpnum = 0;
    Vector3D point = {NAN,NAN,NAN};
    // 第一个点表示这个区域角的数量，---》四面体
    double barycoords[4];
    Vector3D cellvecters[4] = { // 读取四面体的四个顶点
        pv[0].getVelocity(), 
        pv[1].getVelocity(), 
        pv[2].getVelocity(), 
        pv[3].getVelocity()
    };

    // critical point is vertex
    // 顶点查找是否有可用的
    bool isVert = false;
    for(int i=0;i<4;i++){
        if( cellvecters[i].x > -EPS && cellvecters[i].x < EPS && 
            cellvecters[i].y > -EPS && cellvecters[i].y < EPS && 
            cellvecters[i].z > -EPS && cellvecters[i].z < EPS){
            point = pv[i].getPoint();
            criticalPoints[cpnum++] = PointTet(point,tetId,4);
            isVert = true;
        }
    }
    if(isVert) {return true;}

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
    if(xMax < queryvec.x || xMin > queryvec.x){return false;}
    if(yMax < queryvec.y || yMin > queryvec.y){return false;}
    if(zMax < queryvec.z || zMin > queryvec.z){return false;}

    // 判断点是否在四面体内部
    REAL a[3] = {pv[0].getVelocity().x, pv[0].getVelocity().y, pv[0].getVelocity().z};
    REAL b[3] = {pv[1].getVelocity().x, pv[1].getVelocity().y, pv[1].getVelocity().z};
    REAL c[3] = {pv[2].getVelocity().x, pv[2].getVelocity().y, pv[2].getVelocity().z};  
    REAL d[3] = {pv[3].getVelocity().x, pv[3].getVelocity().y, pv[3].getVelocity().z};
    REAL p[3] = {queryvec.x, queryvec.y, queryvec.z};
    REAL F1 = orient3d(a,b,c,d);
    REAL F2 = orient3d(a,b,d,c);
    REAL F3 = orient3d(a,c,d,b);
    REAL F4 = orient3d(b,c,d,a);

    REAL Fo1 = orient3d(a,b,c,p);
    REAL Fo2 = orient3d(a,b,d,p);
    REAL Fo3 = orient3d(a,c,d,p);
    REAL Fo4 = orient3d(b,c,d,p);
    if(F1>0 && Fo1<0){return false;}
    if(F1<0 && Fo1>0){return false;}
    if(F2>0 && Fo2<0){return false;}
    if(F2<0 && Fo2>0){return false;}
    if(F3>0 && Fo3<0){return false;}
    if(F3<0 && Fo3>0){return false;}
    if(F4>0 && Fo4<0){return false;}
    if(F4<0 && Fo4>0){return false;}

    // 四面体退化情况
    if(F1 == 0 && F2 == 0 && F3 == 0 && F4 == 0){
        if(Fo1 == 0 && Fo2 == 0 && Fo3 == 0 && Fo4 == 0){
            // CoPoint abcd
            if(cellvecters[0] == cellvecters[1] &&  cellvecters[0] == cellvecters[2] && cellvecters[0] == cellvecters[3])
                { return false; }

            // Collinear abcd
            Vector3D vab = cellvecters[1]-cellvecters[0];
            Vector3D vac = cellvecters[2]-cellvecters[0];
            Vector3D vad = cellvecters[3]-cellvecters[0];
            Vector3D vaq = queryvec-cellvecters[0];
            // if((cellvecters[1]-cellvecters[0] == cellvecters[2]-cellvecters[0] || cellvecters[1]-cellvecters[0] == cellvecters[0]-cellvecters[2]) &&
            //     (cellvecters[1]-cellvecters[0] == cellvecters[3]-cellvecters[0] || cellvecters[1]-cellvecters[0] == cellvecters[0]-cellvecters[3]))
            // { 
            //     if(cellvecters[1]-cellvecters[0] == queryvec-cellvecters[0] || cellvecters[1]-cellvecters[0] == cellvecters[0]-queryvec){
            if((vab == vac || vab == (vac * -1)) && (vab == vad || vab == (vad * -1))){
                if(vab == vaq || vab == (vaq * -1)){
                    lineInter(pv[0].getVelocity(), pv[1].getVelocity(), pv[0].getPoint(), pv[1].getPoint(), queryvec, point);
                    criticalPoints[cpnum++] = PointTet(point,tetId,300);
                    return true;
                }else{ return false;}
            }

            // Coplanar abcd 面积公式
            int tp[4];
            Vector3D p[4];
            Vector3D v[4];
            tp[0] = planarInterpolation2(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                         pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                         queryvec,p[0]);
            tp[1] = planarInterpolation2(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                         pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                         queryvec,p[1]);
            tp[2] = planarInterpolation2(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                         pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                         queryvec,p[2]);
            tp[3] = planarInterpolation2(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                         pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                         queryvec,p[3]);


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
            if(nearId == -1) {return false;}
            criticalPoints[cpnum++] = PointTet(p[nearId],tetId,200);
            return true;
        }
        return false;
    }


    // 请求点与四面体是否共面
    if(Fo1 == 0){
        int type = 21;
        type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 21;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                  pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                  queryvec, point, tetId);
        }
        criticalPoints[cpnum++] = PointTet(point,tetId,type);
        return true;
    }else if(Fo2 == 0){
        int type = 22;
        type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 22;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                  pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
        }
        criticalPoints[cpnum++] = PointTet(point,tetId,type);
        return true;
    }else if(Fo3 == 0){
        int type = 23;
        type = planarInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 23;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                  pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
        }
        criticalPoints[cpnum++] = PointTet(point,tetId,type);
        return true; 
    }else if(Fo4 == 0){
        int type = 24;
        type = planarInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 24;}
        else if(type == 0){
            type = lineInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                  pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
        }
        criticalPoints[cpnum++] = PointTet(point,tetId,type);
        return true;
    }
    
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

    const Vector3D v1v0 = cellvecters[1] - cellvecters[0];
    const Vector3D v2v0 = cellvecters[2] - cellvecters[0];
    const Vector3D v3v2 = cellvecters[3] - cellvecters[2];
    const Vector3D pv0 = queryvec - cellvecters[0];

    double vol = dot(v2v0, cross(v1v0, v3v2));
    double tetvolumeinv = 1.0f / vol; // 1 / 6Vabcd

    // calculate the barycentric coordinates
    barycoords[0] = dot(cellvecters[2] - queryvec, cross(cellvecters[1] - queryvec, v3v2)) * tetvolumeinv;
    // if(barycoords[0] < 0 || barycoords[0] > 1) {return false;}
    
    barycoords[1] = dot(v2v0, cross(pv0, v3v2)) * tetvolumeinv;
    // if(barycoords[1] < 0 || barycoords[1] > 1) {return false;}

    barycoords[2] = dot(pv0, cross(v1v0, cellvecters[3] - queryvec)) * tetvolumeinv;
    // if(barycoords[2] < 0 || barycoords[2] > 1) {return false;}

    barycoords[3] = 1.0f - barycoords[0] - barycoords[1] - barycoords[2];
    if(barycoords[3] < 0) barycoords[3] = 0;
    // barycoords[3] = dot(v1v0, cross(cellvecters[2]-cellvecters[0], cellvecters[1] - queryvec)) * tetvolumeinv;
    // if(barycoords[3] < 0 || barycoords[3] > 1) {return false;}
    
    // compute barycentric interpolation
    point = pv[0].getPoint() * barycoords[0] +
            pv[1].getPoint() * barycoords[1] +
            pv[2].getPoint() * barycoords[2] +
            pv[3].getPoint() * barycoords[3];
    
    // if(isNaN(fabs(point.x)) || isNaN(fabs(point.y)) || isNaN(fabs(point.z))){return false;}
    if(isinf(point.x) || isinf(point.y) || isinf(point.z) || 
        isnan(point.x) || isnan(point.y) || isnan(point.z)){return false;}

    criticalPoints[cpnum++] = PointTet(point,tetId,1);
    return true;
}

__global__
void cudaSearch(const PointVelocity3D* tets,
                const int tetSize, 
                Vector3D* velocity,
                PointTet* criticalPoints)
{
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t<tetSize; t+=dim){
        findCriticalPoint2(&tets[t*4], *velocity, t, &criticalPoints[t*4]);
    }
}   


__global__
void cudaCompress(PointTet* criticalPoints,
                  const int tetSize, 
                  PointTet* outPoints,
                  int* outPointSize)
{
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;

    for(int t=tid; t<tetSize; t+=dim){
        int pidx = t*4;
        if(criticalPoints[pidx].tetId != -1) {int pIdx = atomicAdd(outPointSize , 1);outPoints[pIdx] =  criticalPoints[pidx];}
        if(criticalPoints[pidx+1].tetId != -1) {int pIdx = atomicAdd(outPointSize , 1);outPoints[pIdx] =  criticalPoints[pidx+1];}
        if(criticalPoints[pidx+2].tetId != -1) {int pIdx = atomicAdd(outPointSize , 1);outPoints[pIdx] =  criticalPoints[pidx+2];}
        if(criticalPoints[pidx+3].tetId != -1) {int pIdx = atomicAdd(outPointSize , 1);outPoints[pIdx] =  criticalPoints[pidx+3];}
    }
}


__global__
void cudaInit(int* outstencil, PointTet* criticalPoints, const int tetSize)
{
    int tid = blockIdx.x * blockDim.x +  threadIdx.x;
    int dim = gridDim.x * blockDim.x;
    for(int t=tid; t < tetSize; t+=dim){
        int idx = t*4;
        outstencil[idx] = -1;
        outstencil[idx+1] = -1;
        outstencil[idx+2] = -1;
        outstencil[idx+3] = -1;
        
        criticalPoints[idx].tetId = -1;
        criticalPoints[idx+1].tetId = -1;
        criticalPoints[idx+2].tetId = -1;
        criticalPoints[idx+3].tetId = -1;
    }
}



__global__ 
void cudaExactinit(REAL *constants)
{
    REAL half;
    REAL check, lastcheck;
    REAL cuda_epsilon = 1.0; 
    REAL cuda_splitter = 1.0; 
    int every_other;

    every_other = 1;
    half = 0.5;
    check = 1.0;
    // Repeatedly divide `cuda_epsilon' by two until it is too small to add to      
    //   one without causing roundoff.  (Also check if the sum is equal to     
    //   the previous sum, for machines that round up instead of using exact   
    //   rounding.  Not that these routines will work on such machines.)       
    do {
        lastcheck = check;
        cuda_epsilon *= half;
        if (every_other) {
            cuda_splitter *= 2.0;
        }
        every_other = !every_other;
        check = 1.0 + cuda_epsilon;
    } while ((check != 1.0) && (check != lastcheck));
    constants[0] = cuda_splitter + 1.0;											//splitter
    constants[1] = cuda_epsilon;												//epsilon											
    /* Error bounds for orientation and incircle tests. */
    constants[2] = (3.0 + 8.0 * cuda_epsilon) * cuda_epsilon;					//resulterrbound
    constants[3] = (3.0 + 16.0 * cuda_epsilon) * cuda_epsilon;					//ccwerrboundA
    constants[4] = (2.0 + 12.0 * cuda_epsilon) * cuda_epsilon;					//ccwerrboundB
    constants[5] = (9.0 + 64.0 * cuda_epsilon) * cuda_epsilon * cuda_epsilon;	//ccwerrboundC
    constants[6] = (10.0 + 96.0 * cuda_epsilon) * cuda_epsilon;					//iccerrboundA
    constants[7] = (4.0 + 48.0 * cuda_epsilon) * cuda_epsilon;					//iccerrboundB
    constants[8] = (44.0 + 576.0 * cuda_epsilon) * cuda_epsilon * cuda_epsilon;	//iccerrboundC
    constants[9] = (7.0 + 56.0 * cuda_epsilon) * cuda_epsilon;					//o3derrboundA
    constants[10] = (3.0 + 28.0 * cuda_epsilon) * cuda_epsilon;					//o3derrboundB
    constants[11] = (26.0 + 288.0 * cuda_epsilon) * cuda_epsilon * cuda_epsilon;//o3derrboundC
    constants[12] = 1.0 / sin(REAL(0.0));										// infinity
    constants[13] = (12.0 + 128.0 * cuda_epsilon) * cuda_epsilon;				// Weighted incircle
}

#endif
