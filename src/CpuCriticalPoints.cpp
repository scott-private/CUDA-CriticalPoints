/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : CpuCriticalPoints.cpp
============================================================================
*/

#include "CriticalPoints.h"
#include "CpuPredicates.h"
#include <cassert>
#include <iostream>
#include <vector>
#include <float.h>


inline Vector3D inverseInterpolation(const PointVelocity3D& p1, 
                                     const PointVelocity3D& p2, 
                                     const REAL velocity[3], 
                                     int direction )
{
    REAL v1 = p1[3+direction];
    REAL v2 = p2[3+direction];

    REAL x, y, z;

    if (v2 == v1) {
        x = (p1.x() + p2.x()) / 2.0f;
        y = (p1.y() + p2.y()) / 2.0f;
        z = (p1.z() + p2.z()) / 2.0f;
    } else {
        /*
         <----+-----+---+----->
              v1    |   v2
                 isolevel
         <----+-----+---+----->
              0     |   1
                  interp
         */
        // interp == 0: vert should be at p1
        // interp == 1: vert should be at p2
        REAL interp = (velocity[direction] - v1) / (v2 - v1);
        REAL oneMinusInterp = 1 - interp;

        x = p1.x() * oneMinusInterp + p2.x() * interp;
        y = p1.y() * oneMinusInterp + p2.y() * interp;
        z = p1.z() * oneMinusInterp + p2.z() * interp;
    }

    return Vector3D(x,y,z);
}

// check is critical point
inline int checkCriticalPoint(const Plane3D p1,
                              const Plane3D p2, 
                              const Plane3D p3,
                              const Vector3D tet[4], 
                              Vector3D& point)
{
    
    if(parallel(p1,p2) ){ // Plane to plane parallel
        if(point_on_plane(p1.a,p2) && point_on_plane(p1.b,p2) && point_on_plane(p1.b,p2)  && !parallel(p1,p3)){ // Plane1 Plane2 Coplanar
            Line3D line = plane_intersection(p1,p3);
            if(dist(line.a,line.b) == 0) {return 0;} 
            point = (line.a+line.b) * 0.5;
            return 3;
        }else {return 0;}
    }else{
        // 两个平面确定一个交线
        Line3D line = plane_intersection(p1,p2);
       
        int sig = line_to_plane(line,p3);
        if(sig == -1){// line on plane
            if(dist(line.a,line.b) == 0) {return 0;};
            point = (line.a+line.b) * 0.5;
            return 2;
        } else if(sig == 0){return 0;}// parallel : no cross point

        // 交线与面确定一个交点
        point = line_plane_intersection(line,p3);
        // 交点是否在四面体内部 可能会有一点点不准确
        // if( same_side(point,tet[3],Plane3D(tet[0],tet[1],tet[2])) && 
        //     same_side(point,tet[1],Plane3D(tet[0],tet[2],tet[3])) && 
        //     same_side(point,tet[0],Plane3D(tet[1],tet[3],tet[2])) && 
        //     same_side(point,tet[2],Plane3D(tet[1],tet[0],tet[3])))
        // {return 1;}
        
        // 体积法
        REAL Vabcd = volume(tet[0],tet[1],tet[2],tet[3]);
        REAL Vabcp = volume(tet[0],tet[1],tet[2],point);
        REAL Vabpd = volume(tet[0],tet[1],point,tet[3]);
        REAL Vapcd = volume(tet[0],point,tet[2],tet[3]);
        REAL Vpbcd = volume(point,tet[1],tet[2],tet[3]);
        return (Vabcd - (Vabcp + Vabpd + Vapcd + Vpbcd)) >= -(EPS);
    }
    return 0;
}

bool findCriticalPoint(const PointVelocity3D pv[4], 
                       const REAL velocity[3],
                       const int  tetId,
                       std::vector<PointTet>& criticalPoints)
{

    /*
     Tetrahedron layout:

           0
           *
          /|\
         / | \
      3 *-----* 1
         \ | /
          \|/
           *
           2
    */

    int state[3] = {-1,-1,-1};
    Vector3D triPos[3][6]; 

    // 这里还需要对3个方向所得的xyz进行比较，如果一样，则说明有关键点
    for(int d = 0; d < 3; d++)
    {
        unsigned char index = 0;
        for (int i = 0; i < 4; ++i)
            if (pv[i][3+d] <= velocity[d])
                index |= (1 << i);

        switch (index) {
            // we don't do anything if everyone is inside or outside
            case 0x00:
            case 0x0F:
                return false;

            // only vert 0 is inside
            case 0x01:
                triPos[d][0] = inverseInterpolation(pv[0], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[0], pv[3], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[0], pv[2], velocity, d);
                state[d] = 1;
                break;

            // only vert 1 is inside
            case 0x02:
                triPos[d][0] = inverseInterpolation(pv[1], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[1], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[3], velocity, d);
                state[d] = 1;
                break;

            // only vert 2 is inside
            case 0x04:
                triPos[d][0] = inverseInterpolation(pv[2], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[2], pv[3], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[2], pv[1], velocity, d);
                state[d] = 1;
                break;

            // only vert 3 is inside
            case 0x08:
                triPos[d][0] = inverseInterpolation(pv[3], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[3], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[3], pv[0], velocity, d);
                state[d] = 1;
                break;

            // verts 0, 1 are inside
            case 0x03:
                triPos[d][0] = inverseInterpolation(pv[3], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[2], pv[0], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[3], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[2], pv[0], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[2], pv[1], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[1], pv[3], velocity, d);
                state[d] = 2;
                break;

            // verts 0, 2 are inside
            case 0x05:
                triPos[d][0] = inverseInterpolation(pv[3], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[1], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[0], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[1], pv[2], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[3], pv[0], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[2], pv[3], velocity, d);
                state[d] = 2;
                break;

            // verts 0, 3 are inside
            case 0x09:
                triPos[d][0] = inverseInterpolation(pv[0], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[0], pv[2], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[3], pv[2], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[0], pv[2], velocity, d);
                state[d] = 2;
                break;

            // verts 1, 2 are inside
            case 0x06:
                triPos[d][0] = inverseInterpolation(pv[0], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[0], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[3], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[0], pv[2], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[3], pv[2], velocity, d);
                state[d] = 2;
                break;

            // verts 2, 3 are inside
            case 0x0C:
                triPos[d][0] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[2], pv[0], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[3], pv[0], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[2], pv[0], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[2], pv[1], velocity, d);
                state[d] = 2;
                break;

            // verts 1, 3 are inside
            case 0x0A:
                triPos[d][0] = inverseInterpolation(pv[3], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[1], pv[0], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[2], velocity, d);

                triPos[d][3] = inverseInterpolation(pv[1], pv[2], velocity, d);
                triPos[d][4] = inverseInterpolation(pv[2], pv[3], velocity, d);
                triPos[d][5] = inverseInterpolation(pv[3], pv[0], velocity, d);
                state[d] = 2;
                break;

            // verts 0, 1, 2 are inside
            case 0x07:
                triPos[d][0] = inverseInterpolation(pv[3], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[3], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[3], pv[1], velocity, d);
                state[d] = 1;
                break;

            // verts 0, 1, 3 are inside
            case 0x0B:
                triPos[d][0] = inverseInterpolation(pv[2], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[2], pv[3], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[2], pv[0], velocity, d);
                state[d] = 1;
                break;

            // verts 0, 2, 3 are inside
            case 0x0D:
                triPos[d][0] = inverseInterpolation(pv[1], pv[0], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[1], pv[3], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[1], pv[2], velocity, d);
                state[d] = 1;
                break;

            // verts 1, 2, 3 are inside
            case 0x0E:
                triPos[d][0] = inverseInterpolation(pv[0], pv[1], velocity, d);
                triPos[d][1] = inverseInterpolation(pv[0], pv[2], velocity, d);
                triPos[d][2] = inverseInterpolation(pv[0], pv[3], velocity, d);
                state[d] = 1;
                break;

            // what is this I don't even
            default:
                assert(false);
                return false;
        }
    }
    
    Vector3D tet[4] = {pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint()};

    bool found = false;
    if(state[0] != -1 && state[1] != -1 && state[2] != -1){
        for(int i=0;i<state[0];i++){
            Plane3D p1(triPos[0][i*3],triPos[0][i*3+1],triPos[0][i*3+2]);
            for(int j=0;j<state[1];j++){
                Plane3D p2(triPos[1][j*3],triPos[1][j*3+1],triPos[1][j*3+2]);
                for(int k=0;k<state[2];k++){
                    Plane3D p3(triPos[2][k*3],triPos[2][k*3+1],triPos[2][k*3+2]);
                    Vector3D point(0,0,0);
                    int check = checkCriticalPoint(p1,p2,p3,tet,point);
                    if(check > 0){
                        found = true; 
                        criticalPoints.push_back(PointTet(point,tetId,check)); 
                        if(check == 2) {std::cout << "TetId:" << tetId << " -> Line on Plane: There are critical segments" << std::endl;}
                        if(check == 3) {std::cout << "TetId:" << tetId << " -> Plane Parallel: There are critical segments" << std::endl;}
                        goto out;
                    }
                }
            }
        } 
        out: 
        ;
    }
    return found;
}


inline void lineInter(Vector3D va, Vector3D vb,
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

inline int lineInterpolation(Vector3D va, Vector3D vb, Vector3D vc, 
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

inline double count_triangle_area(Vector3D a,Vector3D b,Vector3D c){
	double area = -1;
	 
	double side[3];//存储三条边的长度;
 
	side[0] = sqrt(pow(a.x - b.x,2)+pow(a.y - b.y,2) + pow(a.z - b.z,2)); 
	side[1] = sqrt(pow(a.x - c.x,2)+pow(a.y - c.y,2) + pow(a.z - c.z,2));
	side[2] = sqrt(pow(c.x - b.x,2)+pow(c.y - b.y,2) + pow(c.z - b.z,2)); 
 
	//不能构成三角形;
	if(side[0]+side[1]<=side[2] || side[0]+side[2]<=side[1] || side[1]+side[2]<=side[0]){return area;} 
 
	//利用海伦公式。s=sqr(p*(p-a)(p-b)(p-c)); 
	double p = (side[0]+side[1]+side[2])/2; //半周长;
	area = sqrt(p*(p-side[0])*(p-side[1])*(p-side[2])); 
	
	return area;
}

inline int planarInterpolation(Vector3D va, Vector3D vb, Vector3D vc, 
                                Vector3D pa, Vector3D pb, Vector3D pc, 
                                Vector3D queryvec, Vector3D& point, int tetId)
{
    // Vector3D vab = vb-va;
    // Vector3D vac = vc-va;
    // Vector3D vbc = vc-vb;

    double Sabc = count_triangle_area(va,vb,vc);
    double Sbcv = count_triangle_area(vb,vc,queryvec);
    double Sacv = count_triangle_area(va,vc,queryvec);
    double Sabv = count_triangle_area(va,vb,queryvec);
    // if(tetId ==  1377834) std::cout  << "--1-0 " << Sabc << " " << Sbcv << " " << Sacv << " " << Sabv  << std::endl;
    // if(Sabc == -1 || Sbcv == -1 || Sacv == -1 || Sabv == -1){return false;}

    if(Sabc == -1){return 0;}
    if(Sbcv == -1){lineInter(vb,vc,pb,pc,queryvec,point); return 31; }
    if(Sacv == -1){lineInter(va,vc,pa,pc,queryvec,point); return 32; }
    if(Sabv == -1){lineInter(va,vb,pa,pb,queryvec,point); return 33; }

    // double Sabc = lenth(cross(vab,vac)); 
    // if(Sabc < EPS || Sabc > -EPS ) {/*std::cout << "--1*1:"<< Sabc << std::endl;*/return false;}
    // double Sbcv = lenth(cross(vbc,queryvec - vb)); 
    // double Sacv = lenth(cross(vac,queryvec - va)); 
    // double Sabv = lenth(cross(vab,queryvec - va));
    // double Sabc = Sbcv + Sacv + Sabv;

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
    // 预备知识：
    // Wa  = Sbcp / Sabc
    // Wb  = Sacp / Sabc
    // Wc  = Sabp / Sabc
    // Sabc = cross(Vac, Vab) / 2 
    // Ep = Wa*Ea + Wb*Eb + Wc*Ec
    // 1 = Wa + Wb + Wc
}

inline int planarInterpolation2(Vector3D va, Vector3D vb, Vector3D vc, 
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
    // 预备知识：
    // Wa  = Sbcp / Sabc
    // Wb  = Sacp / Sabc
    // Wc  = Sabp / Sabc
    // Sabc = cross(Vac, Vab) / 2 
    // Ep = Wa*Ea + Wb*Eb + Wc*Ec
    // 1 = Wa + Wb + Wc
}


bool CpuCPFinder::findCriticalPoint2(const PointVelocity3D pv[4], 
                                     Vector3D queryvec,
                                     const int  tetId )
{
    Vector3D point = {NAN,NAN,NAN};
    // 第一个点表示这个区域角的数量，---》四面体
    double barycoords[4];
   
    Vector3D cellvecters[4] = { // 读取四面体的四个顶点
        pv[0].getVelocity(), 
        pv[1].getVelocity(), 
        pv[2].getVelocity(), 
        pv[3].getVelocity()
    };
    // if(tetId ==  1377834) std::cout << "+0" << std::endl;
    // critical point is vertex
    // 顶点查找是否有可用的
    bool isVert = false;
    for(int i=0;i<4;i++){
        if( cellvecters[i].x > -EPS && cellvecters[i].x < EPS && 
            cellvecters[i].y > -EPS && cellvecters[i].y < EPS && 
            cellvecters[i].z > -EPS && cellvecters[i].z < EPS){
            point = pv[i].getPoint();
            criticalPoints.push_back(PointTet(point,tetId,4));
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
    if(xMax < 0 || xMin > 0){return false;}
    if(yMax < 0 || yMin > 0){return false;}
    if(zMax < 0 || zMin > 0){return false;}

    // if(tetId ==  1377834) std::cout << "+1" << std::endl;
    // 判断点是否在四面体内部
    REAL a[3] = {pv[0].getVelocity().x, pv[0].getVelocity().y, pv[0].getVelocity().z};
    REAL b[3] = {pv[1].getVelocity().x, pv[1].getVelocity().y, pv[1].getVelocity().z};
    REAL c[3] = {pv[2].getVelocity().x, pv[2].getVelocity().y, pv[2].getVelocity().z};  
    REAL d[3] = {pv[3].getVelocity().x, pv[3].getVelocity().y, pv[3].getVelocity().z};
    REAL p[3] = {queryvec.x, queryvec.y, queryvec.z};
    // if(tetId ==  1377834) std::cout << "+2" << std::endl;
    REAL F1 = orient3d(a,b,c,d);
    REAL F2 = orient3d(a,b,d,c);
    REAL F3 = orient3d(a,c,d,b);
    REAL F4 = orient3d(b,c,d,a);

    REAL Fo1 = orient3d(a,b,c,p);
    REAL Fo2 = orient3d(a,b,d,p);
    REAL Fo3 = orient3d(a,c,d,p);
    REAL Fo4 = orient3d(b,c,d,p);
    // if(tetId ==  1377834) std::cout << "+3"  << F1 << Fo1 << std::endl << F2 << Fo2 << std::endl << F3 << Fo3 << std::endl << F4 << Fo4 << std::endl;
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
            // if(tetId ==  1377834) std::cout << "+3+1" << std::endl;
            // CoPoint abcd
            if(cellvecters[0] == cellvecters[1] &&  cellvecters[0] == cellvecters[2] && cellvecters[0] == cellvecters[3])
                { return false; }

            // if(tetId ==  1377834) std::cout << "+3+2" << std::endl;
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
                    criticalPoints.push_back(PointTet(point,tetId,300));
                    return true;
                }else{ return false;}
            }
            // if(tetId ==  1377834) std::cout << "+3+3" << std::endl;

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


            // if(tetId ==  1377834) std::cout << tp[0] << " " << tp[1] << " " << tp[2] << " " << tp[3] << std::endl;

            REAL minLenth = 1000000000.0;
            int nearId = -1;
            for(int ti=0;ti<4;ti++){
                if(tp[ti] == 1){
                    interpolation(pv,p[ti],v[ti]);
                    // if(tetId ==  1377834)  std::cout << v[ti].x << " " << v[ti].y << " " << v[ti].z << std::endl;
                    REAL len = lenth(v[ti]);
                    if(len < minLenth){
                        minLenth = len;
                        nearId = ti;
                    }
                }
            }

            if(nearId == -1) {return false;}
            criticalPoints.push_back(PointTet(p[nearId],tetId,200));
            return true;
        }
        return false;
    }

    // if(tetId ==  1377834) std::cout << "+4" << std::endl;
    // 请求点与四面体是否共面
    if(Fo1 == 0){
        int type = 21;
        // if(tetId ==  1377834) std::cout << "-1:" << Fo1 << std::endl;
        type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 21;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[2].getVelocity(),
                                  pv[0].getPoint(),pv[1].getPoint(),pv[2].getPoint(),
                                  queryvec, point, tetId);
            // if(type == 0) return false;
        }
        criticalPoints.push_back(PointTet(point,tetId,type));
        return true;
    }else if(Fo2 == 0){
        int type = 22;
        // if(tetId ==  1377834) std::cout << "-2:" << Fo2 << std::endl;
        type = planarInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 22;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[1].getVelocity(),pv[3].getVelocity(),
                                  pv[0].getPoint(),pv[1].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
            // if(type == 0) return false;
        }
        criticalPoints.push_back(PointTet(point,tetId,type));
        return true;
    }else if(Fo3 == 0){
        int type = 23;
        // if(tetId ==  1377834) std::cout << "-3:" << Fo3 << std::endl;
        type = planarInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 23;}
        else if(type == 0){
            type = lineInterpolation(pv[0].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                  pv[0].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
            // if(type == 0) return false;

        }
        // if(tetId ==  1377834) std::cout << "-3:" << point.x<< ","<<point.y <<"," <<point.z<< std::endl;
        criticalPoints.push_back(PointTet(point,tetId,type));
        return true; 
    }else if(Fo4 == 0){
        int type = 24;
        // if(tetId ==  1377834) std::cout << "-4" << std::endl;
        type = planarInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                queryvec, point, tetId);
        if(type == 1){type = 24;}
        else if(type == 0){
            type = lineInterpolation(pv[1].getVelocity(),pv[2].getVelocity(),pv[3].getVelocity(),
                                  pv[1].getPoint(),pv[2].getPoint(),pv[3].getPoint(),
                                  queryvec, point, tetId);
            // if(type == 0) return false;
        }
        criticalPoints.push_back(PointTet(point,tetId,type)); 
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
    // if(tetId ==  1377834) std::cout << "1" << std::endl;

    double vol = dot(v2v0, cross(v1v0, v3v2));
    double tetvolumeinv = 1.0f / vol; // 1 / 6Vabcd
    // if(tetId ==  1377834) std::cout << "2" << std::endl;

    // calculate the barycentric coordinates
    barycoords[0] = dot(cellvecters[2] - queryvec, cross(cellvecters[1] - queryvec, v3v2)) * tetvolumeinv;
    // if(tetId ==  1377834) std::cout << "2-2 tetvolumeinv:" << vol << " barycoords[0]:" << barycoords[0] << std::endl;
    // if(barycoords[0] == 0 || barycoords[0] == 1) {return checkRefine(pv, queryvec, tetId, criticalPoints);}
    // if(barycoords[0] < 0 || barycoords[0] > 1) {return false;}
    // if(tetId ==  1377834) std::cout << "3" << std::endl;
    
    barycoords[1] = dot(v2v0, cross(pv0, v3v2)) * tetvolumeinv;
    // if(tetId ==  1377834) std::cout << "3-2 tetvolumeinv:" << vol << " barycoords[1]:" << barycoords[1] << std::endl;
    // if(barycoords[1] == 0 || barycoords[1] == 1) {return checkRefine(pv, queryvec, tetId, criticalPoints);}
    // if(barycoords[1] < 0 || barycoords[1] > 1) {return false;}

    // if(tetId ==  1377834) std::cout << "4" << std::endl;
    barycoords[2] = dot(pv0, cross(v1v0, cellvecters[3] - queryvec)) * tetvolumeinv;
    // if(tetId ==  1377834) std::cout << "4-2" << " barycoords[2]:" << barycoords[2]  << std::endl;
    // if(barycoords[2] == 0 || barycoords[2] == 1) {return checkRefine(pv, queryvec, tetId, criticalPoints);}
    // if(barycoords[2] < 0 || barycoords[2] > 1) {return false;}

    barycoords[3] = 1.0f - barycoords[0] - barycoords[1] - barycoords[2];
    if(barycoords[3] < 0) barycoords[3] = 0;
    // barycoords[3] = dot(v1v0, cross(cellvecters[2]-cellvecters[0], cellvecters[1] - queryvec)) * tetvolumeinv;
    // if(tetId ==  1377834) std::cout << "5"  << " barycoords[3]:" << barycoords[3] << std::endl;
    // if(barycoords[3] == 0 || barycoords[3] == 1) {return checkRefine(pv, queryvec, tetId, criticalPoints);}
    // if(barycoords[3] < 0 || barycoords[3] > 1) {return false;}
    
    // if(tetId ==  1377834) std::cout << "6" << std::endl;
    // compute barycentric interpolation
    point = pv[0].getPoint() * barycoords[0] +
            pv[1].getPoint() * barycoords[1] +
            pv[2].getPoint() * barycoords[2] +
            pv[3].getPoint() * barycoords[3];
    
    // if(isNaN(fabs(point.x)) || isNaN(fabs(point.y)) || isNaN(fabs(point.z))){return false;}
    if(std::isinf(point.x) || std::isinf(point.y) || std::isinf(point.z) || 
        std::isnan(point.x) || std::isnan(point.y) || std::isnan(point.z)){return false;}

    // if(tetId ==  1377834) std::cout << "7" << std::endl;
    criticalPoints.push_back(PointTet(point,tetId,1)); 
    return true;
}

CpuCPFinder::CpuCPFinder(const char* filename) : BaseCPFinder(filename)
{
    initialization();
}


void CpuCPFinder::search()
{
    Vector3D velocity = {0.0,0.0,0.0};
    for(int t=0; t < seqtets.size()/4; t++){
        if(findCriticalPoint2(&seqtets[t*4], velocity, t)){
            ;
            //std::cout << "Tet [" << t << "] have critical points." << std::endl;
        }
    }
    std::cout << criticalPoints.size() << std::endl;
}  

void CpuCPFinder::initialization()
{
    exactinit();

}




