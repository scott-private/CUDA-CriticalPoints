/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : Geometry3D.h
============================================================================
*/

#ifndef GEOMETRY3D_H
#define GEOMETRY3D_H

#include "common.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm>

struct Vector3D{
    REAL x,y,z;
    __host__ __device__ 
    Vector3D(){}

    // __host__ __device__ 
    // Vector3D(Vector3D& v){
    //     x=v.x; y=v.y; z=v.z;
    // }

    __host__ __device__ 
    Vector3D(REAL _x,REAL _y,REAL _z){
        x=_x;y=_y;z=_z;
    }

    __host__ __device__ 
    bool operator==(const Vector3D &pos){
        if( (x - pos.x < EPS && x - pos.x > -EPS ) && 
            (y - pos.y < EPS && y - pos.y > -EPS ) &&
            (z - pos.z < EPS && z - pos.z > -EPS ) ){
            return true;
        }
        return false;
    }

    __host__ __device__ 
    Vector3D operator-(const Vector3D &ne) {
        return Vector3D(x-ne.x,y-ne.y,z-ne.z);
    }

    __host__ __device__ 
    Vector3D operator+(const Vector3D &ne) {
        return Vector3D(x+ne.x,y+ne.y,z+ne.z);
    }

    __host__ __device__ 
    Vector3D operator*(const REAL t) {
        return Vector3D(x*t,y*t,z*t);
    }
};

struct Line3D{
    Vector3D a,b;
    __host__ __device__ 
    Line3D(){}

    __host__ __device__ 
    Line3D(Vector3D _a,Vector3D _b){
        a=_a;
        b=_b;
    }

};

struct Plane3D{
    Vector3D a,b,c;
    __host__ __device__ 
    Plane3D(){}

    __host__ __device__ 
    Plane3D(Vector3D _a,Vector3D _b,Vector3D _c){
        a=_a;
        b=_b;
        c=_c;
    }

};

__host__ __device__ 
int sign(REAL x);

__host__ __device__ 
void normalize(Vector3D& v);

__host__ __device__ 
Vector3D cross(Vector3D a,Vector3D b);

__host__ __device__ 
REAL dot(Vector3D a,Vector3D b);

__host__ __device__ 
REAL lenth(Vector3D v);

__host__ __device__ 
REAL lenth2(Vector3D v);

__host__ __device__ 
REAL dist(Vector3D a,Vector3D b);

__host__ __device__ 
REAL dist2(Vector3D a,Vector3D b);

//平面法向量
__host__ __device__ 
Vector3D pvec(Plane3D s);

//判定点是否在线段上，包括端点和共线
__host__ __device__ 
bool point_on_seg(Vector3D p,Line3D s);

//判断点在平面上
__host__ __device__ 
bool point_on_plane(Vector3D p,Plane3D s);

//判断点是否在平面的上方
__host__ __device__ 
bool point_above_plane( Vector3D p, Plane3D s);

//判定点是否在空间三角形上，包括边界，三点共线无意义
__host__ __device__ 
bool point_in_triangle( Vector3D p, Plane3D s );

//判定点是否在空间三角形上，不包括边界，三点共线无意义
__host__ __device__ 
int point_in_triangle2( Vector3D p, Plane3D s ) ;

//判定两点在线段同侧，点在线段上返回0，不共面无意义
__host__ __device__ 
bool same_side( Vector3D p1, Vector3D p2, Line3D l ) ;

//判定两点在线段异侧，点在平面上返回0
__host__ __device__ 
bool opposite_side( Vector3D p1, Vector3D p2, Line3D l );

//判定两点在平面同侧，点在平面上返回0
__host__ __device__ 
bool same_side( Vector3D p1, Vector3D p2, Plane3D s ) ;

//判定两点在平面异侧，点在平面上返回0
__host__ __device__ 
bool opposite_side( Vector3D p1, Vector3D p2, Plane3D s ) ;

//判断直线平行
__host__ __device__ 
bool parallel(Line3D u,Line3D v);

//判定两线段相交，不包括端点和部分重合
__host__ __device__ 
bool seg_seg_inter( Line3D u, Line3D v );

//判定线段与空间三角形相交，包括交于边界和（部分）包含
__host__ __device__ 
int seg_triangle_inter( Line3D l, Plane3D s ) ;

//判定线段与空间三角形相交，不包括交于边界和（部分）包含
__host__ __device__ 
int seg_triangle_inter2( Line3D l, Plane3D s );

//面面平行
__host__ __device__ 
bool parallel(Plane3D s1,Plane3D s2);

//判断直线垂直
__host__ __device__ 
bool vertical(Line3D u,Line3D v);

//面面垂直
__host__ __device__ 
bool vertical(Plane3D s1,Plane3D s2);

//判断两直线的位置关系
__host__ __device__ 
int line_to_line(Line3D u,Line3D v);

//直线与平面关系
__host__ __device__ 
int line_to_plane(Line3D u,Plane3D s);

//线面求交
__host__ __device__ 
Vector3D line_plane_intersection(Line3D u,Plane3D s);

//线线求交
__host__ __device__ 
Vector3D line_interseciton(Line3D u,Line3D v);

//面面求交
__host__ __device__ 
Line3D plane_intersection(Plane3D u,Plane3D v);

//点线距离
__host__ __device__ 
REAL dist_point_to_line(Vector3D p,Line3D u);

//点面距离
__host__ __device__ 
REAL dist_point_to_plane(Vector3D p,Plane3D s);

//线线距离
__host__ __device__ 
REAL dist_line_to_line(Line3D u,Line3D v ) ;

//点线垂足
__host__ __device__ 
Vector3D vertical_foot(Vector3D p,Line3D u);

//已知四面体六边求体积
__host__ __device__ 
REAL volume(REAL a,REAL b,REAL c,REAL d,REAL e,REAL f);

//四面体体积
__host__ __device__ 
REAL volume(Vector3D a,Vector3D b,Vector3D c,Vector3D d);

#endif
