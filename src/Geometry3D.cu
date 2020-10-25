
/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : Geometry3D.cpp
============================================================================
*/

#include "Geometry3D.h"
#include <cassert>

__host__ __device__ 
int sign(REAL x) {
    return x > EPS ? 1 : ( x < -EPS ? (-1) : (0));
}

__host__ __device__ 
void normalize(Vector3D& v){
    REAL l = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
    assert(l > 0);
    v.x/=l;
    v.y/=l;
    v.z/=l;
}

__host__ __device__ 
Vector3D cross(Vector3D a,Vector3D b){
    return Vector3D(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

__host__ __device__ 
REAL dot(Vector3D a,Vector3D b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

__host__ __device__ 
REAL lenth(Vector3D v){
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

__host__ __device__ 
REAL lenth2(Vector3D v){
    return v.x*v.x+v.y*v.y+v.z*v.z;
}

__host__ __device__ 
REAL dist(Vector3D a,Vector3D b){
    return lenth(a-b);
}

__host__ __device__ 
REAL dist2(Vector3D a,Vector3D b){
    return dot(a-b,a-b);
}


// 平面法向量
__host__ __device__ 
Vector3D pvec(Plane3D s){
    return cross(s.b-s.a,s.c-s.a);
}


// 判定点是否在线段上，包括端点和共线
__host__ __device__ 
bool point_on_seg(Vector3D p,Line3D s){
    return sign(lenth(cross(p-s.a,s.b-s.a)))==0&&(p.x-s.a.x)*(p.x-s.b.x)<EPS&&(p.y-s.a.y)*(p.y-s.b.y)<EPS&&(p.z-s.a.z)*(p.z-s.b.y)<EPS;
}


// 判断点在平面上
__host__ __device__ 
bool point_on_plane(Vector3D p,Plane3D s){
    return sign(dot(p-s.a,pvec(s)))==0;
}

// 判断点是否在平面的上方
__host__ __device__ 
bool point_above_plane( Vector3D p, Plane3D s) {
	Vector3D l = p - s.a;
	Vector3D norm = pvec(s);
	REAL proj = dot(norm, l);
	return proj > EPS ? true : false;
}

// 判定点是否在空间三角形上，包括边界，三点共线无意义
__host__ __device__ 
bool point_in_triangle( Vector3D p, Plane3D s ) {
    return sign(lenth(cross(s.a-s.b,s.a-s.c))-lenth(cross(p-s.a,p-s.b))-lenth(cross(p-s.b,p-s.c))-lenth(cross(p-s.c,p-s.a)))!=0;
}


// 判定点是否在空间三角形上，不包括边界，三点共线无意义
__host__ __device__ 
int point_in_triangle2( Vector3D p, Plane3D s ) {
    return point_in_triangle(p,s)&&lenth(cross(p-s.a,p-s.b))>EPS&&lenth(cross(p-s.b,p-s.c))>EPS&&lenth(cross(p-s.c,p-s.a))>EPS;
}


// 判定两点在线段同侧，点在线段上返回0，不共面无意义
__host__ __device__  
bool same_side( Vector3D p1, Vector3D p2, Line3D l ) {
    return dot(cross(l.a-l.b,p1-l.b),cross(l.a-l.b,p2-l.b))>EPS;
}

// 判定两点在线段异侧，点在平面上返回0
__host__ __device__  
bool opposite_side( Vector3D p1, Vector3D p2, Line3D l ) {
    return dot(cross(l.a-l.b,p1-l.b),cross(l.a-l.b,p2-l.b))<-EPS;
}

// 判定两点在平面同侧，点在平面上返回0
__host__ __device__ 
bool same_side( Vector3D p1, Vector3D p2, Plane3D s ) {
    return dot(pvec(s),p1-s.a)*dot(pvec(s),p2-s.a)>EPS;
}

// 判定两点在平面异侧，点在平面上返回0
__host__ __device__ 
bool opposite_side( Vector3D p1, Vector3D p2, Plane3D s ) {
    return dot(pvec(s),p1-s.a)*dot(pvec(s),p2-s.a)<-EPS;
}

// 判断直线平行
__host__ __device__ 
bool parallel(Line3D u,Line3D v){
    return sign(lenth(cross(u.b-u.a,v.b-v.a)))==0;
}

// 判定两线段相交，不包括端点和部分重合
__host__ __device__ 
bool seg_seg_inter( Line3D u, Line3D v ) {
    return point_on_plane(u.a,Plane3D(u.b,v.a,v.b))&&opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}

// 判定线段与空间三角形相交，包括交于边界和（部分）包含
__host__ __device__ 
int seg_triangle_inter( Line3D l, Plane3D s ) {
    return !same_side(l.a,l.b,s)&&!same_side(s.a,s.b,Plane3D(l.a,l.b,s.c))&&!same_side(s.b,s.c,Plane3D(l.a,l.b,s.a))&&!same_side(s.c,s.a,Plane3D(l.a,l.b,s.b));
}

// 判定线段与空间三角形相交，不包括交于边界和（部分）包含
__host__ __device__ 
int seg_triangle_inter2( Line3D l, Plane3D s ){
    return opposite_side( l.a, l.b, s ) && opposite_side( s.a, s.b, Plane3D(l.a, l.b, s.c) ) && opposite_side( s.b, s.c, Plane3D(l.a, l.b, s.a) ) && opposite_side( s.c, s.a,Plane3D(l.a, l.b, s.b) );
}

// 面面平行
__host__ __device__ 
bool parallel(Plane3D s1,Plane3D s2){
    return sign(lenth(cross(pvec(s1),pvec(s2))))==0;
}

// 判断直线垂直
__host__ __device__ 
bool vertical(Line3D u,Line3D v){
    return sign(dot(u.b-u.a,v.b-v.a))==0;
}

// 面面垂直
__host__ __device__ 
bool vertical(Plane3D s1,Plane3D s2){
    return sign(dot(pvec(s1),pvec(s2)))==0;
}

// 判断两直线的位置关系
__host__ __device__ 
int line_to_line(Line3D u,Line3D v){
    Plane3D s1(u.a,u.b,v.a),s2(u.a,u.b,v.b);
    if(sign(lenth(cross(pvec(s1),pvec(s2)))))
        return -1;// 异面
    else if(parallel(u,v))
        return 0;// 平行
    else
        return 1;// 相交
}

// 直线与平面关系
__host__ __device__ 
int line_to_plane(Line3D u,Plane3D s){
    if(sign(dot(pvec(s),u.b-u.a))==0){
        if(point_on_plane(u.a,s))
            return -1;// 直线在平面上
        else
            return 0;// 直线平行于平面
    }

    else
        return 1;// 线面相交
}

// 线面求交
__host__ __device__ 
Vector3D line_plane_intersection(Line3D u,Plane3D s){
    Vector3D ret=pvec(s),der=u.b-u.a;
    REAL t=dot(ret,s.a-u.a)/dot(ret,u.b-u.a);
    return u.a+der*t;
}

// 线线求交
__host__ __device__ 
Vector3D line_interseciton(Line3D u,Line3D v){
    Vector3D ret=u.a,v1=cross(u.a-v.a,v.b-v.a),v2=cross(u.b-u.a,v.b-v.a);
    REAL t=lenth(v1)/lenth(v2)*(dot(v1,v2)>0?-1:1);
    return ret+((u.b-u.a)*t);
}

// 面面求交
__host__ __device__ 
Line3D plane_intersection(Plane3D u,Plane3D v){
    Line3D ret;
    ret.a=(line_to_plane(Line3D(v.a,v.b),u)==0)?line_plane_intersection(Line3D(v.b,v.c),u):line_plane_intersection(Line3D(v.a,v.b),u);
    ret.b=(line_to_plane(Line3D(v.c,v.a),u)==0)?line_plane_intersection(Line3D(v.b,v.c),u):line_plane_intersection(Line3D(v.a,v.c),u);
    return ret;
}

// 点线距离
__host__ __device__ 
REAL dist_point_to_line(Vector3D p,Line3D u){
    return lenth(cross(p-u.a,u.b-u.a))/dist(u.a,u.b);
}

// 点面距离
__host__ __device__ 
REAL dist_point_to_plane(Vector3D p,Plane3D s){
    Vector3D pv=pvec(s);
    return fabs(dot(pv,p-s.a))/lenth(pv);
}

// 线线距离
__host__ __device__ 
REAL dist_line_to_line(Line3D u,Line3D v ) {
    Vector3D p=cross(u.a-u.b,v.a-v.b);
    return fabs(dot(u.a-v.a,p))/lenth(p);
}

// 点线垂足
__host__ __device__ 
Vector3D vertical_foot(Vector3D p,Line3D u){
    REAL t=dot(p-u.a,u.b-u.a)/dist2(u.a,u.b);
    Vector3D ret=u.a;
    return ret+((u.b-u.a)*t);
}

// 已知四面体六边求体积
__host__ __device__ 
REAL volume(REAL a,REAL b,REAL c,REAL d,REAL e,REAL f){
    REAL a2=a*a,b2=b*b,c2=c*c,d2=d*d,e2=e*e,f2=f*f;
    REAL tr1=acos((c2+b2-f2)/(2.0*b*c));
    REAL tr2=acos((a2+c2-e2)/(2.0*a*c));
    REAL tr3=acos((a2+b2-d2)/(2.0*a*b));
    REAL tr4=(tr1+tr2+tr3)/2.0;
    REAL temp=sqrt(sin(tr4)*sin(tr4-tr1)*sin(tr4-tr2)*sin(tr4-tr3));
    return a*b*c*temp/3.0;
}

// 四面体体积
__host__ __device__ 
REAL volume(Vector3D a,Vector3D b,Vector3D c,Vector3D d){
    // abc面方向与d一致时为正
    return fabs(dot(cross(b-a,c-a),d-a))/6.0;
}