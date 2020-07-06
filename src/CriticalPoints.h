/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : CriticalPoints.h
============================================================================
*/

#ifndef CRITICAL_POINTS_H
#define CRITICAL_POINTS_H

#include <cstddef> // std::size_t
#include <vector>
#include <vector>
#include <iostream>
#include <fstream>
#include "common.h"
#include "Geometry3D.h"

struct PointTet{
    Vector3D pos;
    int      tetId;
    int      type;

    __host__ __device__
    PointTet(){pos.x  = 0; pos.y  = 0; pos.z  = 0; tetId = 0; type = 0;} 

    __host__ __device__
    PointTet(REAL _x,REAL _y,REAL _z, int _id, int _type=1){ pos.x  = _x; pos.y  = _y; pos.z  = _z; tetId = _id; type = _type;} 

    __host__ __device__
    PointTet(Vector3D _p, int _id, int _type=1){ pos.x  = _p.x; pos.y  = _p.y; pos.z  = _p.z; tetId = _id;  type = _type;} 

    __host__ __device__
    bool operator==(PointTet pt){
        if( (pt.pos.x - pos.x < EPS && pt.pos.x - pos.x > -EPS ) && 
            (pt.pos.y - pos.y < EPS && pt.pos.y - pos.y > -EPS ) &&
            (pt.pos.z - pos.z < EPS && pt.pos.z - pos.z > -EPS ) ){
            return true;
        }
        return false;
    }

    __host__ __device__
    bool operator<(PointTet pt){
        if(pos.x < pt.pos.x ) return true;
        else if (pos.x > pt.pos.x ) return false;
        else{
            if(pos.y < pt.pos.y ) return true;
            else if (pos.y > pt.pos.y ) return false;
            else{
                if(pos.z < pt.pos.z ) return true;
                else if (pos.z > pt.pos.z ) return false;
                else return false;
            }
        }
    }

    void print() const{
        std::cout << "tet:" << tetId << " pos:" << pos.x << " " << pos.y << " " << pos.z << " type:" << type << std::endl;
    }

    void out2file(std::ofstream& f) const{
        f << "tet:" << tetId << " pos:" << pos.x << " " << pos.y << " " << pos.z  << " type:" << type << std::endl;
    }

};

struct PointVelocity3D {
    // REAL x, y, z, vx, vy, vz;
    REAL dat[6];
    int tetId;

    __host__ __device__
    REAL operator[](int i) const {return dat[i];}

    __host__ __device__
    REAL x() const {return dat[0];}
    __host__ __device__
    REAL y() const {return dat[1];}
    __host__ __device__
    REAL z() const {return dat[2];}
    
    __host__ __device__
    REAL vx() const {return dat[3];}
    __host__ __device__
    REAL vy() const {return dat[4];}
    __host__ __device__
    REAL vz() const {return dat[5];}

    __host__ __device__
    Vector3D getPoint() const {return Vector3D(dat[0],dat[1],dat[2]);}
    __host__ __device__
    Vector3D getVelocity() const {return Vector3D(dat[3],dat[4],dat[5]);}

    void print() const{
        std::cout << " x:" << x() << " y:" << y() << " z:" << z() \
            << " vx:" << vx() << " vy:" << vy() << " vz:" << vz() << std::endl;
    }
};


class BaseCPFinder{
private:
    long long getFileSize(const char* filename);
    int loadvtk(const char* filename);
    int loadWriteFile(const char* filename, const std::string outPoints,
					  const std::string outVelocities,const std::string outCells);
    std::vector<Vector3D> lookupPoints(std::vector<unsigned int> pvec);
    std::vector<Vector3D> lookupVelocities(std::vector<unsigned int> vvec);

    Vector3D lookupPoints(const unsigned int off);
    Vector3D lookupVelocities(const unsigned int off);

    unsigned int cellsNum;
    unsigned int pointsNum;
    unsigned int curCellIdx;
    unsigned int tetsNumPerBlock;
public:
    BaseCPFinder(const char* filename);
    virtual ~BaseCPFinder() {};
    
    virtual void search() = 0;
    virtual void search1(){};

    virtual void sortUnique();

    // virtual  void check() = 0;
    void outfile(std::string outfile, bool check=false);

    std::vector<int> criticalTets;
    std::vector<PointTet> criticalPoints;
protected:
    virtual void initialization() = 0;
    void interpolation(const PointVelocity3D pv[4], Vector3D point, Vector3D& velocity);
    int loadDataBlock();

    bool isBigFile;
	std::string outCells;
    std::string outPoints;
	std::string outVelocities;
    std::vector<PointVelocity3D> seqtets;
};


// Critical Points Finder with CPU
class CpuCPFinder : public BaseCPFinder{
private:
    bool findCriticalPoint2(const PointVelocity3D pv[4], Vector3D queryvec, const int tetId);
public:
    CpuCPFinder(const char* filename) ;
    
    virtual ~CpuCPFinder(){};
    
    virtual void search();

protected:
    virtual void initialization();
};

// Using GPU  Critical Points Finder
class GpuCPFinder : public BaseCPFinder{
private:
    int tetnum;
    PointVelocity3D* d_tets1;
    Vector3D* d_queryvec;

    PointTet* d_outpoints;
    PointTet* d_outpoints2;
    int* d_outstencil;
    int* d_pointsize;
public:
    GpuCPFinder(const char* filename) ;
    
    virtual ~GpuCPFinder();

    virtual void search();

protected:
    virtual void initialization();
};


// Using GPU critical point finder by Sub-Step
class GpuSSCPFinder : public BaseCPFinder{
private:
    int tetnum;
    PointVelocity3D* d_tets1;
    PointVelocity3D* d_tets2;
    Vector3D* d_queryvec;

    PointTet* d_outpoints;
    int* d_outstencil;
    int* d_pointsize;
public:
    GpuSSCPFinder(const char* filename) ;
    
    virtual ~GpuSSCPFinder();

    virtual void search();

protected:
    virtual void initialization();
};

#endif
