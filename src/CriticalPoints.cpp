/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : CriticalPoints.cpp
============================================================================
*/
#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include "CriticalPoints.h"

BaseCPFinder::BaseCPFinder(const char* filename)
{
    if(loadvtk(filename)) {
        std::cout << "ERROR: load VTK error!" << std::endl;
        return;
    }
}

int BaseCPFinder::loadvtk(const char* filename)
{
    std::vector<float> points;
    std::vector< std::vector<int> > cells;
    std::vector<float> vectors;
    unsigned int pnum;
	std::ifstream infile;
	infile.open(filename);
	std::string line;
	while (getline(infile, line)){
		std::string first; 	// 将字符串读到input中 
		std::stringstream input;
		input << line;  // 依次输出到result中
		input >> first;
		if(first == "#" || first ==" ") continue;
		else if(first == "DATASET" || first == "dataset"){
			std::string second;
			input >> second;
			if(second == "UNSTRUCTURED_GRID" || second == "unstructured_grid" ){
				continue;
			}else{
				std::cout << "ERROR: Only support unstructured grid" << std::endl;
				return 1;
			}
		}
		else if(first == "BINARY" || first == "binary"){
			std::cout << "ERROR: Not support vtk binary" << std::endl;
			return 1;
		}
		else if(first == "POINTS" || first == "points")
        {
			input >> pnum;		
			while (getline(infile, line)){
				std::stringstream tinput;
				tinput << line;
				float tmp;
				while(tinput >> tmp){ points.push_back(tmp); }
				if(points.size() >= pnum *3) break;
			}
		}
		else if(first == "CELLS" || first == "cells")
        {
			int lnum;
			input >> lnum;
			while (getline(infile, line))
            {
				std::stringstream tinput;
				tinput << line;
				int tmp,trinum;
				tinput >> trinum;
                std::vector<int> cell;
				while(trinum > cell.size()){
					tinput >> tmp;
                    cell.push_back(tmp);
				}
                cells.push_back(cell);
				if(cells.size() >= lnum) break;
			}
		}
        else if(first == "VECTORS" || first == "vectors")
        {	
			while (getline(infile, line)){
				std::stringstream tinput;
				tinput << line;
				float tmp;
				while(tinput >> tmp){ vectors.push_back(tmp); }
				if(vectors.size() >= pnum *3) break;
			}
        }
	}


    for(auto& cell : cells){
        std::vector<PointVelocity3D> tet;
        for(auto& pid : cell){
            PointVelocity3D p = {points[pid*3+0],points[pid*3+1],points[pid*3+2],vectors[pid*3+0],vectors[pid*3+1],vectors[pid*3+2]};
            tet.push_back(p);
            seqtets.push_back(p); 
        }
        tets.push_back(tet);
    }
    return 0;
}

void BaseCPFinder::interpolation(const PointVelocity3D pv[4], Vector3D point, Vector3D& velocity)
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

void BaseCPFinder::sortUnique()
{
	sort(criticalPoints.begin(),criticalPoints.end());
	criticalPoints.erase(unique(criticalPoints.begin(), criticalPoints.end()), criticalPoints.end());

	criticalTets.resize(0);
	for(auto& c : criticalPoints){criticalTets.push_back(c.tetId);}

	sort(criticalTets.begin(),criticalTets.end());
	criticalTets.erase(unique(criticalTets.begin(), criticalTets.end()), criticalTets.end());
}


void BaseCPFinder::outfile(std::string outfile,bool check)
{
	std::ofstream f;
	f.open(outfile);
	if(f) {
		f << "-------------------------------------" << std::endl;
		f << "Number of Critical Tetrahedron: " << criticalTets.size() << std::endl;
    	f << "Number of Critical Points: " << criticalPoints.size() << std::endl;
		f << "-------------------------------------" << std::endl;
		for(auto tet : criticalTets) {f << tet <<  std::endl;} 
		f << "-------------------------------------" << std::endl;
        if(check){
            for(auto pt : criticalPoints) {
                Vector3D velocity={NAN,NAN,NAN};
                interpolation(tets[pt.tetId].data(),pt.pos,velocity);
                f << "tet:" << pt.tetId << " pos:" << pt.pos.x << " " << pt.pos.y << " " << pt.pos.z \
                    << " v( " << velocity.x << " " << velocity.y << " " << velocity.z << " ) type:" << pt.type << std::endl;
		    }
        }else{
            for(auto pt : criticalPoints) { pt.out2file(f);}
        }
	}
	f.close();
}
