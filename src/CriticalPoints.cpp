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
	isBigFile = false;
	long long byteSize = getFileSize(filename);
	if(byteSize / 1024 / 1024 / 1024 > 0){
		// Batch storage loading is required
		isBigFile = true;
		outCells = "outCPFinderCells.cpbin";
		outPoints = "outCPFinderPoints.cpbin";
		outVelocities = "outCPFinderVelocities.cpbin";
     	curCellIdx = 0;
    	tetsNumPerBlock = 5120000;

		if(loadWriteFile(filename,outPoints,outVelocities,outCells)) {
			std::cout << "ERROR: load VTK error!" << std::endl;
			exit(0);
		}
		
		std::cout << "    cellsNum " << cellsNum << std::endl;
		std::cout << "    pointsNum " << pointsNum << std::endl;
		std::cout << "    tetsNumPerBlock " << tetsNumPerBlock << std::endl;

		if(loadDataBlock()) {
			std::cout << "ERROR: load block data error!" << std::endl;
			exit(0);
		}
	}else {
		if(loadvtk(filename)) {
        	std::cout << "ERROR: load VTK error!" << std::endl;
        	exit(0);
		}
    }
}

//假设文件file.txt存在,且在当前目录下
long long BaseCPFinder::getFileSize(const char* filename)
{
	std::ifstream in(filename);
	in.seekg(0, std::ios::end); //设置文件指针到文件流的尾部
	std::streampos ps = in.tellg(); //读取文件指针的位置
	std::cout << "    file size: " << ps << " byte" << std::endl;
	in.close(); //关闭文件流
	return ps;
}

int BaseCPFinder::loadvtk(const char* filename)
{
    std::vector<float> points;
    std::vector< std::vector<int> > cells;
    std::vector<float> vectors;
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
			std::cout << "    load points ... " << std::flush;
			input >> pointsNum;		
			while (getline(infile, line)){
				std::stringstream tinput;
				tinput << line;
				float tmp;
				while(tinput >> tmp){ points.push_back(tmp); }
				if(points.size() >= pointsNum *3) break;
			}
			std::cout << "end" << std::endl;
		}
		else if(first == "CELLS" || first == "cells")
        {
			std::cout << "    load cells ... " << std::flush;
			input >> cellsNum;
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
				if(cells.size() >= cellsNum) break;
			}
			std::cout << "end" << std::endl;
		}
        else if(first == "VECTORS" || first == "vectors")
        {	
			std::cout << "    load vectors ..." << std::flush;
			while (getline(infile, line)){
				std::stringstream tinput;
				tinput << line;
				float tmp;
				while(tinput >> tmp){ vectors.push_back(tmp); }
				if(vectors.size() >= pointsNum *3) break;
			}
			std::cout << "end" << std::endl;
        }
	}
	infile.close();

	std::cout << "    change data ... " << std::flush; 
    for(auto& cell : cells){
		// std::vector<PointVelocity3D> tet;
        for(auto& pid : cell){
            PointVelocity3D p = {points[pid*3+0],points[pid*3+1],points[pid*3+2],vectors[pid*3+0],vectors[pid*3+1],vectors[pid*3+2]};
			p.tetId = seqtets.size()/4;
            seqtets.push_back(p); 
			// tet.push_back(p); 
        }
		// tets.push_back(tet);
    }
	std::cout << "end" << std::endl;
    return 0;
}


int BaseCPFinder::loadWriteFile(const char* filename,
								const std::string outPoints,
								const std::string outVelocities,
								const std::string outCells)
{
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
			std::cout << "    load points ... " << std::flush;
			std::ofstream ofile(outPoints, std::ios::out | std::ios::binary);
			unsigned int linenum = 0;
			input >> pointsNum;		
			while (getline(infile, line)){
				std::stringstream tinput(line);
				Vector3D point;
				tinput >> point.x;
				tinput >> point.y;
				tinput >> point.z;
				ofile.write(reinterpret_cast<char *>(&point),sizeof(Vector3D));
				linenum++;
				if(linenum >= pointsNum) break;
			}
			ofile.close();
			std::cout << "end" << std::endl;
		}
		else if(first == "CELLS" || first == "cells")
        {
			std::cout << "    load cells ... " << std::flush;
			std::ofstream ofile(outCells, std::ios::out | std::ios::binary);
			input >> cellsNum;
			unsigned int linenum = 0;
			while (getline(infile, line)){
				std::stringstream tinput(line);
				unsigned int a,b,c,d;
				tinput >> a;
				tinput >> a;
				tinput >> b;
				tinput >> c;
				tinput >> d;
				ofile.write(reinterpret_cast<char *>(&a),sizeof(unsigned int));
				ofile.write(reinterpret_cast<char *>(&b),sizeof(unsigned int));
				ofile.write(reinterpret_cast<char *>(&c),sizeof(unsigned int));
				ofile.write(reinterpret_cast<char *>(&d),sizeof(unsigned int));
				linenum++;
				if(linenum >= cellsNum) break;
			}
			ofile.close();
			std::cout << "end" << std::endl;
		}
        else if(first == "VECTORS" || first == "vectors")
        {	
			std::cout << "    load vectors ..." << std::flush;
			std::ofstream ofile(outVelocities, std::ios::out | std::ios::binary);
			unsigned int linenum = 0;
			while (getline(infile, line)){
				std::stringstream tinput(line);
				Vector3D velocity;
				tinput >> velocity.x;
				tinput >> velocity.y;
				tinput >> velocity.z;
				ofile.write((char*)(&velocity),sizeof(Vector3D));
				linenum++;
				if(linenum >= pointsNum) break;
			}
			ofile.close();
			std::cout << "end" << std::endl;
        }
	}
	infile.close();
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
        // if(check){
        //     for(auto pt : criticalPoints) {
        //         Vector3D velocity={NAN,NAN,NAN};
        //         interpolation(tets[pt.tetId].data(),pt.pos,velocity);
        //         f << "tet:" << pt.tetId << " pos:" << pt.pos.x << " " << pt.pos.y << " " << pt.pos.z \
        //             << " v( " << velocity.x << " " << velocity.y << " " << velocity.z << " ) type:" << pt.type << std::endl;
		//     }
        // }else{
            for(auto pt : criticalPoints) { pt.out2file(f);}
        // }
	}
	f.close();
}

// int BaseCPFinder::loadDataBlock1()
// {
// 	if(curCellIdx >= cellsNum) 
// 		{return 1;}
// 	seqtets.resize(0);
// 	std::ifstream infile(outCells, std::ios::in | std::ios::binary);
// 	unsigned int offline = curCellIdx * 4 * sizeof(unsigned int) ;
// 	infile.seekg(offline,std::ios::beg);
// 	unsigned int a,b,c,d;
// 	for(unsigned int i=0;i<tetsNumPerBlock;i++){
// 		if(curCellIdx >= cellsNum) {break;}
// 		infile.read(reinterpret_cast<char *>(&a),sizeof(unsigned int));
// 		infile.read(reinterpret_cast<char *>(&b),sizeof(unsigned int));
// 		infile.read(reinterpret_cast<char *>(&c),sizeof(unsigned int));
// 		infile.read(reinterpret_cast<char *>(&d),sizeof(unsigned int));			
// 		Vector3D p = lookupPoints(a);
// 		Vector3D v = lookupVelocities(a);
// 		PointVelocity3D pv1 = {p.x, p.y, p.z, v.x, v.y, v.z};
// 		pv1.tetId = curCellIdx;
// 		p = lookupPoints(b);
// 		v = lookupVelocities(b);
// 		PointVelocity3D pv2 = {p.x, p.y, p.z, v.x, v.y, v.z};
// 		pv2.tetId = curCellIdx;
// 		p = lookupPoints(c);
// 		v = lookupVelocities(c);
// 		PointVelocity3D pv3 = {p.x, p.y, p.z, v.x, v.y, v.z};
// 		pv3.tetId = curCellIdx;
// 		p = lookupPoints(d);
// 		v = lookupVelocities(d);
// 		PointVelocity3D pv4 = {p.x, p.y, p.z, v.x, v.y, v.z};
// 		pv4.tetId = curCellIdx;
// 		seqtets.push_back(pv1);
// 		seqtets.push_back(pv2);
// 		seqtets.push_back(pv3);
// 		seqtets.push_back(pv4);
// 		curCellIdx++;
// 	}
// 	infile.close();
// 	return 0;
// }


int BaseCPFinder::loadDataBlock()
{
	if(curCellIdx >= cellsNum) 
		{return 1;}
	seqtets.resize(0);
	std::ifstream infile(outCells, std::ios::in | std::ios::binary);
	unsigned int offline = curCellIdx * 4 * sizeof(unsigned int) ;
	infile.seekg(offline,std::ios::beg);

	std::vector<unsigned int> tetsBlock(tetsNumPerBlock*4);
	if(curCellIdx+tetsNumPerBlock >= cellsNum) 
	{tetsBlock.resize((cellsNum-curCellIdx)*4);}
	infile.read(reinterpret_cast<char *>(tetsBlock.data()),tetsBlock.size()*sizeof(unsigned int));
	std::vector<Vector3D> points = lookupPoints(tetsBlock);
	std::vector<Vector3D> velocities = lookupVelocities(tetsBlock);
	for(unsigned int i=0; i < tetsBlock.size()/4; i++){
		PointVelocity3D pv1 = {points[i*4].x, points[i*4].y, points[i*4].z,
							   velocities[i*4].x, velocities[i*4].y, velocities[i*4].z};
		PointVelocity3D pv2 = {points[i*4+1].x, points[i*4+1].y, points[i*4+1].z, 
							   velocities[i*4+1].x, velocities[i*4+1].y, velocities[i*4+1].z};
		PointVelocity3D pv3 = {points[i*4+2].x, points[i*4+2].y, points[i*4+2].z, 
							   velocities[i*4+2].x, velocities[i*4+2].y, velocities[i*4+2].z};
		PointVelocity3D pv4 = {points[i*4+3].x, points[i*4+3].y, points[i*4+3].z, 
							   velocities[i*4+3].x, velocities[i*4+3].y, velocities[i*4+3].z};
		pv1.tetId = curCellIdx;
		pv2.tetId = curCellIdx;
		pv3.tetId = curCellIdx;
		pv4.tetId = curCellIdx;
		seqtets.push_back(pv1);
		seqtets.push_back(pv2);
		seqtets.push_back(pv3);
		seqtets.push_back(pv4);
		curCellIdx++;
	}
	infile.close();
	return 0;
}


Vector3D BaseCPFinder::lookupPoints(const unsigned int off)
{
	std::ifstream infile;
	infile.open(outPoints, std::ios::out|std::ios::binary);
	infile.seekg(off * sizeof(Vector3D),std::ios::beg);
	Vector3D point;
	infile.read(reinterpret_cast<char *>(&point),sizeof(Vector3D));
	infile.close();
	return point;
}

Vector3D BaseCPFinder::lookupVelocities(const unsigned int off)
{
	std::ifstream infile;
	infile.open(outVelocities, std::ios::out|std::ios::binary);
	infile.seekg(off * sizeof(Vector3D),std::ios::beg);
	Vector3D velocity;
	infile.read(reinterpret_cast<char *>(&velocity),sizeof(Vector3D));
	infile.close();
	return velocity;
}


std::vector<Vector3D> BaseCPFinder::lookupPoints(std::vector<unsigned int> pvec)
{
	std::ifstream infile;
	infile.open(outPoints, std::ios::out|std::ios::binary);
	std::vector<Vector3D> pointset;
	for(auto& off: pvec){
		Vector3D point;
		infile.seekg(off * sizeof(Vector3D),std::ios::beg);
		infile.read(reinterpret_cast<char *>(&point),sizeof(Vector3D));
		pointset.push_back(point);
	}
	infile.close();
	return pointset;
}

std::vector<Vector3D> BaseCPFinder::lookupVelocities(std::vector<unsigned int> vvec)
{
	std::ifstream infile;
	infile.open(outVelocities, std::ios::out|std::ios::binary);
	std::vector<Vector3D> velocityset;
	for(auto& off: vvec){
		Vector3D vel;
		infile.seekg(off * sizeof(Vector3D),std::ios::beg);
		infile.read(reinterpret_cast<char *>(&vel),sizeof(Vector3D));
		velocityset.push_back(vel);
	}
	infile.close();
	return velocityset;
}

