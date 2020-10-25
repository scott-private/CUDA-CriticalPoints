/*
============================================================================
Author      : Scott Fu
Date        : 17/06/2020
Copyright   : scottfu@foxmail.com
File Name   : main.cu
============================================================================
*/

#include <time.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include "common.h"
#include "CriticalPoints.h"

clock_t start_time_;

void start_timer(const std::string& message = "")
{
    if (!message.empty()) {
        printf("%s", message.c_str());
        fflush(stdout);
    }
    start_time_ = clock();
}

double stop_timer()
{
    return double(clock()-start_time_)/CLOCKS_PER_SEC;
}

int main (int argc, char * argv[])
{
    if(argc < 2) {
        std::cout << "ERROR: No .vtk file path specified!" << std::endl; 
		std::cout << "./CriticalPoints [Vtk Filename] [Default:Use Cpu|1:Use Gpu]" << std::endl;
        return 1;
    }

	bool useGpu = false;
	BaseCPFinder* finder = nullptr;
	if(argc == 2 ){std::cout << "-------------- CPU RUN --------------" << std::endl;}
	else{  useGpu = true; std::cout << "-------------- GPU RUN --------------" << std::endl; }
	
	// 初始化
	start_timer("Init:\n");
	if(useGpu){ finder = new GpuSSCPFinder(argv[1]);}
	else{ finder = new CpuCPFinder(argv[1]);}
	printf("Done (%g seconds)\n", stop_timer());
	
	// 获取关键点
	start_timer("Search:\n");
	finder->search();
	printf("Done (%g seconds)\n", stop_timer());

	// 去除重复
	start_timer("Unique:\n");
	finder->sortUnique();
	printf("Done (%g seconds)\n", stop_timer());
	// flannUnique(criticalPoints);
	
	// 关键点分类 
	start_timer("Classification:\n");
	finder->classification(Vector3D(0.01,0.01,0.01));
	printf("Done (%g seconds)\n", stop_timer());


	// 输出
	start_timer("Outfile:\n");
	std::string outfile(argv[1]);
	finder->outfile(outfile + "out.txt");
	printf("Done (%g seconds)\n", stop_timer());
		
	std::cout << "-------------------------------------" << std::endl;
    std::cout << "Number of Critical Tetrahedron: " << finder->criticalTets.size() << std::endl;
    std::cout << "Number of Critical Points: " << finder->criticalPoints.size() << std::endl;
    std::cout << "-------------------------------------" << std::endl;

	// 销毁
	delete finder;

    return 0;
}
