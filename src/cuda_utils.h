#pragma once
#include <cuda_runtime.h>
#include <thrust/functional.h>
#include <thrust/remove.h>
#include <cooperative_groups.h>
#include "common.h"

//using namespace cooperative_groups;
namespace cg = cooperative_groups;

// #define FLT_MAX         3.402823466e+38F        /* max value */
// #define FLT_MIN         1.175494351e-38F

__device__ __forceinline__ float cuda_time()
{
    return clock();
}

__device__ __forceinline__ float cuda_time(float tim)
{
    //return ((clock() - tim) / (MAX_GPU_CLOCK_MHZ )); //us
    return ((clock() - tim) / (MAX_GPU_CLOCK_MHZ * 1000)); // ms
    //return ((clock() - tim) / (MAX_GPU_CLOCK_MHZ * 1000 * 1000)); // s
}
 
__device__ __forceinline__ void lock(int *mutex) 
{
    while( atomicCAS(mutex, 0, 1) != 0 );
}

__device__ __forceinline__ void unlock(int *mutex) 
{
    atomicExch(mutex, 0); 
}


template<typename value_t, int tile_sz>
__device__ __forceinline__
value_t reduce_sum(cg::thread_block_tile<tile_sz> g, value_t val)
{
    value_t my_val = val;
    #pragma unroll
    for (int i = g.size() / 2; i > 0; i >>= 1){
        const value_t other_val = g.shfl_down(my_val, i);
        my_val += other_val;
    }
    return my_val; // note: only thread 0 will return full sum
}


template<typename value_t, int tile_sz>
__device__ __forceinline__
value_t warp_reduce_min(cg::thread_block_tile<tile_sz> g, value_t val) {
    value_t my_min = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
		const value_t other_min = g.shfl_down(my_min, i);
        my_min = my_min > other_min ? other_min : my_min;
    }
    return my_min;
}

template<typename arg_t, typename value_t, int tile_sz>
__device__ __forceinline__
arg_t warp_reduce_arg_min(cg::thread_block_tile<tile_sz> g, arg_t arg, value_t val) {
    value_t my_min = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
    cg::coalesced_group act = cg::coalesced_threads();
		const value_t other_min = g.shfl_down(my_min, i);
		const arg_t other_arg = g.shfl_down(arg, i);
		arg = my_min > other_min ? other_arg : arg; 
        my_min = my_min > other_min ? other_min: my_min; 
        act.sync();
    }
    return arg;
}

template<typename arg_t, typename value_t, int tile_sz>
__device__ __forceinline__
arg_t warp_reduce_arg_min(cg::thread_block_tile<tile_sz> g, arg_t arg, value_t val, value_t* out_val) {
    value_t my_min = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
    cg::coalesced_group act = cg::coalesced_threads();
		const value_t other_min = g.shfl_down(my_min, i);
		const arg_t other_arg = g.shfl_down(arg, i);
		arg = my_min > other_min ? other_arg : arg; 
        my_min = my_min > other_min ? other_min: my_min; 
        act.sync();
    }
    *out_val = my_min;
    return arg;
}

template<typename arg_t, typename value_t>
__device__ __forceinline__
arg_t warp_reduce_arg_min(arg_t arg, value_t val) {
    cg::coalesced_group g = cg::coalesced_threads();

    value_t my_min = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
		const value_t other_min = g.shfl_down(my_min, i);
		const arg_t other_arg = g.shfl_down(arg, i);
		arg = my_min > other_min ? other_arg : arg; 
        my_min = min(my_min, other_min);
    }
    return arg;
}


template<typename value_t>
__device__ __forceinline__
value_t warp_reduce_arg_max(value_t val) {
    cg::coalesced_group act = cg::coalesced_threads();
    value_t my_max = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = act.size() / 2; i > 0; i >>= 1) {
        cg::coalesced_group g = cg::coalesced_threads();
        const value_t other_max = g.shfl_down(my_max, i);
        my_max = max(my_max, other_max);
    }
    return my_max;
}



template<typename value_t, int tile_sz>
__device__ __forceinline__
value_t warp_reduce_max(cg::thread_block_tile<tile_sz> g, value_t val) {
    value_t my_max = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
		const value_t other_max = g.shfl_down(my_max, i);
        my_max = max(my_max, other_max);
    }
    return my_max;
}

template<typename arg_t, typename value_t, int tile_sz>
__device__ __forceinline__
arg_t warp_reduce_arg_max(cg::thread_block_tile<tile_sz> g, arg_t arg, value_t val) {
    value_t my_max = val;
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
		const value_t other_max = g.shfl_down(my_max, i);
		const arg_t other_arg = g.shfl_down(arg, i);
		arg = my_max < other_max ? other_arg : arg; 
        my_max = max(my_max, other_max);
    }
    return arg;
}


template<typename arg_t, typename value_t, int tile_sz>
__device__ __forceinline__
arg_t warp_reduce_arg_max(cg::thread_block_tile<tile_sz> g, arg_t arg, value_t val, int& out_idx) 
{
    value_t my_max = val;
	int idx = g.thread_rank();
    // Find nn node with max total value using warp level shuffling
	#pragma unroll
	for (int i = g.size() / 2; i > 0; i >>= 1) {
		const value_t other_max = g.shfl_down(my_max, i);
		const arg_t other_arg = g.shfl_down(arg, i);
		const int other_idx = g.shfl_down(idx, i);
		if(my_max < other_max){
			arg = other_arg;
			my_max = other_max;
			idx = other_idx;
		}
    }
	out_idx = idx;
    return arg;
}


__device__ 
void sort(int* data, int len, const cg::thread_block_tile<WARP_SIZE>& tile32)
{
    unsigned int counter = len;
    int8_t isOdd = 0;
    
    cg::sync(tile32);
    do{
        isOdd = (++isOdd) % 2;
        cg::sync(tile32);
        for(int idx=tile32.thread_rank(); idx < len / 2; idx+=tile32.size()){
            if(isOdd+2*idx+1>=len) {continue;}    //boundary
            int d0 = data[isOdd+2*idx];
            int d1 = data[isOdd+2*idx+1];
            
            if(d0 > d1){
                data[isOdd+2*idx] = d1;
                data[isOdd+2*idx+1] = d0;
            }
        }
    }while(counter--);
   
    cg::sync(tile32);
}
