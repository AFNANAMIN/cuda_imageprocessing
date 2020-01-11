#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
using namespace std;

#define  i_size 6//image size
#define o_size 6
#define k_size 3//kernel size
int input[i_size][i_size];
int kernel[k_size][k_size];
int output[i_size][i_size];
typedef int itype[i_size];
typedef int ktype[k_size];
void fill_image(int m[i_size][i_size]) {
	static int n = 0; int i, j;
	for (i = 0; i < i_size; i++)
		for (j = 0; j < i_size; j++)
			m[i][j] = n++;
}

void fill_kernel(int m[k_size][k_size]) {
	static int n = 0; int i, j;
	for (i = 0; i < k_size; i++)
		for (j = 0; j < k_size; j++)
			m[i][j] = n++;
}
void fill_output(int m[i_size][i_size]) {
	int i, j;

	for (i = 0; i < o_size; i++) {
		cout << "\n \t\t |";
		for (j = 0; j < o_size; j++)
			cout << "\t\t" << m[i][j];
		cout << "|";
	}

}
__global__ void add_arrays_gpu(int* a, int *b, int* c)
{
	c[threadIdx.x] = a[threadIdx.x] + b[threadIdx.x];
}
__global__ void processing(itype  *a,ktype *kernel, itype *o)
{
	
	int r = 0;
	int i = (blockIdx.y*blockDim.y + threadIdx.y)+1;
	int j = (blockIdx.x*blockDim.x + threadIdx.x)+1;
	
	
	for (int k = -1; k < 2; k++)
	{
		for (int m = -1; m < 2; m++)
		{

			r += a[i + k][j + m] * kernel[k + 1][m + 1];
			o[i][j] = r;
		}

	}
}




int main()
{fill_image(input);
fill_kernel(kernel);
cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);


	itype *device_a, *device_c;
	ktype *device_b;
	const int size = i_size*i_size * sizeof(int);
	
	size_t pitch;
	cudaMallocPitch((void**)&device_a, &pitch, i_size * sizeof(float), i_size);
	cudaMallocPitch((void**)&device_b, &pitch, k_size * sizeof(float), k_size);
	cudaMallocPitch((void**)&device_c, &pitch, i_size * sizeof(float), i_size);
	dim3 blockspergrid(2,2,1);
	dim3 threadperblock(2,2,1);
	cudaMemcpy(
		device_a, input,
		size,
		cudaMemcpyHostToDevice
	);

	cudaMemcpy(
		device_b, kernel,
		size,
		cudaMemcpyHostToDevice
	);
	int r = 0;
	//add_arra<< <1, count >> > (device_a, device_b, device_c);
	cudaEventRecord(start);
	//processing <<<blockspergrid, threadperblock >> > (device_a,device_b,device_c);
	for (int i = 1; i < i_size - 1; i++)
	{
		for (int j = 1; j < i_size - 1; j++)
			for (int k = -1; k < 2; k++)
			{
				for (int m = -1; m < 2; m++)
				{

					r += input[i + k][j + m] * kernel[k + 1][m + 1];
					output[i][j] = r;
				}

			}
	}




	cudaEventRecord(stop);
	
	/*cudaMemcpy(
		output, device_c,
		size,
		cudaMemcpyDeviceToHost
	);*/
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	cout << "\n total time to execute operation is " << milliseconds << "\n";
	fill_output(output);
	

	//getchar();

	return 0;
}