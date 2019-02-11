/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/
//__kernel void benchmarking() { 
//	return;
//}
//
//
//void matrixMultiplication(__global float* A, __global float* B, __global float* C, int widthA, int widthB) {
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	float value = 0;
//	for (int k = 0; k < widthA; k++)
//	{
//		value = value + A[k + j * widthA] * B[k*widthB + i];
//	}
//	C[i + widthA * j] = value;
//}
// TODO it should be removed to other ocl file 
typedef struct {
	int width;
	int height;
	__global float* elements;
} Matrix;
// Thread block size
#define BLOCK_SIZE 16
// Matrix multiplication function called by MatMulKernel()
void matrixMul(Matrix A, Matrix B, Matrix C)
{
	float Cvalue = 0;
	int row = get_global_id(1);
	int col = get_global_id(0);
	for (int e = 0; e < A.width; ++e)
		Cvalue += A.elements[row * A.width + e]
		* B.elements[e * B.width + col];
	C.elements[row * C.width + col] = Cvalue;
}
// Matrix multiplication kernel called by MatMulHost()
__kernel void MatMulKernel(
	int Awidth, int Aheight, __global float* Aelements,
	int Bwidth, int Bheight, __global float* Belements,
	int Cwidth, int Cheight, __global float* Celements,
	int factor)
{
	Matrix A = { Awidth, Aheight, Aelements };
	Matrix B = { Bwidth, Bheight, Belements };
	Matrix C = { Cwidth, Cheight, Celements };
	for (int i = 0; i < factor; ++i) { 
		matrixMul(A, B, C);
	}
}