/*
 * FiniteElement.h
 *
 *  Created on: 2021年1月24日
 *      Author: sds
 */

#ifndef FINITEELEMENT_H_
#define FINITEELEMENT_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

#include <math.h>
#include <stdlib.h>
#include<iomanip>

using namespace std;
#define N 15   //单元
#define M 16 //三角形，
#define T 6
//static int E = 1;//假设厚度t=1，弹性模量E=1
//static int t = 1;
//三角形单元的刚度矩阵k,设厚度t=1,弹性模量E=1；
static double E = 2e+11;//假设厚度t=1，弹性模量E=1
static int t = 1;
class FiniteElement {
public:
	FiniteElement();
	void doFE(int ix,int iy,int jx,int jy,int mx,int my,double u=0.25);
	void doB1S( ifstream &infile);
	void showSUMB1S(ofstream &);
	void doSumFE();
	void showB1S(ofstream &);
	void doInit(ifstream & bmIn);
	virtual ~FiniteElement();
//辅助计算
	double _area;
	double **B1S;//只在分配内存用到过
//参与计算
	double **SUM_B1S;
private:
	int BM[3][N]={0};
};

#endif /* FINITEELEMENT_H_ */
