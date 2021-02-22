//============================================================================
// Name        : myFiniteElement.cpp
// Author      : 黄俊想 2021 1 24
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "FiniteElement.h"
//三角形单元的刚度矩阵k,设厚度t=1,弹性模量E；
int main()
{
	FiniteElement obj;
	const string infile="in.txt";
	const string bmIn="BM.txt";
	const string outDate="OutData.txt";
	ifstream in(infile);
	ifstream bm(bmIn);
	ofstream OD(outDate);
	obj.doInit(bm);
	obj.doB1S(in);
	obj.doSumFE();
	obj.showB1S(OD);
	obj.showSUMB1S(OD);
	in.close();
	bm.close();
	OD.close();
	return 0;
}
