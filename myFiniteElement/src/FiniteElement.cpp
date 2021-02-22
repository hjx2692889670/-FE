/*
 * FiniteElement.cpp
 *
 *  Created on: 2021年1月24日
 *      Author: sds
 */
#include "FiniteElement.h"
FiniteElement::FiniteElement():_area(0),SUM_B1S(nullptr),B1S(nullptr)
{
	// TODO Auto-generated constructor stub
	//一个整体
	SUM_B1S=new double*[2*N];
 	 for(int i=0;i<2*N;i++)
 	 {
 		SUM_B1S[i]=new  double[2*N] {0};
 	 }
 		B1S=new double *[T];
 		for(int j=0;j<T;j++)
 		{
 			B1S[j]=new double[T] {0};
 		}
}

FiniteElement::~FiniteElement() {
	// TODO Auto-generated destructor stub
	 for(int i=0;i<T;i++)
	 {
		 delete []B1S[i];
	 }
	 for(int i=0;i<T;i++)
	 {
		 delete []SUM_B1S[i];
	 }
}
void FiniteElement::doInit(ifstream & bmIn)
{
	   string line;
	   int t=0;
	    if(!bmIn.is_open())
	    	cerr<<"open failed"<<endl;
	    while(getline(bmIn,line)&&t<3)
	    {
	        istringstream record(line);
	        string temp;
	        int n=0;
	        while(record>>temp){
	        	BM[t][n]=atoi(temp.c_str());
	        	if(n>N)break;
	        	n++;
	        }
	        t++;
	    }
}

void FiniteElement::doB1S( ifstream &infile)
{
	int arrxy[6]={0};
	   string line;
	    if(!infile.is_open()&&!infile.is_open())
	    	cerr<<"open failed"<<endl;
	    	getline(infile,line);
	        istringstream record(line);
	        string temp;
	        int i=0;
	        while(record>>temp){
	        		arrxy[i]=atoi(temp.c_str());
	        	i++;
	        }
		doFE(arrxy[0], arrxy[1], arrxy[2], arrxy[3], arrxy[4], arrxy[5]);
}
void FiniteElement::doFE(int ix,int iy,int jx,int jy,int mx,int my,double u)
{
//先创建三个节点的x、y坐标
//结点1
	int x1=ix, y1=iy;//1,1         0,0
//结点2
	int x2=jx, y2=jy;//0,0         1,1
//结点3
	int x3=mx, y3=my;//1,0         0,1

	int b1 = y2 - y3;//0     0
	int b2 = y3 - y1;//-1   1
	int b3 = y1 - y2;//1   -1

	int c1 = x3 - x2;//1   -1
	int c2 = x1 - x3;//0   0
	int c3 = x2 - x1;//-1  1

//算出该三角形的边长length1、2、3；
	double length1 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	double length2 = sqrt((x2 - x3)*(x2 - x3) + (y2 - y3)*(y2 - y3));
	double length3 = sqrt((x1 - x3)*(x1 - x3) + (y1 - y3)*(y1 - y3));
//算出该三角形的周长per/2
	double per = (length1 + length2 + length3) / 2;
//算出该三角形的面积area；//海伦公式
	double area = sqrt(per*(per - length1)*(per - length2)*(per - length3));
	this->_area=area;
//建立几何矩阵（用二维数组实现）
	double B[3][6]={0.0};

//二维数组初始化
	double cd = 1 / (2 * area);
//根据课本中的3-9公式//就是矩阵B
//第一行
	B[0][0] = b1*cd;
	B[0][2] = b2*cd;
	B[0][4] = b3*cd;
//第二行
	B[1][1] = c1*cd;
	B[1][3] = c2*cd;
	B[1][5] = c3*cd;
//第三行
	B[2][0] = c1*cd;
	B[2][1] = b1*cd;
	B[2][2] = c2*cd;
	B[2][3] = b2*cd;
	B[2][4] = c3*cd;
	B[2][5] = b3*cd;
//完成几何矩阵的创建；

//计算应力弹性矩阵
double D[3][3] = {
		{ 1.0, u, 0.0 }, \
		{u, 1.0, 0.0}, \
		{0.0, 0.0, (1 - u) / 2
	}
};
for (int i = 0; i < (sizeof(D) / sizeof(D[0])); i++){
	D[0][i] *= (E / (1 - u*u));
	D[1][i] *= (E / (1 - u*u));
	D[2][i] *= (E / (1 - u*u));
}//求得最终的应力弹性矩阵；
//计算应力转换矩阵
	double S[3][6]={0.0};
//S=D*B;🆗
for (int j = 0; j < (sizeof(D) / sizeof(D[0])); j++){
	for (int i = 0; i < (sizeof(D) / sizeof(D[0])); i++){
		S[j][0] += D[j][i] * B[i][0];
		S[j][1] += D[j][i] * B[i][1];
		S[j][2] += D[j][i] * B[i][2];
		S[j][3] += D[j][i] * B[i][3];
		S[j][4] += D[j][i] * B[i][4];
		S[j][5] += D[j][i] * B[i][5];
	}
}
//将B、S及t，面积area带入公式得单元刚度矩阵Ke；
//求得B的转置矩阵B1；
	double B1[6][3]={0.0};

for (int i = 0; i < (sizeof(B1) / sizeof(B1[0])); i++){
	B1[i][0] = B[0][i];
	B1[i][1] = B[1][i];
	B1[i][2] = B[2][i];
}
//6:k=B逆*D*B
for (int j = 0; j < 6; j++)
{
		for (int i = 0;  i< 3; i++)
		{
			B1S[j][0] += B1[j][i] * S[i][0];
			B1S[j][1] += B1[j][i] * S[i][1];
			B1S[j][2] += B1[j][i] * S[i][2];
			B1S[j][3] += B1[j][i] * S[i][3];
			B1S[j][4] += B1[j][i] * S[i][4];
			B1S[j][5] += B1[j][i] * S[i][5];
		}
}
for (int j = 0; j < T; j++)
	{
		for (int i = 0; i < T; i++)
		{
			B1S[j][i]*= _area*t;
		}
	}
	//cout<<n<<endl;//检查输出是否正确
}
//    组装劲度矩阵
void FiniteElement::doSumFE()
{
	int i,j,m;
			for(int k=0;k<M;k++)
			{
				i=2*BM[0][k]-2;
				j=2*BM[1][k]-2;
				m=2*BM[2][k]-2;

				SUM_B1S[i][i]+=B1S[0][0];
				SUM_B1S[i][i+1]+=B1S[0][1];
				SUM_B1S[i+1][i]+=B1S[1][0];
				SUM_B1S[i+1][i+1]+=B1S[1][1];//k ii

				SUM_B1S[i][j]+=B1S[0][2];
				SUM_B1S[i][j+1]+=B1S[0][3];
				SUM_B1S[i+1][j]+=B1S[1][2];
				SUM_B1S[i+1][j+1]+=B1S[1][2];//k ij

				SUM_B1S[i][m]+=B1S[0][4];
				SUM_B1S[i][m+1]+=B1S[0][5];
				SUM_B1S[i+1][m]+=B1S[1][4];
				SUM_B1S[i+1][m+1]+=B1S[1][5];//k im

				SUM_B1S[j][i]+=B1S[2][0];
				SUM_B1S[j][i+1]+=B1S[2][1];
				SUM_B1S[j+1][i]+=B1S[3][0];
				SUM_B1S[j+1][i+1]+=B1S[3][1];//k ji

				SUM_B1S[j][j]+=B1S[2][2];
				SUM_B1S[j][j+1]+=B1S[2][3];
				SUM_B1S[j+1][j]+=B1S[3][2];
				SUM_B1S[j+1][j+1]+=B1S[3][3];//k jj

				SUM_B1S[j][m]+=B1S[2][4];
				SUM_B1S[j][m+1]+=B1S[2][5];
				SUM_B1S[j+1][m]+=B1S[3][4];
				SUM_B1S[j+1][m+1]+=B1S[3][5];//k jm

				SUM_B1S[m][i]+=B1S[4][0];
				SUM_B1S[m][i+1]+=B1S[4][1];
				SUM_B1S[m+1][i]+=B1S[5][0];
				SUM_B1S[m+1][i+1]+=B1S[5][1];//k mi

				SUM_B1S[m][j]+=B1S[4][2];
				SUM_B1S[m][j+1]+=B1S[4][3];
				SUM_B1S[m+1][j]+=B1S[5][2];
				SUM_B1S[m+1][j+1]+=B1S[5][3];//k mj

				SUM_B1S[m][m]+=B1S[4][4];
				SUM_B1S[m][m+1]+=B1S[4][5];
				SUM_B1S[m+1][m]+=B1S[5][4];
				SUM_B1S[m+1][m+1]+=B1S[5][5];//k mm
			}
}

void FiniteElement::showSUMB1S(ofstream &os)
{
	cout<<"ok Ke"<<endl;
	os<<"------------------------------------------"<<endl;
		for (int i = 0; i < 2*N; i++)
		{
			for (int j = 0; j <2*N; j++)
			{
				os<<setw(8)<<SUM_B1S[i][j]<<" ";//为最终结果！
			}
		os<<endl;
	}
}
void FiniteElement::showB1S(ofstream &os)
{
	cout<<"ok unit Ke "<<endl;
	os<<"------------------------------------------"<<endl;
		for (int i = 0; i < T; i++)
		{
			for (int j = 0; j < T; j++)
			{
				os<<setw(8)<<B1S[i][j]<<" ";//为最终结果！
			}
		os<<endl;
	}
}
