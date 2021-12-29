#include<iostream>
using namespace std;
#include"arc.h"
#include"Line.h"
#include <cmath>
#define _USE_MATH_DEFINES

Line::Line()//默认构造函数
{

}

//构造函数，用初始化直线的两个端点，pt1赋值给起点，pt2赋值给终点
Line::Line(Point pt1, Point pt2)
{
	m_ptf = pt1;
	m_pts = pt2;
}

//输出直线的一些信息，做测试用
void Line::PrintLine()
{
	cout << "第一个点:";
	m_ptf.PrintPoint();
	cout << endl;
	cout << "第二个点:";
	m_pts.PrintPoint();
	cout << endl;
	cout << "直线的长度是:" << Length() << endl;
}

//计算并返回直线的长度，返回值为直线的长度
double Line::Length()
{
	double a = (m_pts.m_x - m_ptf.m_x);
	double b = (m_pts.m_y - m_ptf.m_y);
	return sqrt(a*a + b * b);
}

//计算直线某比例处的点坐标，比例值是从起点开始计算
//参数a表示比例，在（0—1）范围内
Point Line::proportionPt(double a)
{
	if (a == 0)
	{
		return this->m_ptf;
	}
	else if(a == 1)
	{
		return this->m_pts;
	}
	else
	{
		double len = this->Length();//直线长度
		double x1 = len * a;
		double k = x1 / len;//比率

		Point fs;//表示直线起点和端点构成的向量
		fs.m_x = (m_pts.m_x - m_ptf.m_x);
		fs.m_y = (m_pts.m_y - m_ptf.m_y);

		Point fp;//表示直线起点和比例处点构成的向量
		fp.m_x = (k * fs.m_x);
		fp.m_y = (k * fs.m_y);

		Point of = m_ptf;//表示原点和起点构成的向量
		Point op;//表示原点和比例处p点构成的向量,即比例处点的坐标
		op.m_x = (of.m_x + fp.m_x);
		op.m_y = (of.m_y + fp.m_y);
		return op;
	}
	
}

/*判断点是否在线段L上,参数pt表示要确定的点
返回值为index，当值为-1时，表示点不在线上，当值为1时，表示点在线上*/
int Line::IsInLine(Point pt)const
{
	int index = 0;
	if ((pt.m_x - this->m_ptf.m_x) * (this->m_pts.m_y - this->m_ptf.m_y) == (this->m_pts.m_x - this->m_ptf.m_x) * (pt.m_y - this->m_ptf.m_y)  //叉乘 
		//保证Q点坐标在ps,pf之间 
		&& min(this->m_ptf.m_x, this->m_pts.m_x) <= pt.m_x && pt.m_x <= max(this->m_ptf.m_x, this->m_pts.m_x)
		&& min(this->m_ptf.m_y, this->m_pts.m_y) <= pt.m_y && pt.m_y <= max(this->m_ptf.m_y, this->m_pts.m_y))
	{
		index = 1;
		//cout << "点在线段上" << endl;
		return index;
	}
	else
	{
		index = -1;
		//cout << "点不在线段上" << endl;
		return index;
	}
}

/*线段自我平移，参数传入一个点Dis（x,y)，表示平移的方向和距离。
数值为正时表示向上（向右）平移，为负时相反*/
void Line::Translation(Point Dis)
{
	this->m_ptf.m_x += Dis.m_x;
	this->m_ptf.m_y += Dis.m_y;
	this->m_pts.m_x += Dis.m_x;
	this->m_pts.m_y += Dis.m_y;
}

/*线段平移后生成一条新线，参数传入一个点Dis（x,y)，表示平移的方向和距离。
数值为正时表示向上（向右）平移，为负时相反
返回值是新生成的线*/
Line* Line::NewLine(Point Dis)
{
	if (Dis.m_x, Dis.m_y >= -EPSINON && Dis.m_x, Dis.m_y <= EPSINON)
	{
		return NULL;
	}
	Line* L = new Line;
	L->m_ptf.m_x = this->m_ptf.m_x + Dis.m_x;
	L->m_ptf.m_y = this->m_ptf.m_x + Dis.m_x;
	L->m_pts.m_x = this->m_pts.m_x + Dis.m_x;
	L->m_pts.m_y = this->m_pts.m_x + Dis.m_x;
	return L;
}

Line* Line::CopyFun()//克隆直线
{
	if (this->m_ptf.m_x - this->m_pts.m_x >= -EPSINON && this->m_ptf.m_x - this->m_pts.m_x <= EPSINON
		&& this->m_ptf.m_y - this->m_pts.m_y >= -EPSINON && this->m_ptf.m_y - this->m_pts.m_y <= EPSINON)
	{
		return NULL;
	}
	Line* L = new Line;
	L->m_ptf = this->m_ptf;
	L->m_pts = this->m_pts;
	return L;
}

void Line::ReverseLine()//线段自我反向
{
	Point temp(0, 0);
	temp = this->m_ptf;
	this->m_ptf = this->m_pts;
	this->m_pts = temp;
}



//线段自我旋转degree(degree为弧度值）当弧度为负时表示顺时针，为正时表示逆时针；
//参数pt表示线段绕点pt旋转
void Line::RotationLine(double degree, Point pt)
{
	//逆时针
	/*this->m_ptf.m_x = (this->m_ptf.m_x - pt.m_x)*cos(degree) - (this->m_ptf.m_y - pt.m_y)*sin(degree) + pt.m_x;
	this->m_ptf.m_y = (this->m_ptf.m_y - pt.m_y)*cos(degree) + (this->m_ptf.m_x - pt.m_x)*sin(degree) + pt.m_y;

	this->m_pts.m_x = (this->m_pts.m_x - pt.m_x)*cos(degree) - (this->m_pts.m_y - pt.m_y)*sin(degree) + pt.m_x;
	this->m_pts.m_y = (this->m_pts.m_y - pt.m_y)*cos(degree) + (this->m_pts.m_x - pt.m_x)*sin(degree) + pt.m_y;*/
	
	
	double x0 = pt.m_x;
	double y0 = pt.m_y;
	double x1 = this->m_ptf.m_x;
	double x2 = this->m_pts.m_x;
	double y1 = this->m_ptf.m_y;
	double y2 = this->m_pts.m_y;

	if (degree < 0)//顺时针
	{
		this->m_ptf.m_x = (x1 - x0) * cos(degree) - (y1 - y0) * sin(degree) + x0;
		if (this->m_ptf.m_x >= -EPSINON && this->m_ptf.m_x <= EPSINON)
		{
			this->m_ptf.m_x = 0;
		}
		this->m_ptf.m_y = (x1 - x0) * sin(degree) + (y1 - y0) * cos(degree) + y0;
		if (this->m_ptf.m_y >= -EPSINON && this->m_ptf.m_y <= EPSINON)
		{
			this->m_ptf.m_y = 0;
		}

		this->m_pts.m_x = (x2 - x0) * cos(degree) - (y2 - y0) * sin(-degree) + x0;
		if (this->m_pts.m_x >= -EPSINON && this->m_pts.m_x <= EPSINON)
		{
			this->m_pts.m_x = 0;
		}
		this->m_pts.m_y = (x2 - x0) * sin(-degree) + (y2 - y0) * cos(degree) + y0;
		if (this->m_pts.m_y >= -EPSINON && this->m_pts.m_y <= EPSINON)
		{
			this->m_pts.m_y = 0;
		}
	}
	else if (degree > 0)//逆时针旋转
	{
		this->m_ptf.m_x = (x1 - x0)*cos(degree) - (y1 - y0)*sin(degree) + x0;
		if (m_ptf.m_x >= -EPSINON && this->m_ptf.m_x <= EPSINON)
		{
			this->m_ptf.m_x = 0;
		}
		this->m_ptf.m_y = (x1 - x0)*sin(degree) + (y1 - y0)*cos(degree) + y0;
		if (this->m_ptf.m_y >= -EPSINON && this->m_ptf.m_y <= EPSINON)
		{
			this->m_ptf.m_y = 0;
		}

		this->m_pts.m_x = (x2 - x0)*cos(degree) - (y2 - y0)*sin(degree) + pt.m_x;
		if (this->m_pts.m_x >= -EPSINON && this->m_pts.m_x <= EPSINON)
		{
			this->m_pts.m_x = 0;
		}
		this->m_pts.m_y = (x2 - x0)*sin(degree) + (y2 - y0)*cos(degree) + pt.m_y;
		if (this->m_pts.m_y >= -EPSINON && this->m_pts.m_y <= EPSINON)
		{
			this->m_pts.m_y = 0;
		}
	}
	//}
	//else if (choice == -1)//逆时针旋转
	//{
		/*this->m_ptf.m_x = (x1 - x0)*cos(degree) - (y1 - y0)*sin(degree) + x0;
		this->m_ptf.m_y = (x1 - x0)*sin(degree) + (y1 - y0)*cos(degree) + y0;

		this->m_pts.m_x = (x2 - x0)*cos(degree) - (y2 - y0)*sin(degree) + pt.m_x;
		this->m_pts.m_y = (x2 - x0)*sin(degree) + (y2 - y0)*cos(degree) + pt.m_y;*/
	//}
	//else
	//{
	//	cout << "输入有误" << endl;
	//}
}

//线段自我延伸,参数a为延长长度,参数p判定延长点为那个，true表示起点，false表示终点
void Line::Extent(bool p,double a)
{	
	double len = 0;
	len = this->Length();
	
	if (p == true)//待延长的点为起点
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//起点在终点的左边，延长方向为左
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//减直线
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //起点在终点的右边，延长方向为右
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//减直线
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
	else if (p == false)//待延长的点为终点
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//起点在终点的左边，延长方向为右
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//减直线
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //起点在终点的右边，延长方向为左
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//减直线
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
}

//线段延伸生成新线;
	//参数p用来判定延长点为哪个,true表示直线的起点，false表示直线的终点
	//参数a表示延长距离
	//返回值为新线
Line Line::NewExtent(bool p, double a)
{
	Line NewLine;
	double len = 0;
	len = this->Length();

	if (p == true)//待延长的点为起点
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//起点在终点的左边，延长方向为左
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//减直线
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //起点在终点的右边，延长方向为右
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//减直线
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
	else if (p == false)//待延长的点为终点
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//起点在终点的左边，延长方向为右
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//减直线
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //起点在终点的右边，延长方向为左
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//增（或水平）直线
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//减直线
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
	return NewLine;
}

/*对线段进行镜像，参数l是一条直线，表示线段相对直线l镜像
返回值是镜像后的直线*/
Line* Line::NewMirrorLine(Line l)
{
	if ((l.m_ptf.m_x - l.m_pts.m_x),(l.m_ptf.m_y - l.m_pts.m_y) >= -EPSINON && 
		(l.m_ptf.m_x - l.m_pts.m_x), (l.m_ptf.m_y - l.m_pts.m_y) <= EPSINON)
	{
		return NULL;
	}

	/*int ret1 = this->IsInLine(l.m_ptf);
	int ret2 = this->IsInLine(l.m_pts);
	if (ret1 == 1 || ret2 == 1)
	{
		if (ret1 == 1)
		{
			if()
		}
	}*/

	double cross1 = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_ptf.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_ptf.m_y - l.m_ptf.m_y); //|AB*AP|：矢量乘
	double d2 = (l.m_pts.m_x - l.m_ptf.m_x) * (l.m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (l.m_pts.m_y - l.m_ptf.m_y); //|AB|^2：矢量AB的大小的平方
	Point c1;
	double r1 = cross1 / d2;  //相似三角形原理求出c1点的坐标
	c1.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r1;
	c1.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r1;
	double dis1 = sqrt((this->m_ptf.m_x - c1.m_x) * (this->m_ptf.m_x - c1.m_x) + (c1.m_y - this->m_ptf.m_y) * (c1.m_y - this->m_ptf.m_y));//|CP|的大小
	Line L1(this->m_ptf, c1);
	L1.Extent(false, dis1);


	double cross2 = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_pts.m_y - l.m_ptf.m_y); //|AB*AP|：矢量乘
	Point c2;
	double r2 = cross2 / d2;  //相似三角形原理求出c2点的坐标
	c2.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r2;
	c2.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r2;
	double dis2 = sqrt((this->m_pts.m_x - c2.m_x) * (this->m_pts.m_x - c2.m_x) + (c2.m_y - this->m_pts.m_y) * (c2.m_y - this->m_pts.m_y));//|CP|的大小
	Line L2(this->m_pts, c2);
	L2.Extent(false, dis2);

	Line* newLine = new Line(L1.m_pts, L2.m_pts);
	return newLine;
	
	//Line *L = new Line;
	//double A = l.m_ptf.m_x - l.m_pts.m_x;
	//double B = l.m_ptf.m_y - l.m_pts.m_y;
	//Point F (-B, A);//传入直线l的法向量
	///*设直线为ax + by - c = 0,则由已知条件可得
 //   -Bx + Ay - C = 0;且由已知条件可求得C值*/
	//double C = A * l.m_ptf.m_x + B * l.m_ptf.m_y;
	//
	////表示直线起点到传入直线l的距离
	//double Dis1 = ((-B) * this->m_ptf.m_x + A * this->m_ptf.m_y + (-C)) / sqrt((-B)*(-B) + A * A);
	////表示直线终点到传入直线l的距离
	//double Dis2 = ((-B) * this->m_pts.m_x + A * this->m_pts.m_y + (-C)) / sqrt((-B)*(-B) + A * A);

	///*L->m_ptf.m_x = ((A*A - (-B)*(-B))*this->m_ptf.m_x - 2 * (-B)*A*this->m_ptf.m_y - 2 * (-B)*(-C)) / ((-B)*(-B) + A * A);
	//L->m_ptf.m_y = (((-B)*(-B) - A * A)*this->m_ptf.m_y - 2 * (-B)*A*this->m_ptf.m_x - 2 * A*(-C)) / (((-B)*(-B) + A * A));

	//L->m_pts.m_x = ((A*A - (-B)*(-B))*this->m_pts.m_x - 2 * (-B)*A*this->m_pts.m_y - 2 * (-B)*(-C)) / ((-B)*(-B) + A * A);
	//L->m_pts.m_y = (((-B)*(-B) - A * A)*this->m_pts.m_y - 2 * (-B)*A*this->m_pts.m_x - 2 * A*(-C)) / (((-B)*(-B) + A * A));*/


	//L->m_ptf.m_x = this->m_ptf.m_x - (((-B) / sqrt((-B)*(-B) + A * A)) * 2 * Dis1);
	//L->m_ptf.m_y = this->m_ptf.m_y - ((A / sqrt((-B)*(-B) + A * A)) * 2 * Dis1);

	//L->m_pts.m_x = this->m_pts.m_x - (((-B) / sqrt((-B)*(-B) + A * A)) * 2 * Dis2);
	//L->m_pts.m_y = this->m_pts.m_y - ((A / sqrt((-B)*(-B) + A * A)) * 2 * Dis2);
	//return L;
}

bool Line::IsOverlap(Line l)//判断两条线段是否相交
{
	double a1 = this->m_pts.m_y - this->m_ptf.m_y;
	double b1 = this->m_ptf.m_x - this->m_pts.m_x;
	double c1 = this->m_ptf.m_x*this->m_pts.m_y - this->m_pts.m_x*this->m_ptf.m_y;
	double a2 = l.m_pts.m_y - l.m_ptf.m_y;
	double b2 = l.m_ptf.m_x - l.m_pts.m_x;
	double c2 = l.m_ptf.m_x*l.m_pts.m_y - l.m_pts.m_x*l.m_ptf.m_y;
	
	//两线段的向量积，为零表示平行，不为零表示不平行
	double det = a1 * b2 - a2 * b1;
	//cout << det << endl;
	if (det >= -EPSINON && det <= EPSINON)//线段平行
	{
		int ret1 = this->IsInLine(l.m_ptf);
		int ret2 = this->IsInLine(l.m_pts);
		if (ret1 == 1 || ret2 == 1) //线段重合
		{
			//cout << "重合" << endl;
			return true;
		}
		else //没有重合
		{
			return false;
		}
	}
	else //对应直线相交
	{
		Point pt;
		pt.m_x = (c1*b2 - c2 * b1) / det;
		pt.m_y = (a1*c2 - a2 * c1) / det;
		pt.PrintPoint();

		int index = this->IsInLine(pt);
		int ret = l.IsInLine(pt);
		//pt.PrintPoint();
		//判断交点在不在两条线上，同时在表示线段有交点，任意不在表示无交点
		//this->PrintLine();
		//int ret1 = this->IsInLine(pt);
		//cout << ret1 << endl;
		//l.PrintLine();
		//int ab = l.IsInLine(pt);
		
		//pt.PrintPoint();
		//cout << ab << endl;
		if (ret == 1 && index == 1)//有交点
		{
			return true;
		}
		else//无交点
		{
			return false;
		}
	}
		

	//double a1 = this->m_pts.m_y - this->m_ptf.m_y;
	//double b1 = this->m_ptf.m_x - this->m_pts.m_x;
	//double c1 = this->m_ptf.m_x*this->m_pts.m_y - this->m_pts.m_x*this->m_ptf.m_y;
	//double a2 = l.m_pts.m_y - l.m_ptf.m_y;
	//double b2 = l.m_ptf.m_x - l.m_pts.m_x;
	//double c2 = l.m_ptf.m_x*l.m_pts.m_y - l.m_pts.m_x*l.m_ptf.m_y;
	//double det = a1 * b2 - a2 * b1;

	//if (det == 0)//不相交
	//{
	//	return false;
	//}
	//else//对应直线相交，且得到对应直线的相交点
	//{
	//	pt.m_x = (c1*b2 - c2 * b1) / det;
	//	pt.m_y = (a1*c2 - a2 * c1) / det;
	//}

	////判断交点在不在两条线上，同时在表示线段有交点，任意不在表示无交点
	//int ret1 = this->IsInLine(pt);
	//int ret2 = l.IsInLine(pt);
	//if (ret1 == -1 || ret2 == -1)//无交点
	//{
	//	return false;
	//}
	//else//有交点
	//{
	//	return true;
	//}
	
	/*Point pt;
	double x1 = this->m_ptf.m_x;
	double x2 = this->m_pts.m_x;
	double x3 = l.m_ptf.m_x;
	double x4 = l.m_pts.m_x;
	double y1 = this->m_ptf.m_y;
	double y2 = this->m_pts.m_y;
	double y3 = l.m_ptf.m_y;
	double y4 = l.m_pts.m_y;
	pt.m_x = ((x2 - x1) * (x3 - x4) * (y3 - y1) -
		x3 * (x2 - x1) * (y3 - y4) + x1 * (y2 - y1) * (x3 - x4)) /
		((y2 - y1) * (x3 - x4) - (x2 - x1) * (y3 - y4));

	pt.m_y = ((y2 - y1) * (y3 - y4) * (x3 - x1) -
		y3 * (y2 - y1) * (x3 - x4) + y1 * (x2 - x1) * (y3 - y4)) /
		((y2 - y1) * (y3 - y4) - (y2 - y1) * (x3 - x4));
	int ret1 = this->IsInLine(pt);
	int ret2 = l.IsInLine(pt);
	if (ret1 == -1 || ret2 == -1)
	{
		return false ;
	}
	else
	{
		return true;
	}*/
}

//线段求交点，并返回交点坐标。
//参数是另一条直线，和点pt，返回值pt为交点坐标。
//如果存在交点，则返回交点的坐标；如果不存在交点，则返回false
bool Line::LinePoint(Line l, Point& pt)
{
	bool ret = this->IsOverlap(l);
	if (ret == false)
	{
		
		cout << "两条线段没有交点" << endl;
		return false;
	}
	else
	{
		cout << "两直线存在交点" << endl;

		double a1 = this->m_pts.m_y - this->m_ptf.m_y;
		double b1 = this->m_ptf.m_x - this->m_pts.m_x;
		double c1 = this->m_ptf.m_x*this->m_pts.m_y - this->m_pts.m_x*this->m_ptf.m_y;
		double a2 = l.m_pts.m_y - l.m_ptf.m_y;
		double b2 = l.m_ptf.m_x - l.m_pts.m_x;
		double c2 = l.m_ptf.m_x*l.m_pts.m_y - l.m_pts.m_x*l.m_ptf.m_y;
		double det = a1 * b2 - a2 * b1;
		/*double x1 = this->m_ptf.m_x;
		double x2 = this->m_pts.m_x;
		double x3 = l.m_ptf.m_x;
		double x4 = l.m_pts.m_x;
		double y1 = this->m_ptf.m_y;
		double y2 = this->m_pts.m_y;
		double y3 = l.m_ptf.m_y;
		double y4 = l.m_pts.m_y;
		pt.m_x = ((x2 - x1) * (x3 - x4) * (y3 - y1) -
			x3 * (x2 - x1) * (y3 - y4) + x1 * (y2 - y1) * (x3 - x4)) /
			((y2 - y1) * (x3 - x4) - (x2 - x1) * (y3 - y4));

		pt.m_y = ((y2 - y1) * (y3 - y4) * (x3 - x1) -
			y3 * (y2 - y1) * (x3 - x4) + y1 * (x2 - x1) * (y3 - y4)) /
			((y2 - y1) * (y3 - y4) - (y2 - y1) * (x3 - x4));*/
		pt.m_x = (c1*b2 - c2 * b1) / det;
		pt.m_y = (a1*c2 - a2 * c1) / det;
		return true;
	}
}
//Point Line::LinePoint(Line l)//线段求交点，并返回交点
//{
//	Point p(0.00, 0.00);
//	bool ret = this->IsOverlap(l);
//	if (ret == false)
//	{
//		cout << "两条线段没有交点" << endl;
//		return p;
//	}
//	else
//	{
//		cout << "两直线存在交点" << endl;
//		return pt;
//	}
//}