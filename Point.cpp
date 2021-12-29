#include<iostream>
using namespace std;
#include"Point.h"

Point::Point()//默认构造函数
{
	m_x = 0;
	m_y = 0;
}

Point::Point(double x, double y)//拷贝构造函数，初始化点的x、y坐标
{
	this->m_x = x;
	this->m_y = y;
}

void Point::operator=(Point &pt)//赋值操作符重载
{
	this->m_x = pt.m_x;
	this->m_y = pt.m_y;
}

bool Point::operator==(Point &pt)//==重载
{
	if(this->m_x == pt.m_x &&this->m_y == pt.m_y)
		{
			return true;
		}
		return false;
}

//“-”运算符重载
//void  Point::operator-(Point p)
//{
//
//}

void Point::PrintPoint()////输出点
{
	cout << "(" << m_x << "," << m_y << ")" << endl;
}

//计算点绕某点旋转后的坐标
//参数pt表示被绕的点；degree表示旋转的弧度，逆时针为正，反正为负
Point Point::RotationPoint(double degree, Point pt)
{
	double x0 = pt.m_x;
	double y0 = pt.m_y;
	Point newPoint;
	double x1 = this->m_x;
	double y1 = this->m_y;

	if (degree < 0)//顺时针
	{
		newPoint.m_x = (x1 - x0) * cos(degree) - (y1 - y0) * sin(degree) + x0;
		if (newPoint.m_x >= -EPSINON && this->m_x <= EPSINON)
		{
			newPoint.m_x = 0;
		}
		newPoint.m_y = (x1 - x0) * sin(degree) + (y1 - y0) * cos(degree) + y0;
		if (newPoint.m_y >= -EPSINON && newPoint.m_y <= EPSINON)
		{
			newPoint.m_y = 0;
		}
	}
	else if (degree > 0)//逆时针旋转
	{
		newPoint.m_x = (x1 - x0)*cos(degree) - (y1 - y0)*sin(degree) + x0;
		if (newPoint.m_x >= -EPSINON && newPoint.m_x <= EPSINON)
		{
			newPoint.m_x = 0;
		}
		newPoint.m_y = (x1 - x0)*sin(degree) + (y1 - y0)*cos(degree) + y0;
		if (newPoint.m_y >= -EPSINON && newPoint.m_y <= EPSINON)
		{
			newPoint.m_y = 0;
		}

	}
	return newPoint;
}
////求点关于某条直线对称后点的坐标
////传入的参数有待镜像的点pt和一条直线l，表示点pt相对直线l镜像
////返回值是点镜像后新的点
//Point Point::NewMirrorPoint(Line l)
//{
//	if ((l.m_ptf.m_x - l.m_pts.m_x), (l.m_ptf.m_y - l.m_pts.m_y) >= -EPSINON &&
//		(l.m_ptf.m_x - l.m_pts.m_x), (l.m_ptf.m_y - l.m_pts.m_y) <= EPSINON)
//	{
//		return;
//	}
//
//	double cross = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_y - l.m_ptf.m_y); //|AB*AP|：矢量乘
//	double d = (l.m_pts.m_x - l.m_ptf.m_x) * (l.m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (l.m_pts.m_y - l.m_ptf.m_y); //|AB|^2：矢量AB的大小的平方
//	Point c;
//	double r = cross / d;  //相似三角形原理求出c1点的坐标
//	c.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r;
//	c.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r;
//	double dis = sqrt((this->m_x - c.m_x) * (this->m_x - c.m_x) + (c.m_y - this->m_y) * (c.m_y - this->m_y));//|CP|的大小
//
//	Point pt;
//	pt.m_x = this->m_x;
//	pt.m_y = this->m_y;
//
//	Line L(pt, c);
//	L.Extent(false, dis);
//
//	return L.m_pts;
//}

// 浮点数判同
int Point::double_equal_p(double a, double b)
{
	static const double ZERO = 1e-9;
	return fabs(a - b) < ZERO;
}



Point::~Point()
{

}

