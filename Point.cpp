#include<iostream>
using namespace std;
#include"Point.h"

Point::Point()//Ĭ�Ϲ��캯��
{
	m_x = 0;
	m_y = 0;
}

Point::Point(double x, double y)//�������캯������ʼ�����x��y����
{
	this->m_x = x;
	this->m_y = y;
}

void Point::operator=(Point &pt)//��ֵ����������
{
	this->m_x = pt.m_x;
	this->m_y = pt.m_y;
}

bool Point::operator==(Point &pt)//==����
{
	if(this->m_x == pt.m_x &&this->m_y == pt.m_y)
		{
			return true;
		}
		return false;
}

//��-�����������
//void  Point::operator-(Point p)
//{
//
//}

void Point::PrintPoint()////�����
{
	cout << "(" << m_x << "," << m_y << ")" << endl;
}

//�������ĳ����ת�������
//����pt��ʾ���Ƶĵ㣻degree��ʾ��ת�Ļ��ȣ���ʱ��Ϊ��������Ϊ��
Point Point::RotationPoint(double degree, Point pt)
{
	double x0 = pt.m_x;
	double y0 = pt.m_y;
	Point newPoint;
	double x1 = this->m_x;
	double y1 = this->m_y;

	if (degree < 0)//˳ʱ��
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
	else if (degree > 0)//��ʱ����ת
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
////������ĳ��ֱ�߶Գƺ�������
////����Ĳ����д�����ĵ�pt��һ��ֱ��l����ʾ��pt���ֱ��l����
////����ֵ�ǵ㾵����µĵ�
//Point Point::NewMirrorPoint(Line l)
//{
//	if ((l.m_ptf.m_x - l.m_pts.m_x), (l.m_ptf.m_y - l.m_pts.m_y) >= -EPSINON &&
//		(l.m_ptf.m_x - l.m_pts.m_x), (l.m_ptf.m_y - l.m_pts.m_y) <= EPSINON)
//	{
//		return;
//	}
//
//	double cross = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_y - l.m_ptf.m_y); //|AB*AP|��ʸ����
//	double d = (l.m_pts.m_x - l.m_ptf.m_x) * (l.m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (l.m_pts.m_y - l.m_ptf.m_y); //|AB|^2��ʸ��AB�Ĵ�С��ƽ��
//	Point c;
//	double r = cross / d;  //����������ԭ�����c1�������
//	c.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r;
//	c.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r;
//	double dis = sqrt((this->m_x - c.m_x) * (this->m_x - c.m_x) + (c.m_y - this->m_y) * (c.m_y - this->m_y));//|CP|�Ĵ�С
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

// ��������ͬ
int Point::double_equal_p(double a, double b)
{
	static const double ZERO = 1e-9;
	return fabs(a - b) < ZERO;
}



Point::~Point()
{

}

