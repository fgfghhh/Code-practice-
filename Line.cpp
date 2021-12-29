#include<iostream>
using namespace std;
#include"arc.h"
#include"Line.h"
#include <cmath>
#define _USE_MATH_DEFINES

Line::Line()//Ĭ�Ϲ��캯��
{

}

//���캯�����ó�ʼ��ֱ�ߵ������˵㣬pt1��ֵ����㣬pt2��ֵ���յ�
Line::Line(Point pt1, Point pt2)
{
	m_ptf = pt1;
	m_pts = pt2;
}

//���ֱ�ߵ�һЩ��Ϣ����������
void Line::PrintLine()
{
	cout << "��һ����:";
	m_ptf.PrintPoint();
	cout << endl;
	cout << "�ڶ�����:";
	m_pts.PrintPoint();
	cout << endl;
	cout << "ֱ�ߵĳ�����:" << Length() << endl;
}

//���㲢����ֱ�ߵĳ��ȣ�����ֵΪֱ�ߵĳ���
double Line::Length()
{
	double a = (m_pts.m_x - m_ptf.m_x);
	double b = (m_pts.m_y - m_ptf.m_y);
	return sqrt(a*a + b * b);
}

//����ֱ��ĳ�������ĵ����꣬����ֵ�Ǵ���㿪ʼ����
//����a��ʾ�������ڣ�0��1����Χ��
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
		double len = this->Length();//ֱ�߳���
		double x1 = len * a;
		double k = x1 / len;//����

		Point fs;//��ʾֱ�����Ͷ˵㹹�ɵ�����
		fs.m_x = (m_pts.m_x - m_ptf.m_x);
		fs.m_y = (m_pts.m_y - m_ptf.m_y);

		Point fp;//��ʾֱ�����ͱ������㹹�ɵ�����
		fp.m_x = (k * fs.m_x);
		fp.m_y = (k * fs.m_y);

		Point of = m_ptf;//��ʾԭ�����㹹�ɵ�����
		Point op;//��ʾԭ��ͱ�����p�㹹�ɵ�����,���������������
		op.m_x = (of.m_x + fp.m_x);
		op.m_y = (of.m_y + fp.m_y);
		return op;
	}
	
}

/*�жϵ��Ƿ����߶�L��,����pt��ʾҪȷ���ĵ�
����ֵΪindex����ֵΪ-1ʱ����ʾ�㲻�����ϣ���ֵΪ1ʱ����ʾ��������*/
int Line::IsInLine(Point pt)const
{
	int index = 0;
	if ((pt.m_x - this->m_ptf.m_x) * (this->m_pts.m_y - this->m_ptf.m_y) == (this->m_pts.m_x - this->m_ptf.m_x) * (pt.m_y - this->m_ptf.m_y)  //��� 
		//��֤Q��������ps,pf֮�� 
		&& min(this->m_ptf.m_x, this->m_pts.m_x) <= pt.m_x && pt.m_x <= max(this->m_ptf.m_x, this->m_pts.m_x)
		&& min(this->m_ptf.m_y, this->m_pts.m_y) <= pt.m_y && pt.m_y <= max(this->m_ptf.m_y, this->m_pts.m_y))
	{
		index = 1;
		//cout << "�����߶���" << endl;
		return index;
	}
	else
	{
		index = -1;
		//cout << "�㲻���߶���" << endl;
		return index;
	}
}

/*�߶�����ƽ�ƣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴*/
void Line::Translation(Point Dis)
{
	this->m_ptf.m_x += Dis.m_x;
	this->m_ptf.m_y += Dis.m_y;
	this->m_pts.m_x += Dis.m_x;
	this->m_pts.m_y += Dis.m_y;
}

/*�߶�ƽ�ƺ�����һ�����ߣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴
����ֵ�������ɵ���*/
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

Line* Line::CopyFun()//��¡ֱ��
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

void Line::ReverseLine()//�߶����ҷ���
{
	Point temp(0, 0);
	temp = this->m_ptf;
	this->m_ptf = this->m_pts;
	this->m_pts = temp;
}



//�߶�������תdegree(degreeΪ����ֵ��������Ϊ��ʱ��ʾ˳ʱ�룬Ϊ��ʱ��ʾ��ʱ�룻
//����pt��ʾ�߶��Ƶ�pt��ת
void Line::RotationLine(double degree, Point pt)
{
	//��ʱ��
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

	if (degree < 0)//˳ʱ��
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
	else if (degree > 0)//��ʱ����ת
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
	//else if (choice == -1)//��ʱ����ת
	//{
		/*this->m_ptf.m_x = (x1 - x0)*cos(degree) - (y1 - y0)*sin(degree) + x0;
		this->m_ptf.m_y = (x1 - x0)*sin(degree) + (y1 - y0)*cos(degree) + y0;

		this->m_pts.m_x = (x2 - x0)*cos(degree) - (y2 - y0)*sin(degree) + pt.m_x;
		this->m_pts.m_y = (x2 - x0)*sin(degree) + (y2 - y0)*cos(degree) + pt.m_y;*/
	//}
	//else
	//{
	//	cout << "��������" << endl;
	//}
}

//�߶���������,����aΪ�ӳ�����,����p�ж��ӳ���Ϊ�Ǹ���true��ʾ��㣬false��ʾ�յ�
void Line::Extent(bool p,double a)
{	
	double len = 0;
	len = this->Length();
	
	if (p == true)//���ӳ��ĵ�Ϊ���
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//������յ����ߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//��ֱ��
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //������յ���ұߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//��ֱ��
			{
				this->m_pts.m_x = this->m_pts.m_x;
				this->m_pts.m_y = this->m_pts.m_y;
				this->m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
	else if (p == false)//���ӳ��ĵ�Ϊ�յ�
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//������յ����ߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//��ֱ��
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //������յ���ұߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//��ֱ��
			{
				this->m_ptf.m_x = this->m_ptf.m_x;
				this->m_ptf.m_y = this->m_ptf.m_y;
				this->m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				this->m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
}

//�߶�������������;
	//����p�����ж��ӳ���Ϊ�ĸ�,true��ʾֱ�ߵ���㣬false��ʾֱ�ߵ��յ�
	//����a��ʾ�ӳ�����
	//����ֵΪ����
Line Line::NewExtent(bool p, double a)
{
	Line NewLine;
	double len = 0;
	len = this->Length();

	if (p == true)//���ӳ��ĵ�Ϊ���
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//������յ����ߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//��ֱ��
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x - a * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //������յ���ұߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y + a * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//��ֱ��
			{
				NewLine.m_pts.m_x = this->m_pts.m_x;
				NewLine.m_pts.m_y = this->m_pts.m_y;
				NewLine.m_ptf.m_x = this->m_ptf.m_x + a * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_ptf.m_y = this->m_ptf.m_y - a * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
		}
	}
	else if (p == false)//���ӳ��ĵ�Ϊ�յ�
	{
		if (this->m_ptf.m_x < this->m_pts.m_x)//������յ����ߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y + (a) * ((this->m_pts.m_y - this->m_ptf.m_y) / len);
			}
			else//��ֱ��
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x + (a) * ((this->m_pts.m_x - this->m_ptf.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
		}
		else //������յ���ұߣ��ӳ�����Ϊ��
		{
			if (this->m_pts.m_y >= this->m_ptf.m_y)//������ˮƽ��ֱ��
			{
				NewLine.m_ptf.m_x = this->m_ptf.m_x;
				NewLine.m_ptf.m_y = this->m_ptf.m_y;
				NewLine.m_pts.m_x = this->m_pts.m_x - (a) * ((this->m_ptf.m_x - this->m_pts.m_x) / len);
				NewLine.m_pts.m_y = this->m_pts.m_y - (a) * ((this->m_ptf.m_y - this->m_pts.m_y) / len);
			}
			else//��ֱ��
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

/*���߶ν��о��񣬲���l��һ��ֱ�ߣ���ʾ�߶����ֱ��l����
����ֵ�Ǿ�����ֱ��*/
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

	double cross1 = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_ptf.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_ptf.m_y - l.m_ptf.m_y); //|AB*AP|��ʸ����
	double d2 = (l.m_pts.m_x - l.m_ptf.m_x) * (l.m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (l.m_pts.m_y - l.m_ptf.m_y); //|AB|^2��ʸ��AB�Ĵ�С��ƽ��
	Point c1;
	double r1 = cross1 / d2;  //����������ԭ�����c1�������
	c1.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r1;
	c1.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r1;
	double dis1 = sqrt((this->m_ptf.m_x - c1.m_x) * (this->m_ptf.m_x - c1.m_x) + (c1.m_y - this->m_ptf.m_y) * (c1.m_y - this->m_ptf.m_y));//|CP|�Ĵ�С
	Line L1(this->m_ptf, c1);
	L1.Extent(false, dis1);


	double cross2 = (l.m_pts.m_x - l.m_ptf.m_x) * (this->m_pts.m_x - l.m_ptf.m_x) + (l.m_pts.m_y - l.m_ptf.m_y) * (this->m_pts.m_y - l.m_ptf.m_y); //|AB*AP|��ʸ����
	Point c2;
	double r2 = cross2 / d2;  //����������ԭ�����c2�������
	c2.m_x = l.m_ptf.m_x + (l.m_pts.m_x - l.m_ptf.m_x) * r2;
	c2.m_y = l.m_ptf.m_y + (l.m_pts.m_y - l.m_ptf.m_y) * r2;
	double dis2 = sqrt((this->m_pts.m_x - c2.m_x) * (this->m_pts.m_x - c2.m_x) + (c2.m_y - this->m_pts.m_y) * (c2.m_y - this->m_pts.m_y));//|CP|�Ĵ�С
	Line L2(this->m_pts, c2);
	L2.Extent(false, dis2);

	Line* newLine = new Line(L1.m_pts, L2.m_pts);
	return newLine;
	
	//Line *L = new Line;
	//double A = l.m_ptf.m_x - l.m_pts.m_x;
	//double B = l.m_ptf.m_y - l.m_pts.m_y;
	//Point F (-B, A);//����ֱ��l�ķ�����
	///*��ֱ��Ϊax + by - c = 0,������֪�����ɵ�
 //   -Bx + Ay - C = 0;������֪���������Cֵ*/
	//double C = A * l.m_ptf.m_x + B * l.m_ptf.m_y;
	//
	////��ʾֱ����㵽����ֱ��l�ľ���
	//double Dis1 = ((-B) * this->m_ptf.m_x + A * this->m_ptf.m_y + (-C)) / sqrt((-B)*(-B) + A * A);
	////��ʾֱ���յ㵽����ֱ��l�ľ���
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

bool Line::IsOverlap(Line l)//�ж������߶��Ƿ��ཻ
{
	double a1 = this->m_pts.m_y - this->m_ptf.m_y;
	double b1 = this->m_ptf.m_x - this->m_pts.m_x;
	double c1 = this->m_ptf.m_x*this->m_pts.m_y - this->m_pts.m_x*this->m_ptf.m_y;
	double a2 = l.m_pts.m_y - l.m_ptf.m_y;
	double b2 = l.m_ptf.m_x - l.m_pts.m_x;
	double c2 = l.m_ptf.m_x*l.m_pts.m_y - l.m_pts.m_x*l.m_ptf.m_y;
	
	//���߶ε���������Ϊ���ʾƽ�У���Ϊ���ʾ��ƽ��
	double det = a1 * b2 - a2 * b1;
	//cout << det << endl;
	if (det >= -EPSINON && det <= EPSINON)//�߶�ƽ��
	{
		int ret1 = this->IsInLine(l.m_ptf);
		int ret2 = this->IsInLine(l.m_pts);
		if (ret1 == 1 || ret2 == 1) //�߶��غ�
		{
			//cout << "�غ�" << endl;
			return true;
		}
		else //û���غ�
		{
			return false;
		}
	}
	else //��Ӧֱ���ཻ
	{
		Point pt;
		pt.m_x = (c1*b2 - c2 * b1) / det;
		pt.m_y = (a1*c2 - a2 * c1) / det;
		pt.PrintPoint();

		int index = this->IsInLine(pt);
		int ret = l.IsInLine(pt);
		//pt.PrintPoint();
		//�жϽ����ڲ����������ϣ�ͬʱ�ڱ�ʾ�߶��н��㣬���ⲻ�ڱ�ʾ�޽���
		//this->PrintLine();
		//int ret1 = this->IsInLine(pt);
		//cout << ret1 << endl;
		//l.PrintLine();
		//int ab = l.IsInLine(pt);
		
		//pt.PrintPoint();
		//cout << ab << endl;
		if (ret == 1 && index == 1)//�н���
		{
			return true;
		}
		else//�޽���
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

	//if (det == 0)//���ཻ
	//{
	//	return false;
	//}
	//else//��Ӧֱ���ཻ���ҵõ���Ӧֱ�ߵ��ཻ��
	//{
	//	pt.m_x = (c1*b2 - c2 * b1) / det;
	//	pt.m_y = (a1*c2 - a2 * c1) / det;
	//}

	////�жϽ����ڲ����������ϣ�ͬʱ�ڱ�ʾ�߶��н��㣬���ⲻ�ڱ�ʾ�޽���
	//int ret1 = this->IsInLine(pt);
	//int ret2 = l.IsInLine(pt);
	//if (ret1 == -1 || ret2 == -1)//�޽���
	//{
	//	return false;
	//}
	//else//�н���
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

//�߶��󽻵㣬�����ؽ������ꡣ
//��������һ��ֱ�ߣ��͵�pt������ֵptΪ�������ꡣ
//������ڽ��㣬�򷵻ؽ�������ꣻ��������ڽ��㣬�򷵻�false
bool Line::LinePoint(Line l, Point& pt)
{
	bool ret = this->IsOverlap(l);
	if (ret == false)
	{
		
		cout << "�����߶�û�н���" << endl;
		return false;
	}
	else
	{
		cout << "��ֱ�ߴ��ڽ���" << endl;

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
//Point Line::LinePoint(Line l)//�߶��󽻵㣬�����ؽ���
//{
//	Point p(0.00, 0.00);
//	bool ret = this->IsOverlap(l);
//	if (ret == false)
//	{
//		cout << "�����߶�û�н���" << endl;
//		return p;
//	}
//	else
//	{
//		cout << "��ֱ�ߴ��ڽ���" << endl;
//		return pt;
//	}
//}