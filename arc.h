#pragma once
#include<iostream>
#include"Point.h"
using namespace std;
#include"Line.h"
#include<cmath>
//#include"double_equals.h"
#include<vector>
#include<algorithm>
//����һ��Բ���࣬ʹ�� Բ����Բ�ģ�Բ����ʼ����x��������ļн�(��ʱ��Ϊ������
//Բ����㵽�յ�ɨ���ĽǶ�(��ʱ��Ϊ������Բ���뾶R�ķ�ʽ��Բ����ʵ��һЩԲ���Ļ�������
class Arc
{
public:
	Arc();

	//���캯������ʼ��Բ���Ĳ�����a��ֵ�������x��������нǣ�
	//b��ֵ����㵽�յ�ɨ���ĽǶȣ�r��ֵ��Բ�İ뾶
	Arc(Point pt1, double a, double b, double r);
	
	//����Բ����㲢�����������
	Point Gain_Star();

	//����Բ���յ㲢�����յ�����
	Point Gain_Fin();

	/*����Բ��ĳ�������ĵ����꣬����ֵ�Ǵ���㿪ʼ����
	����a��ʾ��������[0��1]��Χ��
	����ֵΪ�˱������������*/
	Point proportionPt(double a);

	//���㲢����Բ���ĳ���
	double Length(); 

	//���Բ���������˵��Բ�ļ�Բ�����ȳ��Ȳ�����
	void cPrintArc(); 

	//��¡Բ��,����ֵ����Բ����ָ��
	Arc* CopyFun();

	//�жϵ�pt�Ƿ���Բ����,������ʾ���жϵĵ㣻����ֵ����1��-1.
    //�㷵��ֵΪ - 1ʱ��ʾ�㲻��Բ���ϣ�����ֵΪ1ʱ����ʾ����Բ����
	int IsInArc(Point pt);

	/*Բ������ƽ�ƣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
	��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴*/
	void Translation(Point Dis);

	/*Բ��ƽ��������Բ������������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
    ��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴������ֵΪ��Բ����ָ��*/
	Arc* NewArc(Point Dis);


	/*��Բ�����о��񣬲���l��һ��ֱ�ߣ���ʾԲ�����ֱ��l����
    ����ֵ�Ǿ�����Բ��ָ��*/
	Arc *NewMirrorArc(Line l);

	//Բ�����ҷ���
	void ReverseArc();

	//Բ��������תdegree(degreeΪ����ֵ��������Ϊ��ʱ��ʾ˳ʱ�룬Ϊ��ʱ��ʾ��ʱ�룻
    //����pt��ʾ�߶��Ƶ�pt��ת
	void RotationArc(double degree, Point pt);

	//��������Բ���Ľ��㣬����Ϊ��һ��Բ����
	//����ֵΪ�������꣬��������ڽ��㣬���ؿ�ֵ��
	//�������һ�����㣬����һ���������ꣻ
	//��������������㣬�����������������ֵ
	vector<Point> ArcPoint(Arc c1);

	//��Բ��ĳ�˵�����һ�ξ���
    //����pΪ����ֵ�����pΪtrue���ӳ���Ϊ��㣻��֮Ϊ�յ㣻
    //����aΪ�ӳ�����
    //����ֵΪ������Բ��ָ��
	void Extent(bool p, double a);

	// �ж�������������ͬ
	int double_equal_c(double a, double b);

	~Arc();

public:
	Point m_center;//Բ����Բ��
	double m_Stangle;//Բ����ʼ����x��������ļн�(��ʱ��Ϊ����
	double m_Swangle;//Բ����㵽�յ�ɨ���ĽǶ�(��ʱ��Ϊ����
	double m_R;//Բ���뾶R
};