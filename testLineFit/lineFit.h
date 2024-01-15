/*
 * @Author: niew nie_w@reachauto.com
 * @Date: 2023-04-20 15:09:32
 * @LastEditors: niew nie_w@reachauto.com
 * @LastEditTime: 2023-04-20 15:16:22
 * @FilePath: \roadFusion\includes\lineFit.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A
 */
#ifndef _LINE_FITH_
#define _LINE_FITH_

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <cstring>
#include "Eigen/Dense"
#include "spline.h"
using namespace std;

class LineFit
{
public:
	LineFit();
	~LineFit();
	void GetLinecoeff(const vector<double> &Vx, const vector<double> &Vy, int order);
	void fitCircle(const vector<double> &Vx, const vector<double> &Vy);
	void spineLine(const vector<double> &tmp_x, const vector<double> &tmp_y, std::vector<double> &outVecPointX, std::vector<double> &outVecPointY, const float &step);
	double GetCoeffA0() { return A0; }
	double GetCoeffA1() { return A1; }
	double GetCoeffA2() { return A2; }
	double GetCoeffA3() { return A3; }

	double GetCoeffB0() { return b0; }
	double GetCoeffB1() { return b1; }
	double ClaMutiA(std::vector<double> &Vx, int ex);
	double ClaMutiB(std::vector<double> &Vx,std::vector<double> &Vy, int ex);

private:
	void fitWithLagrange(std::vector<double> &X, std::vector<double> &Y, int orders);
	void FitterLeastSquareMethod(vector<double> &X, vector<double> &Y, int orders);
	double A0;
	double A1;
	double A2;
	double A3;
	double b0;
	double b1;
};

#endif /* _LINE_FITH_ */