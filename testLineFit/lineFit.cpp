/*
 * @Author: niew nie_w@reachauto.com
 * @Date: 2023-04-20 15:09:32
 * @LastEditors: niew nie_w@reachauto.com
 * @LastEditTime: 2024-01-12 16:36:56
 * @FilePath: \roadFusion\includes\lineFit.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A
 */

#include "lineFit.h"
LineFit::LineFit()
{
	A0 = 0;
	A1 = 0;
	A2 = 0;
	A3 = 0;
}

void LineFit::GetLinecoeff(const vector<double> &Vx, const vector<double> &Vy, int order)
{
	A0 = 0;
	A1 = 0;
	A2 = 0;
	A3 = 0;
	b0 = 0;
	b1 = 0;
	std::vector<double> outVecPointX, outVecPointY;
	outVecPointX.clear();
	outVecPointY.clear();
	if (Vx.size() < 4 && Vx.size() > 1)
	{
		tk::spline s;
		// s.set_boundary(tk::spline::first_deriv, 0.0,
		//	tk::spline::first_deriv, 0.0);
		s.set_points(Vx, Vy, tk::spline::linear);

		for (double n = Vx[0]; n < Vx[Vx.size() - 1]; n = n + 0.1)
		{
			double x = n;
			outVecPointX.push_back(x);
			outVecPointY.push_back(s(x));
		}

		double x = Vx[Vx.size() - 1];
		outVecPointX.push_back(Vx[Vx.size() - 1]);
		outVecPointY.push_back(Vy[Vy.size() - 1]);
	}
	else
	{

		outVecPointX = Vx;
		outVecPointY = Vy;
	}

	if (outVecPointX.size() > 3)
	{
		fitWithLagrange(outVecPointX, outVecPointY, order);
	}
}

LineFit::~LineFit()
{
}

double LineFit::ClaMutiA(std::vector<double> &Vx, int ex)
{
	double ReSum = 0;
	for (int i = 0; i < Vx.size(); i++)
	{
		ReSum += pow(Vx[i] - Vx[0], ex);
	}
	return ReSum;
}

double LineFit::ClaMutiB(std::vector<double> &Vx, std::vector<double> &Vy, int ex)
{
	double dReMultiSum = 0;
	for (int i = 0; i < Vx.size(); i++)
	{
		dReMultiSum += pow(Vx[i] - Vx[0], ex) * (Vy[i] - Vy[0]);
	}
	return dReMultiSum;
}

void LineFit::FitterLeastSquareMethod(vector<double> &X, vector<double> &Y, int orders)
{
	Eigen::Map<Eigen::VectorXd> xvals(X.data(), X.size());
	Eigen::Map<Eigen::VectorXd> yvals(Y.data(), Y.size());
	Eigen::MatrixXd A(xvals.size(), orders + 1);
	for (int i = 0; i < xvals.size(); i++)
	{
		A(i, 0) = 1.0;
	}
	for (int j = 0; j < xvals.size(); j++)
	{
		for (int i = 0; i < orders; i++)
		{
			A(j, i + 1) = A(j, i) * xvals(j);
		}
	}
	Eigen::VectorXd result = A.householderQr().solve(yvals);
	if (result.size() == 4)
	{
		A0 = result[0];
		A1 = result[1];
		A2 = result[2];
		A3 = result[3];
	}
}

void LineFit::fitWithLagrange(std::vector<double> &X, std::vector<double> &Y, int orders)
{

	int size = orders + 3;

	Eigen::VectorXd b(size, 1);

	b.setZero();

	Eigen::MatrixXd A(size, size);

	A.setZero();

	for (int i = 0; i < orders + 1; i++)
	{
		A(i, orders + 1) = pow(X[0] - X[0], i);
		A(i, orders + 2) = pow(X[X.size() - 1] - X[0], i);
		b[i] = ClaMutiB(X, Y, i);
	}

	A(0, 0) = X.size();
	for (int i = 1; i < orders + 1; i++)
	{
		A(0, i) = ClaMutiA(X, i);
	}

	for (int i = 0; i < orders; i++)
	{
		A(1, i) = A(0, i + 1);
	}
	A(1, orders) = ClaMutiA(X, 4);

	for (int i = 0; i < orders; i++)
	{
		A(2, i) = A(1, i + 1);
	}
	A(2, orders) = ClaMutiA(X, 5);

	for (int i = 0; i < orders; i++)
	{
		A(3, i) = A(2, i + 1);
	}
	A(3, orders) = ClaMutiA(X, 6);



	b[orders + 1] = Y[0] - Y[0];
	b[orders + 2] = Y[Y.size() - 1] - Y[0];

	//for (int j = 0; j < orders + 1; j++)
	//{
	//	A(orders + 1, j) = pow(X[0] - X[0], j);
	//}
	//for (int j = 0; j < orders + 1; j++)
	//{
	//	A(orders + 2, j) = pow(X[X.size() - 1] - X[0], j);
	//}

	//std::cout << A << std::endl;


	A.block(orders + 1, 0, 1, orders + 1) = A.block(0, orders + 1,  orders + 1, 1).transpose();
	A.block(orders + 2, 0, 1, orders + 1) = A.block(0, orders + 2,  orders + 1, 1).transpose();


	std::cout << A << std::endl;

	Eigen::VectorXd result = A.householderQr().solve(b);



	double coeffA, coeffB, coeffC, coeffD = 0;
	if (result.size() == orders + 3)
	{
		coeffA = result[0];
		coeffB = result[1];
		coeffC = result[2];
		if (orders == 3)
		{
			coeffD = result[3];
		}

		double &a0 = X[0];
		A0 = coeffA - coeffD * a0 * a0 * a0 + coeffC * a0 * a0 - coeffB * a0 + Y[0];
		A1 = 3 * coeffD * a0 * a0 - 2 * coeffC * a0 + coeffB;
		A2 = -3 * coeffD * a0 + coeffC;
		A3 = coeffD;
	}
}


void LineFit::fitCircle(const vector<double> &Vx, const vector<double> &Vy)
{

	A0 = 0;
	A1 = 0;
	A2 = 0;
	A3 = 0;
	b0 = 0;
	b1 = 0;
	Eigen::MatrixXd A(Vx.size(), 3);
	Eigen::VectorXd B(Vx.size());
	for (int i = 0; i < Vx.size(); i++)
	{
		A(i, 0) = Vx[i];
		A(i, 1) = Vy[i];
		A(i, 2) = 1;
		B[i] = Vx[i] * Vx[i] + Vy[i] * Vy[i];
	}

	Eigen::VectorXd result = A.householderQr().solve(B);
	if (result.size() == 3)
	{
		A0 = result[0] / 2.0;
		A1 = result[1] / 2.0;
		A2 = std::sqrt(A0 * A0 + A1 * A1 - result[1]);
	}
}


void LineFit::spineLine(const vector<double> &tmp_x, const vector<double> &tmp_y, std::vector<double> &outVecPointX, std::vector<double> &outVecPointY, const float &step)
{

	if (tmp_x.size() < 1)
		return;

	outVecPointX.clear();
	outVecPointY.clear();
	tk::spline s;
	// s.set_boundary(tk::spline::first_deriv, 0.0,
	//	tk::spline::first_deriv, 0.0);
	s.set_points(tmp_x, tmp_y, tk::spline::linear);
	// printf("fitWithSpline  tmp_x1 0 %d %f %f \n", tmp_x.size(), tmp_x[0], tmp_x[tmp_x.size() - 1]);
	for (double n = tmp_x[0]; n < tmp_x[tmp_x.size() - 1]; n = n + step)
	{
		double x = n;
		outVecPointX.push_back(x);
		outVecPointY.push_back(s(x));
	}
	if (outVecPointX.size() > 0 && tmp_x.size() > 0)
	{
		if (tmp_x[tmp_x.size() - 1] - outVecPointX[outVecPointX.size() - 1] > 0.5 * step)
		{
			outVecPointX.push_back(tmp_x[tmp_x.size() - 1]);
			outVecPointY.push_back(tmp_y[tmp_y.size() - 1]);
		}
		else
		{
			outVecPointX[outVecPointX.size() - 1] = tmp_x[tmp_x.size() - 1];
			outVecPointY[outVecPointY.size() - 1] = tmp_y[tmp_y.size() - 1];
		}
	}
}