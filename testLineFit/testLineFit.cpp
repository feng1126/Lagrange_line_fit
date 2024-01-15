// testLineFit.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "lineFit.h"
int main()
{

    double A0 = 2.3367111683;
    double A1 = 0.0273767686;
    double A2 = -0.0028188142;
    double A3 = -0.0000539880;

    std:vector<double > x, y;

    for (int i = 100; i < 20000; i++) 
    {

        double ram = rand() / (RAND_MAX * 1.f) * 10.0 ;

       // std::cout << ram << std::endl;

        double tmp = i / 10.0 ;
        x.push_back(tmp);

        double ytmp = A3 * tmp * tmp * tmp + A2 * tmp * tmp + A1 * tmp + A0 + ram;
        y.push_back(ytmp);
    }

    //A.block(orders + 1, 0, 1, orders + 1) = A.block(0, orders + 1,  orders + 1, 1).transpose();
//A.block(orders + 2, 0, 1, orders + 1) = A.block(0, orders + 2,  orders + 1, 1).transpose();


    LineFit pLineFit;
    pLineFit.GetLinecoeff(x, y, 3);

    printf("%.10f,  %.10f  ,%.10f , %.10f \n" ,pLineFit.GetCoeffA0(), pLineFit.GetCoeffA1(), pLineFit.GetCoeffA2(), pLineFit.GetCoeffA3());

    return 0;

    std::cout << "Hello World!\n";
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
