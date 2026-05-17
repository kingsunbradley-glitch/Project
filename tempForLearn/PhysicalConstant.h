#pragma once // 防止头文件被重复包含
#include "PhysicalConstantPlus.h" // for PhysicalConstant::pi

// 1. 声明自定义命名空间
namespace PhysicalConstant {
    // 2. 在命名空间中声明常量
    const double h=1.0545718e-34; // 约化普朗克常数，单位为焦耳·秒
    //const double pi=3.14159265358979323846; // 圆周率
    const double c=299792458; // 真空中的光速，单位为米

    double division(double x, double y);

    const double hbar = division(h, 2 * pi); // 约化普朗克常数，单位为焦耳·秒

}