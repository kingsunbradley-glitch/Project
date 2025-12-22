#include "Calibrator.h"
#include "TApplication.h"
#include <iostream>

int main(int argc, char** argv) {
    // TApplication 用于保证 ROOT 的图形库正常初始化，即使在批处理模式下
    TApplication app("app", &argc, argv);
    
    try {
        std::cout << "=== Step 2: DSSD Calibration Program ===" << std::endl;
        
        // 创建实例并运行
        Calibrator cal;
        cal.Run();
        
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}