#include "Calibrator.h"
#include "TApplication.h"
#include "Config.h"

int main(int argc, char** argv) {
    // 即使是批处理，加载某些库时初始化 TApplication 也是个好习惯
    TApplication app("app", &argc, argv);

    Calibrator cal;
    cal.Run(); // 内部根据 Config.h 读取文件列表

    return 0;
}