这个文件夹就是你说的“function 目录”：所有 .cpp/.h + Makefile + mulicore_run.sh 都在同一层。

常用：
1) 手动编译：
   make -j

2) 并行跑（推荐）：
   # 假设 function 在数据目录的子目录里，数据在上一层：
   ./mulicore_run.sh 8 ../ 1 385

输出默认在 function/out/ 里：
  out/summary.root
  out/qa_plots.root
  out/pdf/
  out/jobs/
