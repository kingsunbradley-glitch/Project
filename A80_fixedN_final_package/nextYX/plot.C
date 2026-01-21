void plot() {
    // 1. 定义 Ntuple，结构要和列顺序对应
    TNtuple *nt = new TNtuple("nt", "Data from Excel", "D1:Q1:Q2:N");

    // 2. 读取文本文件
    // 如果是 CSV (逗号分隔)，通常 ReadFile 能自动识别，或者指定分隔符
    nt->ReadFile("data.txt"); 

    // 3. 绘图分析
    TCanvas *c1 = new TCanvas("c1", "Scan Results", 800, 600);
    
    // 技巧：只画出粒子数 N > 某个阈值的点，过滤掉无效数据
    nt->Draw("Q2:Q1:D1:N", "N > 10", "glcolz"); 
    
    // 或者画出 2D 投影：看看 D1 和 Q1 的关系，N作为颜色
    // new TCanvas();
    // nt->Draw("Q1:D1", "N", "colz");
}