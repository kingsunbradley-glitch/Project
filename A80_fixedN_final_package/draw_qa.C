void draw_qa() {
    // 1. 打开 ROOT 文件
    TFile *file = TFile::Open("qa_plots.root");
    if (!file || file->IsZombie()) {
        printf("错误: 无法打开文件 qa_plots.root\n");
        return;
    }

    // 2. 获取 run00176 目录 (分支)
    // 注意：如果 run00176 是由 TTree::Branch 创建的，逻辑不同。
    // 但看截图这些是直方图，run00176 应该是一个 TDirectory (文件夹)。
    TDirectory *dir = (TDirectory*)file->Get("run00115");
    
    // 如果找不到该目录，尝试直接在顶层查找（防止路径理解偏差）
    if (!dir) {
        printf("提示: 未找到 run00176 目录，尝试在文件顶层查找...\n");
        dir = gDirectory; 
    }

    // 3. 创建画板 (Canvas)
    // 宽 1200, 高 900
    TCanvas *c1 = new TCanvas("c1", "QA Plots Layout", 1200, 900);
    
    // 分割画板：4列 (Columns) x 3行 (Rows)
    c1->Divide(4, 3); 

    // 4. 定义后缀数组，对应 all, p1, p2, p3
    const char* suffixes[] = {"all", "p1", "p2", "p3"};
    
    // 5. 循环绘制
    // i=0 -> all, i=1 -> p1, i=2 -> p2, i=3 -> p3
    for (int i = 0; i < 4; ++i) {
        // --- 第一行：绘制 hXY (2D图) ---
        //位置计算: 第1行的pad编号是 1, 2, 3, 4
        c1->cd(i + 1); 
        TH2 *hXY = (TH2*)dir->Get(Form("hXY_%s", suffixes[i]));
        if (hXY) {
            hXY->SetStats(0); // 可选：去掉统计框让画面干净
            hXY->Draw("colz"); // 二维图推荐用 colz 选项
        } else {
            printf("丢失直方图: hXY_%s\n", suffixes[i]);
        }

        // --- 第二行：绘制 hX (1D图) ---
        //位置计算: 第2行的pad编号是 5, 6, 7, 8 (即上一行编号 + 4)
        c1->cd(i + 1 + 4); 
        TH1 *hX = (TH1*)dir->Get(Form("hX_%s", suffixes[i]));
        if (hX) hX->Draw();

        // --- 第三行：绘制 hY (1D图) ---
        //位置计算: 第3行的pad编号是 9, 10, 11, 12 (即上一行编号 + 8)
        c1->cd(i + 1 + 8); 
        TH1 *hY = (TH1*)dir->Get(Form("hY_%s", suffixes[i]));
        if (hY) hY->Draw();
    }

    // 更新画板显示
    c1->Update();
    
    // 如果需要保存图片，取消下面这行的注释
    // c1->SaveAs("qa_layout.png");
}