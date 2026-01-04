// 辅助函数：将直方图转换为旋转后的 TGraph (模拟 hist 效果)
TGraph* CreateRotatedGraph(TH1* h) {
    if (!h) return nullptr;
    
    TGraph* gr = new TGraph();
    int nbins = h->GetNbinsX();
    int pointIndex = 0;
    
    for (int i = 1; i <= nbins; ++i) {
        double content = h->GetBinContent(i);
        double edgeLow = h->GetBinLowEdge(i);
        double edgeHigh = h->GetBinLowEdge(i+1);
        
        // 关键：交换 X 和 Y
        // 正常是 (edge, content)，这里变成 (content, edge)
        // 为了画出阶梯状 (Step)，每个 bin 需要两个点
        gr->SetPoint(pointIndex++, content, edgeLow);
        gr->SetPoint(pointIndex++, content, edgeHigh);
    }
    
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    return gr;
}

void draw_run00359() {
    // 1. 打开文件
    TFile *file = TFile::Open("qa_plots.root");
    if (!file || file->IsZombie()) {
        printf("错误: 无法打开文件 qa_plots.root\n");
        return;
    }
    TDirectory *dir = (TDirectory*)file->Get("run00324");
    if (!dir) dir = gDirectory; 

    // 2. 创建画板
    TCanvas *c = new TCanvas("c_rot", "Run 00359 Rotated Line", 900, 900);
    c->Divide(2, 2); 
    
    Color_t lineColor = kBlue;

    // --- Pad 1 (1,1): XY 主图 ---
    c->cd(1);
    gPad->SetRightMargin(0.12);
    TH2 *hXY = (TH2*)dir->Get("hXY_all");
    if (hXY) {
        hXY->SetStats(0);
        hXY->SetTitle("XY Distribution");
        hXY->Draw("colz");
    }

    // --- Pad 2 (1,2): Y 投影 (旋转的线条) ---
    c->cd(2);
    TH1 *hY = (TH1*)dir->Get("hY_all");
    if (hY) {
        // 准备工作
        hY->SetLineColor(lineColor);
        
        // 1. 创建旋转后的 Graph
        TGraph *grY = CreateRotatedGraph(hY);
        
        // 2. 这里的坐标轴需要手动设置，因为我们是在画 Graph 而不是 Hist
        //    X轴现在表示计数 (Counts)，Y轴表示通道 (Channel)
        double maxCounts = hY->GetMaximum() * 1.1; // 稍微留点空隙
        double yMin = hY->GetXaxis()->GetXmin();
        double yMax = hY->GetXaxis()->GetXmax();
        
        // 3. 画一个空的坐标系框架 (Frame)
        //    参数: X_min, Y_min, X_max, Y_max
        TH1F *frame = c->cd(2)->DrawFrame(0, yMin, maxCounts, yMax);
        frame->SetTitle("Y Projection;Counts;DSSD Y Channel");
        
        // 4. 画出我们的折线
        if (grY) grY->Draw("L same"); // "L" 表示只画线
    }

    // --- Pad 3 (2,1): X 投影 (标准线条) ---
    c->cd(3);
    TH1 *hX = (TH1*)dir->Get("hX_all");
    if (hX) {
        hX->SetStats(0);
        hX->SetTitle("X Projection");
        hX->SetLineColor(lineColor);
        hX->Draw("hist"); // 标准 hist，没有柱子，只有线
    }

    // --- Pad 4 (2,2): 能量谱 ---
    c->cd(4);
    TH1 *hE = (TH1*)dir->Get("hE_5000_6000");
    if (hE) {
        hE->SetStats(0);
        hE->Draw("hist");
    }

    c->Update();
}