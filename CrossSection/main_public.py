import ROOT as root
import os
import numpy as np
import time
from calendar import month_abbr
import calendar
import matplotlib as mpl
import matplotlib.pyplot as plt

# v1.1 - Fixed & Updated

class Input():
    def __init__(self):
        # Fitting Parameters
        self.fit = False
        self.fit_num = "all"
        self.k = 1.0615 # MACCT = k * BPMSum + b
        self.b = 1.9582
        
        # Experiment Information
        self.auto = True
        self.exp_num = "SS032"
        self.file_path = os.path.abspath("/home/evalie2/Project/document/273Ds/inter_map/")
        
        # --- FIX 2: 自动创建 history 文件夹，防止报错 ---
        self.history_path = os.getcwd() + "/history"
        if not os.path.exists(self.history_path):
            os.makedirs(self.history_path)
            print(f"已创建文件夹: {self.history_path}")
        # ---------------------------------------------

        self.charge = 13 # Beam charge (电荷态)
        # self.thickness = 0.469424  #1
        self.thickness = 0.604  # 2              # Target thickness in mg/cm^2  
        self.duty_ratio = 0.91 # Duty ratio 
        self.eff = [1, 0.52, 1] # Transmission efficiency  差分段传输效率，谱仪传输效率，探测器探测效率
        self.target = 238 # Target atomic mass
        self.BPM_threshold = 3 # BPM threshold
        self.MACCT_threshold = 0 # MACCT threshold
        self.reset = [ ] 
        
        # --- FIX 1: 确保 self.skip 是列表 (List) ---
        # 如果只想跳过一个文件，也要写成 [22]
        self.skip = list(range(28, 35))  # Skip file numbers
        # self.skip = [22] 
        # -----------------------------------------
        
        self.target_limit = 1   # 靶位限制
        # self.runstart = 2  #1
        self.runstart = 21  # 2  # Start file number
        # self.runstop = 81  #1
        self.runstop = 50   # 2  # Stop file number
        # self.D1_mag =2.21  #1
        self.D1_mag = 2.212 # 2   # 磁刚度
        self.TKE = 225 # Beam kinetic energy in MeV
        # self.central_energy = 212.7    #1
        self.central_energy = 212.2 # 2  # Central energy in MeV
        # self.excitation_energy = 49.5 #1
        self.excitation_energy = 49.1 # 2    # Excitation energy in MeV
        
        self.file_history, self.perhour_history, self.total_time = self.read_history()
        self.num_list = self.read_numder_list() 
        
        # Paint Parameters
        self.drawset = True

    def read_numder_list(self):
        result = []
        if not os.path.exists(self.file_path):
            print(f"错误: 找不到数据文件夹 {self.file_path}")
            return []
            
        for filename in os.listdir(self.file_path):
            if filename.startswith(self.exp_num) and filename.endswith("map.root"):
                try:
                    num = int(filename[5:10])
                    if num >= self.runstart and num <= (self.runstop if "runstop" in self.__dict__ else float("inf")):
                        result.append(num)
                except ValueError:
                    continue
        result.sort()
        return result

    def read_history(self):
        result1, result2, result3 = [], np.array([]), 0
        
        # 读取主历史文件
        hist_file = self.history_path + "/{}.dat".format(self.exp_num)
        if os.path.exists(hist_file):
            with open(hist_file, "r") as File:
                for lines in File.readlines():
                    tmp = lines.strip().split("\t")
                    if len(tmp) > 2:
                        result1.append(tmp)
                        try:
                            result3 += int(tmp[2].split(" s")[0])
                        except:
                            pass
                            
        # 读取每小时历史文件
        hour_file = self.history_path + "/{}_perhour.dat".format(self.exp_num)
        if os.path.exists(hour_file):
            try:
                result2 = np.loadtxt(hour_file, dtype=float, usecols=(0, 1, 2))
                if result2.ndim == 1 and len(result2) > 0: # 处理只有一行数据的情况
                    result2 = np.array([result2])
            except:
                result2 = np.array([])
                
        return result1, result2, result3
    
class StatisticPV(Input):
    def __init__(self):
        super().__init__()
        self.year = int(time.asctime().split()[-1])
        
    def main(self):
        # 获取上次处理到的文件号，如果是新运行则从头开始
        last_file_num = int(self.file_history[-1][3].split("-")[0]) if self.file_history else self.runstart - 1
        last_file_dose = float(self.file_history[-1][6]) if self.file_history else 0
        
        tmp_MACCT, flag = np.array([]), False
        
        # 恢复绘图数据
        if len(self.perhour_history) != 0:
            pre_time = self.perhour_history[-1][0]
            # 处理 reshape 问题，确保维度正确
            if self.perhour_history.ndim == 1:
                systime = np.array([self.perhour_history[0]])
                MACCT_perhour = np.array([self.perhour_history[1]])
                total_dose_perhour = np.array([self.perhour_history[2]])
            else:
                systime = self.perhour_history[:,0]
                MACCT_perhour = self.perhour_history[:,1]
                total_dose_perhour = self.perhour_history[:,2]
        else:
            pre_time, systime, MACCT_perhour, total_dose_perhour = 0, np.array([]), np.array([0]), np.array([0])
            
        for num in self.num_list:
            total_time = 0
            
            # --- 逻辑判断：是否跳过文件 ---
            if num <= last_file_num or num in self.skip:
                continue
            # --------------------------
            
            if num in self.reset:
                last_file_dose = 0
                flag = True
                
            print("正在读取{0}{1:0>5d}_map.root".format(self.exp_num, num))
            
            try:
                f = root.TFile(self.file_path + "/{0}{1:0>5d}_map.root".format(self.exp_num, num))
                if f.IsZombie():
                    print(f"警告: 文件 {num} 损坏或无法读取")
                    continue
                trPV = f["tr_PV"]
            except:
                print(f"读取文件 {num} 失败")
                continue

            if num == self.runstart:
                trPV.GetEntry(0)
                systime = np.append(systime, trPV.SysTime)
                pre_time = trPV.SysTime
            
            entries = trPV.GetEntries()
            MACCT, dose = np.array([]), 0
            start_time = ""
            end_time = ""

            for entry in range(entries):
                trPV.GetEntry(entry)
                if entry == 0:
                    start_time = self.transform_time(trPV.SysTime)
                if entry == entries - 1:
                    end_time = self.transform_time(trPV.SysTime)
                
                # 处理时间间隔过大的情况 (补零)
                if trPV.SysTime - pre_time > 3600:
                    while trPV.SysTime - pre_time >= 7200:
                        pre_time += 3600
                        MACCT_perhour = np.append(MACCT_perhour, 0)
                        total_dose_perhour = np.append(total_dose_perhour, total_dose_perhour[-1])
                        systime = np.append(systime, pre_time)
                    
                    MACCT_ave = np.mean(tmp_MACCT) if len(tmp_MACCT) != 0 else 0
                    tmp_MACCT = np.array([])
                    MACCT_perhour = np.append(MACCT_perhour, MACCT_ave)
                    
                    if flag:
                        tmp_dose = 0
                        flag = False
                    else:
                        tmp_dose = total_dose_perhour[-1]
                        
                    dose_perhour = (
                        tmp_dose
                        + MACCT_ave
                        * 10**-6
                        / (self.charge * 1.6 * 10**-19)
                        * (trPV.SysTime - pre_time)
                        * self.duty_ratio
                        * self.eff[0]
                        * self.target_limit
                    )
                    total_dose_perhour = np.append(total_dose_perhour, dose_perhour)
                    pre_time = trPV.SysTime
                    systime = np.append(systime, pre_time)
                
                # 有效束流判断
                if self.conditions([trPV.BPMSum, trPV.Dump, [trPV.Q1, trPV.Q2, trPV.Q3, trPV.D1, trPV.D2]]):
                    total_time += 10
                    if self.fit and (self.fit_num == "all" or num in self.fit_num):
                        tmp = trPV.BPMSum * self.k + self.b
                    else:
                        tmp = trPV.MACCT02
                    
                    MACCT = np.append(MACCT, max(tmp,0))
                    tmp_MACCT = np.append(tmp_MACCT, max(tmp,0))
                    
                    # 剂量累加
                    dose += (
                        tmp
                        * 10**-6
                        / (self.charge * 1.6 * 10**-19)
                        * 10
                        * self.duty_ratio
                        * self.eff[0]
                        * self.target_limit
                    )
                else:
                    tmp_MACCT = np.append(tmp_MACCT, 0)
            
            last_file_dose += dose
            
            # 计算反应截面
            if last_file_dose > 0:
                cross_section = (
                    1
                    / (
                        last_file_dose
                        * self.thickness
                        * 10**-3
                        / self.target
                        * 6.022
                        * 10**23
                        * self.eff[1]
                        * self.eff[2]
                    )
                    * 10**36
                )
            else:
                cross_section = 0
                
            f.Close()
            
            # 计算文件平均 euA
            file_MACCT_ave = np.mean(MACCT) if len(MACCT) != 0 else 0
            
            # --- FIX 3: 打印时转换为 puA 单位 ---
            # puA = euA * 效率 / 电荷态
            current_puA = (file_MACCT_ave * self.eff[0]) / self.charge
            # --------------------------------
            
            print(
                    "{0}——{1}束流总时间:{2} s({3:.2f} h),文件号{4}{5:0>5d}.root-{4}{6:0>5d}.root, SSACCT流强(puA):{7:.2f}, 束流剂量:{8:.3e},累计束流剂量:{9:.3e},1个事件反应截面:{10:.3f} pb".format(
                        start_time,
                        end_time,
                        total_time,
                        total_time / 3600,
                        self.exp_num,
                        num,
                        num,
                        current_puA, # 这里显示 puA
                        dose,
                        last_file_dose,
                        cross_section
                    )
                )
            self.total_time += total_time
            if self.auto:
                self.write_file([start_time, end_time, total_time, num, file_MACCT_ave, dose, last_file_dose, cross_section])
        
        # 循环结束
        print("有效束流总时间{0:.2f} h".format(self.total_time/3600))
        
        # --- FIX 4: 计算全实验平均 puA ---
        if self.total_time > 0:
            # I_avg(puA) = Total_Q / Total_t
            # last_file_dose 是总粒子数, 乘 e 得到总库仑量, 除以时间(s)得到安培, 再乘 1e6 得到 uA
            total_avg_puA = (last_file_dose * 1.60217663e-19) / (self.total_time * 1e-6)
            print("-" * 50)
            print("【统计结果】")
            print("实验全过程总平均流强: {0:.2f} puA".format(total_avg_puA))
            print("累计总粒子数 (Dose): {0:.3e}".format(last_file_dose))
            print("-" * 50)
        # ------------------------------

        if self.drawset:
            self.Draw(MACCT_perhour, total_dose_perhour, systime)
            self.write_perhour(MACCT_perhour, total_dose_perhour, systime)
    
    def transform_time(self, systime):
        timearray = time.localtime(systime)
        return "{0}月{1}日 {2:0>2d}:{3:0>2d}".format(timearray.tm_mon, timearray.tm_mday, timearray.tm_hour, timearray.tm_min)
    
    def conditions(self, state):
        BPM, DUMP, currentlis = state[0], state[1], state[2]
        result = True
        if BPM < self.BPM_threshold:
            result = False
        return result
    
    def write_file(self, output):
        # 这里的 output[4] 依然保持 euA 写入文件，保持历史一致性
        with open(self.history_path + "/{}.dat".format(self.exp_num), "a") as File:
            File.write(
                "{0}\t{1}\t{2} s({3:.2f} h)\t{4}-{5}\t{6:.2f}\t{7:.3e}\t{8:.3e}\t{9:.3f}\t{10}\n".format(
                    output[0],
                    output[1],
                    output[2],
                    output[2] / 3600,
                    output[3],
                    output[3],
                    output[4], # 存入文件的是 file_MACCT_ave (euA)
                    output[5],
                    output[6],
                    output[7],
                    self.duty_ratio
                )
            )
        with open(self.history_path + "/{}_wiki.dat".format(self.exp_num), "a") as File:
            File.write("|-\n")
            File.write(
                "| '''{0} -- {1}''' || {2} s({3:.2f} h) || {4}{5:0>5d}.root-{4}{6:0>5d}.root || 3n || {7} || {8} || {9} || {10:.3e} || {11:.3e} || {12:.3f} || {13} \n".format(
                    output[0],
                    output[1],
                    output[2],
                    output[2] / 3600,
                    self.exp_num,
                    output[3],
                    output[3],
                    self.TKE,
                    self.central_energy,
                    self.excitation_energy,
                    output[5],
                    output[6],
                    output[7],
                    self.D1_mag
                )
            )
            
    def Draw(self, MACCT, dose, systime):
        if len(systime) == 0:
            return
            
        SSACCT = MACCT * self.eff[0] / self.charge
        
        # 修复 Draw 函数中的文件读取，防止找不到文件
        try:
            f = root.TFile(self.file_path + "/{0}{1:0>5d}_map.root".format(self.exp_num, self.runstart))
            if f.IsZombie():
                # 如果 runstart 文件坏了，尝试直接用 systime[0]
                timearray = time.localtime(systime[0])
            else:
                tmp_PV = f["tr_PV"]
                tmp_PV.GetEntry(0)
                timearray = time.localtime(tmp_PV.SysTime)
                f.Close()
        except:
             timearray = time.localtime(systime[0])

        month, day, h, m, s = timearray.tm_mon, timearray.tm_mday, timearray.tm_hour, timearray.tm_min, timearray.tm_sec
        first_time = systime[0] - int(h) * 3600 - int(m) * 60 - int(s)
        my_ticks, my_labels=[first_time], [list(calendar.month_abbr)[month] + " " + str(day)]
        interval = len(systime) // 480 + 1
        
        loop_time = first_time
        while loop_time <= systime[-1]:
            loop_time += 24 * 3600 * interval
            my_ticks.append(loop_time)
            
        for i in range(1, len(my_ticks)):
            my_labels.append(time.strftime("%b %d", time.localtime(my_ticks[i])).replace(" 0"," "))
            
        tmp = my_ticks[0]
        normalize_time = systime - tmp
        my_ticks = [x - tmp for x in my_ticks]
        
        mpl.use("tkagg")
        fig, (ax1, ax2)=plt.subplots(2, 1)
        ax1.step(normalize_time, SSACCT, where="post", color="black")
        ax1.set_ylim(0, 1.05 * max(SSACCT) if len(SSACCT)>0 else 1)
        ax1.set_xlim(0, normalize_time[-1])
        if len(SSACCT) > 0:
            ax1.vlines(normalize_time[-1] + 3600, 0, SSACCT[-1], color="black")
            ax1.hlines(SSACCT[-1], normalize_time[-1], normalize_time[-1] + 3600, color="black")
            
        ax1.set_xticks(my_ticks)
        ax1.set_xticklabels(my_labels)
        ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(6*3600*interval))
        ax1.tick_params(which="major", length=7)
        ax1.tick_params(which="minor", length=5)
        ax1.set_ylabel(r"SS-ACCT (p$\mu$A)")
        
        ax2.step(normalize_time, dose / 1e18, where="post", color="black")
        ax2.set_ylim(0, 1.05 * max(dose / 1e18) if len(dose)>0 else 1)
        ax2.set_xlim(0, normalize_time[-1])
        if len(dose) > 0:
            ax2.vlines(normalize_time[-1] + 3600, 0, dose[-1] / 1e18, color="black")
            ax2.hlines(dose[-1] / 1e18, normalize_time[-1], normalize_time[-1] + 3600, color="black")
            
        ax2.set_xticks(my_ticks)
        ax2.set_xticklabels(my_labels)
        ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(6*3600*interval))
        ax2.tick_params(which="major", length=7)
        ax2.tick_params(which="minor", length=5)
        ax2.set_xlabel("Date")
        ax2.set_ylabel(r"Total Beam Dose ($\times$10$^{18}$ ions)")
        plt.show()
        
    def write_perhour(self, MACCT, total_dose, systime):
        with open(self.history_path + "/{}_perhour.dat".format(self.exp_num), "w") as F:
            for i in range(len(MACCT)):
                F.write(
                        "{0}\t{1:.2f}\t{2:.5e}\t{3}\n".format(
                            systime[i],
                            MACCT[i],
                            total_dose[i],
                            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(systime[i]))
                            )
                        )
        
if __name__ == "__main__":
    Sta = StatisticPV()
    Sta.main()