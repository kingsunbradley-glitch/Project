程序用来计算束流剂量和反应截面，使用的是python语言+root。执行命令为python3 XX.py

由于root在6.22版本才支持python，如需在自己目录下使用需要在主目录的.bashrc文件修改root的位置
将/usr/local/root-5.34/....修改为/root-6.22/....
修改后可能会对别的程序造成影响，也可以在目录/home/zhulin/public下使用

输入参数在Class InPut函数中，各项参数如下：
fit,fit_num-拟合选项及需要拟合的文件号
k,b-T2-FC或MACCT02和BPMSum的拟合参数
auto-自动计算每个文件流强
exp_num-实验编号
file_path-map文件路径
history_path-历史记录文件位置（一般不变）
charge-束流电荷量
thickness-靶厚
duty_ratio-占空比
eff-分别为差分传输效率，靶到DSSD探测器的传输效率，DSSD探测器的效率
target-靶的质量数
BPM_threshold,DumpThreshold,MACCT_threshold-BPMSum、Dump和MACCT02的阈值，需要根据map文件设置
skip-需要跳过的文件
reset-需要重新计算剂量的文件号
target_limit-靶位限制，根据屏蔽靶的数量设置，10块屏蔽1块就是0.9
runstart,runstop-计算的文件编号（runstop可以注释掉,会选取从runstart往后的所有文件）
D1_mag,TKE,central_energy_excitation_energy-D1磁刚度，束流能量，靶中心能量，激发能
drawset-画图选项

v1.0(auto)更新说明：
整合了按文件计算和按小时计算的代码（测试中），使得程序运行更快。

v0.7(auto)更新说明：
整合DrawDose.py进main.py,新增拟合选项，用来拟合BPMSum和MACCT02或T2FC(需要修改传输效率)。
新增函数备注

各项参数如下：
auto-自动计算每日流强
exp_num-实验编号
path-map文件路径
fit,fit_num-拟合选项及需要拟合的文件号
k,b-T2-FC或MACCT02和BPMSum的拟合参数
charge-束流电荷量
thickness-靶厚
duty_ratio-占空比
eff-分别为差分传输效率，靶到DSSD探测器的传输效率，DSSD探测器的效率
target-靶的质量数
BPM_threshold,DumpThreshold,MACCT_threshold-BPMSum、Dump和MACCT02的阈值，需要根据map文件设置
runlog_file-runlog.txt文件位置
history_path-历史记录文件位置
skip-需要跳过的文件
reset-需要重新计算剂量的文件号
runstart,runstop-计算的文件编号
total_dose-之前的总计束流剂量
drawset-画图选项
dose_start,dose_stop-文件范围

v0.6(auto)更新说明：
新增每小时平均流强，每小时累积剂量统计图

v0.5(auto)更新说明：
除需要修改以下信息外，可以直接运行test.py给出每天的束流剂量和截面
大部分实验信息（包括实验编号，文件路径，靶信息等等）以及160行reset重置束流剂量,161行skip跳过的文件号
手动计算则需要更改159行auto为False，再修改runstart,runstop,total_dose
删除cal==1的代码。
新增BPMSum和MACCT02异常检测，目前只判断了卡死的情况。
电流更改次数太多，并且不容易掉电，备注掉了电流判断。

实验备注：SS018直接计算靶前流强所以差分效率更改为1；82号文件开始3号靶20%屏蔽，4号靶30%屏蔽，dose计算时乘以系数(1-1/4*0.5),SS019实验已去掉。
SS020加速器拟合没做好，继续使用BPMSUM与T2FC拟合，差分效率更改为0.9，更换能量后重新计算靶前流强,SS020实验开始记录MACCT02。
修改：24.5.7将束流剂量按平均流强计算，减少重复计算。
SS024使用19块U靶，1块Yb靶，只计算了U的剂量，*0.95

history 每列信息
start	end	total_time	file	SS-ACCT	dose	total_dose	cross_section	1.84*cross_section	duty_ratio
