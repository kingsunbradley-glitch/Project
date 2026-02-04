import glob
import os

def main():
    # 1. 查找当前目录下所有符合 *all_nuclei-new.dat 模式的文件
    input_files = glob.glob("*all_nuclei-new.dat")
    
    if not input_files:
        print("未找到符合 *all_nuclei-new.dat 的文件。")
        return

    print(f"找到 {len(input_files)} 个输入文件，开始处理...")

    for input_file in input_files:
        # 2. 提取模型名称作为输出文件名的一部分
        # 例如: UNEDF1_all_nuclei-new.dat -> Q@a_UNEDF1.dat
        base_name = input_file.replace("_all_nuclei-new.dat", "")
        output_file = f"Q@a_{base_name}.dat"
        
        print(f"正在转换: {input_file}  ->  {output_file}")
        
        count = 0
        try:
            with open(input_file, 'r', encoding='utf-8') as fin, \
                 open(output_file, 'w', encoding='utf-8') as fout:
                
                # 3. 写入表头
                fout.write("#  Z    N    A    Qalpha(MeV)\n")
                
                for line in fin:
                    line = line.strip()
                    # 跳过空行、注释行或标题行
                    if not line or line.startswith('#') or line.startswith('Symbol'):
                        continue
                    
                    parts = line.split()
                    
                    if len(parts) < 10:
                        continue
                    
                    # 提取关键数据
                    z_str = parts[1]
                    n_str = parts[2]
                    a_str = parts[3]
                    q_val = parts[9]
                    
                    # 4. 数据清洗：跳过 "No_Data"
                    if q_val == "No_Data":
                        continue
                    
                    # 5. 【新增条件】只输出 Z >= 104 的数据
                    try:
                        z_val = int(z_str)
                        if z_val < 104:
                            continue
                    except ValueError:
                        continue  # 如果Z不是数字，跳过
                    
                    # 6. 格式化输出: Z N A Q
                    fout.write(f"{z_str:>4} {n_str:>4} {a_str:>4} {q_val:>12}\n")
                    count += 1
            
            print(f"  - 已生成 {count} 条数据 (Z >= 104)")
                    
        except Exception as e:
            print(f"处理文件 {input_file} 时出错: {e}")

    print("所有文件处理完成。")

if __name__ == "__main__":
    main()