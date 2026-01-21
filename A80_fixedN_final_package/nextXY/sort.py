import pandas as pd

def process_data(input_file, output_file):
    # 1. 读取数据
    try:
        df = pd.read_csv(input_file, sep='\s+')
    except FileNotFoundError:
        print(f"错误: 找不到文件 {input_file}")
        return

    # 定义自定义聚合逻辑
    def custom_agg(x):
        # 如果最大值和最小值的差在 1500 以内 (<= 1500)
        if (x.max() - x.min()) <= 1500:
            return x.mean()
        else:
            # 差值大于 1500，取最大值
            return x.max()

    # 2. 分组并应用自定义逻辑
    # 使用 apply 来应用上面的 custom_agg 函数
    df_processed = df.groupby(['D1', 'Q1', 'Q2'])['N'].apply(custom_agg).reset_index()

    # 3. 排序
    # 第一排降序，第一排相同第二排降序，第二排一样第三排降序
    df_sorted = df_processed.sort_values(by=['D1', 'Q1', 'Q2'], ascending=[False, False, False])

    # 4. 输出结果
    print("处理后的前20行数据预览:")
    print(df_sorted.head(20).to_string(index=False))
    
    # 保存到文件
    df_sorted.to_csv(output_file, sep='\t', index=False)
    print(f"\n完整结果已保存至: {output_file}")

if __name__ == "__main__":
    process_data('data.txt', 'processed_result.txt')