import uproot
import pandas as pd
import numpy as np

root_path = "summary.root"
tree_name = "summary"
out_xlsx  = "summary.xlsx"

with uproot.open(root_path) as f:
    t = f[tree_name]
    arr = t.arrays(library="np")  # 不用 pandas 适配层，避免 awkward-pandas 依赖

# arr 是 dict: {branch_name: numpy array}
df = pd.DataFrame()

for k, v in arr.items():
    v = np.asarray(v)
    # 标量列
    if v.ndim == 1:
        df[k] = v
    # 固定长度数组列，比如 muE[3] -> muE_0 muE_1 muE_2
    elif v.ndim == 2:
        for i in range(v.shape[1]):
            df[f"{k}_{i}"] = v[:, i]
    else:
        # 更高维/不规则的先跳过
        print(f"Skip {k}: ndim={v.ndim}")

df.to_excel(out_xlsx, index=False)
print("Wrote", out_xlsx, "rows=", len(df), "cols=", len(df.columns))
