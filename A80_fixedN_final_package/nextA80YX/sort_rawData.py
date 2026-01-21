import pandas as pd

def process_data(input_file: str, output_file: str):
    # rawData.txt is TSV (Tab-separated) because headers contain spaces
    df = pd.read_csv(input_file, sep="\t")

    # type conversion
    for c in ["A80eq_fixN_2", "D1 (T.m)", "Q1 (T/m)", "Q2 (T/m)"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["runnum"] = pd.to_numeric(df["runnum"], errors="coerce").astype("Int64")

    group_cols = ["D1 (T.m)", "Q1 (T/m)", "Q2 (T/m)"]

    def agg_group(g: pd.DataFrame) -> pd.Series:
        # treat non-positive (<=0) as invalid (e.g., -999)
        valid_mask = g["A80eq_fixN_2"] > 0

        if valid_mask.any():
            maxv = g.loc[valid_mask, "A80eq_fixN_2"].max()
            # keep only values within 50 of the max (inclusive)
            in_mask = valid_mask & (g["A80eq_fixN_2"] >= (maxv - 50))
            avg = g.loc[in_mask, "A80eq_fixN_2"].mean()
        else:
            # no valid values in this group
            in_mask = pd.Series(False, index=g.index)
            avg = -999.0

        in_runs = g.loc[in_mask, "runnum"].dropna().astype(int).tolist()
        out_runs = g.loc[~in_mask, "runnum"].dropna().astype(int).tolist()

        in_runs.sort()
        out_runs.sort()

        return pd.Series({
            "A80eq_fixN_2": avg,
            "runnum_in_avg": ",".join(map(str, in_runs)),
            "runnum_outside_50": ",".join(map(str, out_runs)),
        })

    df_processed = df.groupby(group_cols, dropna=False).apply(agg_group).reset_index()

    # sort: D1 desc, then Q1 desc, then Q2 desc (like your sort.py)
    df_sorted = df_processed.sort_values(by=group_cols, ascending=[False, False, False])

    print("Preview (top 20 rows):")
    print(df_sorted.head(20).to_string(index=False))

    df_sorted.to_csv(output_file, sep="\t", index=False)
    print(f"\nSaved to: {output_file}")

if __name__ == "__main__":
    process_data("rawData.txt", "sortedData.txt")
