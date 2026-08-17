import pandas as pd
import numpy as np
from pathlib import Path

HERE = Path(__file__).parent # Expression_analysis Folder

# Load and Read the data file into DataFrame and inspect
def data_load_inspection():
    counts = pd.read_csv(HERE/"counts.csv")
    meta = pd.read_csv(HERE/"metadata.csv")

    print(counts.head())
    print("counts shape:",counts.shape)
    print("columns:",list(counts.columns))
    print("dtype:\n", counts.dtypes)

    counts.info()

    print(counts.describe())

    print("metadata:\n", meta.tail())

    return counts, meta


# Normalise the data with numpy
def normalise_the_data(counts):
    genes = counts["gene"].to_numpy()
    sample_cols = [c for c in counts.columns if c != "gene"]
    mat =counts[sample_cols].to_numpy(dtype=float, copy=True)

    print("matrix shape:",mat.shape, "dtype:", mat.dtype)

    lib_size = mat.sum(axis=0)

    print("library sizes per sample:",lib_size.astype(int))

    #Counts-per-million: divide each column by its total
    cpm = mat/lib_size * 1e6
    log_cpm = np.log1p(cpm)

    #boolean mask of "high" value and labeling
    threshold = 11.0
    mask = log_cpm > threshold
    print("high-expression cells:", int(mask.sum()),"of", mask.size)
    labels = np.where(mask,"high","low")

    #Selecting the most expressed(argmax DOWN the row)
    top_idx = np.argmax(mat, axis=0)

    for sample, gi in zip(sample_cols, top_idx):
        print(f" {sample}: top gene = {genes[gi]}")

    print("transposed shape(sample x genes):", mat.T.shape)
    print("flattened legth:", mat.ravel().shape[0])
    print("sorted unique library sizes:", np.unique(np.sort(lib_size)).astype(int))
    print("samples ordered smallest->largest library:",[sample_cols[i] for i in np.argsort(lib_size)])

    return log_cpm, labels


#Reshaping the wide to long & back
def reshape_long(counts):
    long =counts.melt(id_vars="gene", var_name="sample", value_name="count")
    print(long.head())
    print("long shape:", long.shape)

    wide_again = long.pivot_table(index="gene",columns="sample", values="count")
    print("pivot back shape:",wide_again.shape)
    print(wide_again.head(5))

    return long


#Add the metadata and fixing the missing values
def merge_and_clean(long, meta):
    merged = long.merge(meta,on="sample", how ="left")
    print("missing values per column:\n", merged.isna().sum())
    print("row if we dropped NaNs:",merged.dropna().shape[0],"of",merged.shape[0])

    merged["batch"] = merged["batch"].fillna(1)
    merged["batch"] = merged["batch"].astype(int)

    print("missing after fill:",int(merged.isna().sum().sum()))
    return merged


#Explore
def explore(merged):
    hits = merged[(merged["condition"] == "treated") & (merged["count"]>200)]

    print("treated rows with count. > 200:",hits.shape[0])

    print("loc[0, 'gene']:", merged.loc[0, "gene"])
    print("iloc[0,0]:",merged.iloc[0,0])

    print("condition counts:\n", merged["condition"].value_counts())
    print("distinct condition:", merged["condition"].unique(),"| n =",merged["condition"].nunique())

    print("all sample IDs start with 'S':",merged["sample"].str.startswith("S").all())

    merged["sample_num"] = merged["sample"].str.replace("S","",regex=False).astype(int)
    merged["seq_date"] = pd.to_datetime(merged["seq_date"])
    merged["month"] = merged["seq_date"].dt.month

    merged["log_count"] = np.log1p(merged["count"])

    merged["level"] = merged["count"].apply(lambda c: "high" if c >= 200 else "low")

    print(merged[["gene","sample","count","log_count", "level", "month"]].head())

    return merged


#Group, Combine and Summarise
def summarise(merged):
    by_condition = (merged.groupby("condition").agg(mean_count=("count","mean"),n=("count","size")).reset_index())
    print("per condition:\n",by_condition)

    by_gene = (merged.groupby("gene").agg(mean_count = ("count","mean"),n=("count","size")).reset_index().sort_values("mean_count",ascending=False))
    print("per gene (top 3):\n", by_gene.head(3))

    ctrl = merged[merged["condition"]=="control"].groupby("gene")["count"].mean()
    trt  = merged[merged["condition"]=="treated"].groupby("gene")["count"].mean()

    combined = pd.concat([ctrl.rename("control_mean"),trt.rename("treated_mean")], axis=1).reset_index()

    combined["fold_change"] = (combined["treated_mean"] / combined["control_mean"]).round(2)

    print("control vs treated:\n", combined.sort_values("fold_change",ascending = False))

    stacked =np.vstack([ctrl.to_numpy(), trt.to_numpy()])
    print("vstacked shape (2 x genes):", stacked.shape)

    return by_gene


#Save results
def save_outputs(merged, by_gene):
    merged.to_csv(HERE/"tidy_expression.csv", index=False)
    by_gene.to_csv(HERE/"gene_summary.csv", index=False)

    print(f"Wrote {len(merged)} row to tidy_expression.csv " f"and {len(by_gene)} row to gene_summary.csv")


#Run the pipeline
if __name__ == "__main__":
    counts,meta = data_load_inspection()
    normalise_the_data(counts)
    long =reshape_long(counts)
    merged = merge_and_clean(long,meta)
    merged = explore(merged)
    by_gene = summarise(merged)
    save_outputs(merged,by_gene)
    print("\nDone.")
