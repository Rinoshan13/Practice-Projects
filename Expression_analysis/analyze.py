import pandas as pd
import numpy as np
from pathlib import Path

HERE = Path(__file__).parent # Expression_analysis Folder

# Load and Read the data file into DataFrame and inspect
def data_load_inspection():
    counts = pd.read_csv(HERE/"counts.csv")
    meta = pd.read_csv(HERE/"metadata.csv")

    # print(counts.head())
    # print("counts shape:",counts.shape)
    # print("columns:",list(counts.columns))
    # print("dtype:\n", counts.dtypes)

    # counts.info()

    # print(counts.describe())

    # print("metadata:\n", meta.tail())

    return counts, meta

counts,meta= data_load_inspection()


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







normalise_the_data(counts)








