# Expression Analysis

A small project that takes a gene-count table, cleans it, does some maths on it, and saves a tidy summary. Built with pandas and numpy to practise data work and Git.

## Run it

```
python3 analyze.py
```

This prints the results and saves `tidy_expression.csv` and `gene_summary.csv`.

## Files

- `analyze.py` — the main program
- `counts.csv` — gene counts (genes x samples)
- `metadata.csv` — sample info (condition, batch, date)
- `tidy_expression.csv`, `gene_summary.csv` — the output (made when you run it)

## What it does

- Loads the counts and the sample info
- Normalises the counts with numpy (so samples can be compared)
- Reshapes the table from wide to long
- Joins in the sample info and fills a missing value
- Groups by condition and gene to get averages
- Saves the tidy table and a summary

## What I practised

pandas (read_csv, filtering, groupby, merge, pivot/melt, missing data) and numpy (arrays, broadcasting, axis maths, masks), plus Git commits — one per step.
