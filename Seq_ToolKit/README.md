# Seq ToolKit

A small tool that reads DNA sequences from a FASTA file and works out each one's length and GC content. Built with plain Python to practise core coding and Git.

## Run it

```
python3 analyze.py
```

This prints the results and saves `results.tsv` and `results.json`.

## Files

- `analyze.py` - the main program
- `sequences.fasta` - the input sequences
- `results.tsv`, `results.json` - the output (made when you run it)

## What it does

- Reads each sequence from the FASTA file
- Finds its length and GC content
- Skips empty sequences instead of crashing
- Keeps only the sequences I ask for
- Saves the results to a table and a JSON file

## What I practised

Lists, dicts, sets, tuples, loops, functions, a class, file reading/writing, error handling, and Git commits.
