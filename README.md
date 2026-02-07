**Digital Comms Lab 2**
# Huffman Coding for Bernoulli Sequences and Erdős–Rényi Graphs

## Overview
This repository contains MATLAB scripts for implementing Huffman coding on:

1. **Bernoulli sequences** – simple binary sequences with specified probabilities.
2. **Erdős–Rényi (ER) random graphs** – adjacency matrices of graphs with varying edge probability.

The scripts perform:
- Probability estimation
- Entropy calculation
- Huffman encoding and decoding (lossless)
- Compression ratio calculation
- Graphical representation of ER graphs and symbol probabilities

## Files

| File | Description |
|------|-------------|
| `Ex1L2.m` | Implements Huffman coding on a Bernoulli sequence |
| `Ex2L2.m` | Generates ER graphs, computes adjacency, compresses via Huffman |
| `Figures/` | Contains all output figures for network plots, probability charts, and compression analysis |
| `README.md` | This file |

## Requirements
- MATLAB R2022a or newer
- Communications Toolbox (for `huffmandict`, `huffmanenco`, `huffmandeco`)
- Optional: MATLAB plotting functions (`graph`, `bar`, `plot`)  

## How to Run
1. Open MATLAB and set the current folder to this repository.
2. Run `Ex1L2.m` to see Bernoulli sequence compression.
3. Run `Ex2L2.m` to generate ER graph compression results and figures.
