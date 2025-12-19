# AminoHeartz
Enjoy the sound of proteins, based on science.
> **Protein sonification based on infrared spectroscopy and biophysics.**
> *Listen to the raw physics of life.*

[VISIT THE WIKI FOR IN DEEP INFORMATION](https://github.com/ampinzonv/aminoheartz/wiki/The-theory)

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Experimental-orange)]()

## What is AminoHeartz?

**AminoHeartz** is a CLI (Command Line Interface) tool written in Python that translates protein sequences (FASTA files) into audio.

Unlike other sonification algorithms that assign arbitrary musical notes to amino acids (e.g., A = C, G = D), **AminoHeartz uses molecular physics principles**. The sound you hear is an audible transposition of the real vibrational frequencies of chemical bonds (Amide I, Amide II, Disulfide Bridges, etc.) that make up the protein.

## Key Features

* **Real Physics:** Frequencies derived from Infrared Spectroscopy (IR) data.
* **Scanner Mode (Sliding Window):** Generates evolving textures that traverse the protein from the N-terminal to the C-terminal.
* **Additive Synthesis:** Builds timbre by summing pure sine waves for each chemical bond.
* **Verbose Mode:** Visualizes in the console which "notes" and frequencies are dominating your protein.

## Installation

1.  Clone this repository:
    ```bash
    git clone [https://github.com/your-username/aminoheartz.git](https://github.com/your-username/aminoheartz.git)
    cd aminoheartz
    ```

2.  Install the necessary dependencies:
    ```bash
    pip install numpy scipy
    ```

## Usage

The script takes a `.fasta` file as input and generates a `.wav` file.

### 1. Evolving Sonification (Recommended)
Creates a texture that changes as it "reads" the sequence.
```bash
python bio_sonificator.py insulin.fasta -o insulin_evo.wav -s 4 -w 10
