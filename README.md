# `vdjrearranger`
#### Jamie Heather | 2026

`vdjrearranger` is a python tool to generate synthetic FASTQ sequence files in the format produced by the 10X V(D)J kit, for testing and validating different analyses involved in processing such data with `cellranger vdj` - particularly the testing of different germline reference sets.

It takes a TSV of paired chain clonotypes, containing at least three columns specifying the name of each clonotype, and the complete nucleotide sequence of the complete, spliced transcript for each of the receptor chains (e.g. as produced by [`thimble` in the `stitchr` package](https://jamieheather.github.io/stitchr/)), and *rearranges* it to generate paired end FASTQs, mocking up what cells expressing those receptors should look like.

Current it only covers `SC5P-R2` (Single Cell 5' R2-only) chemistry, and uses the [internal primer sequences available from the 10X website](https://kb.10xgenomics.com/s/article/360047454291-Which-genes-does-each-V-D-J-specific-primer-map-to) to trim the unused constant region sequences. Note this currently only covers human loci for which there are primers available.

`vdjrearranger` then goes through each clonotype, and then using a variety of parameters to assign:
* A number of different cell barcodes (drawn by default from the `cellranger`-provided `737K-august-2016.txt` file), to mimic a single cell in a GEM droplet.
* For each barcode, a number of different UMIS...
* ... with a range of numbers of reads...
* ... several of which can be either fixed or drawn from a different distribution (normal, Poisson, or power law).

Note that this tool doesn't aim to produce perfectly indistinguishable synthetic data, but the combination of parameters allows users to bias their data to replicate different types of data. Concessions towards realism include:

* Having FASTQ quality scores dip from ~Q40 to ~Q30 across reads
* Having a small proportion of noisy droplets, which have lower average quality
* The capacity to vary the transcriptional ratio of the two chains

Note that this has primarily only been tested with TRA/TRB data, and so some tweaking might be required for other receptors.

## Installation

Either locally:

```bash
git clone https://github.com/JamieHeather/vdjrearranger.git
pip install -e .
```

Or from PyPi:

```bash
pip install vdjrearranger
```

## Basic Usage

Two example clonotype files have been included in the `examples/` directory on the GitHub repo associated with this tool:
* `tcr_clonotypes.tsv` was produced by applying `thimble` on a bunch of TCR chains pulled out of random [VDJdb entries](https://github.com/antigenomics/vdjdb-db/)
* `bcr_clonotypes.tsv` was made by picking 100 paired light and heavy chain clonotype sequences from the `airr_rearrangement.tsv` output from a `cellranger` v10 run on the [10X-provided healthy human B cell tutorial dataset](https://www.10xgenomics.com/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0) 

### Coverage modes

There are several different modes that influence how 'reads' are subsampled from the length of the provided chains (post constant region trimming):

* `tiling`: Generates overlapping, bidirectional reads from both 5' and 3' directions to ensure comprehensive and uniform coverage across the full transcript. Will generate as many UMIs as needed to cover the whole transcript in both directions; if the requested value is set higher, it'll then repeat the process with random jitter around start points.
* `random`: Extracts fragments starting at uniformly distributed random positions across the entire transcript length.
* `3prime`: Concentrates fragment starts toward the 3' end of the transcript using an exponential decay distribution to mimic standard poly-A capture, aiming to mimic the 5'RACE bias.
* `normal`: Places fragment starts according to a normal distribution centered at a specified proportion of the transcript. Recommended to try and place the center around the approximate CDR3 location (so ~0.9 in a trimmed constant-region trimmed transcript).

### Minimal example

The following command can be run on the TCR example data (what a lot of the defaults are geared towards), to produce a dataset in which each clonotype occurs in a single cell, with fixed numbers of UMIs and reads, tiling across every transcript to get high uniform coverage.

```bash
vdjrearranger \
    --input-tsv examples/tcr_clonotypes.tsv \
    --outdir tcr_fixed 
```

Similarly this simplistic simulation can be applied to the test BCR data by providing the necessary column details:

```bash
vdjrearranger \
    --input-tsv examples/bcr_clonotypes.tsv \
    --outdir bcr_fixed \
    --chain1-col light_chain_seq \
    --chain2-col heavy_chain_seq \
    --clonotype-col clone_ID
```

Alternatively, users can simulate more complex scenarios potentially more representative of real data, with varying numbers of cells, UMIs, and reads per clonotype (with differing ratios of expression between the two chains). The following command demonstrates some of the possibilities, that might serve as a starting point for users to find the parameters that best suit their needs. 

```bash
vdjrearranger \
    --input-tsv examples/tcr_clonotypes.tsv \
    --outdir tcr_complex \
    --coverage-mode 3prime \
    --chain-ratio 0.3 \
    --cells-mode powerlaw \
    --cells-value 50 \
    --umis-mode normal \
    --umis-value 200 \
    --umis-variance 50 \
    --reads-mode poisson \
    --reads-value 100 \
    --reads-variance 20 \
    --noise-cells 200 \
    --noise-umis 2 
```

#### Parameter considerations

* The basic default parameters try to ensure that all cells and clonotypes will be detected in the output of `cellranger vdj` when run on the FASTQ files produced, not necessarily that these simulated data best reflect true wet lab results.
* It's not expected that every user will necessarily explore every parameter: they are merely provided in case it's useful in certain contexts, depending on their use case. 
* Generally the larger the values and the more complexity involve, the longer `vdjrearranger` will take to run, and the larger the output files will become.
* Curiously increasing the UMI and read count values doesn't always increase the likelihood of getting complete cell/clonotype detection.
* When I've tried to mimic ranges and distributions informed from real datasets, I often get shockingly low numbers of cell and clonotype recognition by `cellranger vdj` (at least when using small test datasets, which may not be representative of even typical smaller biological datasets). YMMV. 


## Outputs
Upon completion, the specified `--outdir` will contain:
1. `[sample_name]_S1_L001_R1/R2_001.fastq.gz`: The simulated sequence files, ready for input to `cellranger vdj`.
2. `truth_reads.tsv`: A complete mapping of every generated read ID back to its originating clonotype, chain, UMI, and barcode.
3. `rank_umi_plot.png`: A log-log 'knee plot' visualizing the UMI count separation between true cells and background noise.
4. `summary_stats.txt`: Summary file containing analogues of values output in the `cellranger` summary, including estimated cell count, reads per cell, and median UMIs per chain.
5. `run_command.txt`: The exact shell command executed to generate the run, for reproducibility purposes.

## Limitations

Currently `vdjrearranger` doesn't:
* Allow for cells expressing >2 recombined receptors
* Generate any PCR or sequencing errors
* Work in species other than human, or in gamma/delta chains 
  * Some combination of adding suitable primer sequences to the TSV or manually pre-trimming sequences should get you most the way there however.
