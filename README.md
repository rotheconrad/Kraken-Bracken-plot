# Plot Kraken, Kraken 2 and Bracken results
 Creates stacked bar plots from Kraken 2 and Bracken output. The plot is a vectorized pdf file that can be opened in Adobe Illustrator or other for additional modification.

![Kraken-Bracken plot](https://github.com/rotheconrad/Kraken-Bracken-plot/blob/main/Kraken-Bracken-plot_example_genus.png)

Builds bar plots from Bracken results with optional adjustment for
unclassified reads from Kraken results.

## Motivation

Bracken modifies Kraken 2 results for a more acurate abundance estimate
which is great but Bracken does not include the unclassified reads. For
more a thorough analysis, this scripts incorporates the fraction of
the metagenome that is not classified compared to what is classified.
This is important since environmental samples are frequently >50% unclassified.

Therefore, this script has the option of accepting the Kraken reports
that the Bracken reports were generated from. It will read the total
number of unclassified reads from the Kraken report and adjust the
"Fraction of Total Reads" reported by Bracken because this "Fraction of
Total Reads" is only the fractional total of reads assigned to the
"Level ID" specified (the unclassified reads are not counted). With this
optional adjustment, our plots will show the fraction of total reads
(including unclassified).

## Usage

```bash
python Kraken-Bracken-plot.py -h
```

Input files are all two column tab ("\t") separated files (tsv).

input file: List of Bracken report files, SampleName

-i, input_file_specificiations:
    This is a .txt file with one entry per line in the order you want.
    The output plot will follow the order in this file.

    Example:
                Path/to/BrackenFilename\tSampleA
                Path/to/BrackenFilename\tSampleB
                Path/to/BrackenFilename\tSampleC

optional input:

-k, kraken_file_specifications:
    This is a .txt file with one entry per line. The Sample names must
    match the input file sample neames and the Kraken files should 
    correspond to those use for the Bracken files. The order of this
    file does not matter but I find it is easiest if the order is the
    same as the input file.

    Example:
                Path/to/KrakenFilename\tSampleA
                Path/to/KrakenFilename\tSampleB
                Path/to/KrakenFilename\tSampleC

-l, custom_legend_file:
    This is a .txt file with one entry per line of classification, color.
    The classification must match those used in the input file.
    The colors can be whatever you want. I prefer to use hex colors.

    Example: 
                Alteromonadales\t#2ca25f
                Sphingobacteriales\t#8856a7
                Rickettsiales\t#43a2ca

output: stacked bar plots in vector pdf format. Publication ready.

## Running Kraken2 and Bracken

Example batch processing scripts for an HPC running a slurm scheduler are
included for reference under `sbatch/`. Kraken and Bracken both rely on a 
Kraken database which can be obtained by the tool directly via `kraken2-build`
as either a pre-built database or one created *ad hoc* by the user.

## Requirements and References

Requires:

* Python 3.6+
* Matplotlib
* Pandas
* Kraken2
* Bracken

Ref Links:

* https://ccb.jhu.edu/software/kraken2
* https://github.com/DerrickWood/kraken2
* https://ccb.jhu.edu/software/bracken
* https://github.com/jenniferlu717/Bracken

# Customized figure example:

<img width="930" alt="image" src="https://github.com/rotheconrad/Kraken-Bracken-plot/assets/36962040/67d1abae-0ac8-48f3-825e-7a7aeecf9b85">

