# Kracken-Bracken-plot
 Creates stacked bar plots from Kracken and Bracken output. The plot is a vectorized pdf file that can be opened in Adobe Illustrator or other for additional modification.



Builds bar plots from Bracken results with optional adjustment for
unclassified reads from Kraken results.

Bracken modifies Kraken results for a more acurate abundance estimate
which is great but Bracken does not include the unclassified reads. I
think it is important for a thorough analysis to report the fraction of
the metagenome that is not classified compared to what is classified.
Environmental samples are frequently >50% unclassified.

Therefore, this script has the option of accepting the Kraken reports
that the Bracken reports were generated from. It will read the total
number of unclassified reads from the Kraken report and adjust the
"Fraction of Total Reads" reported by Bracken because this "Fraction of
Total Reads" is only the fractional total of reads assigned to the
"Level ID" specified (the unclassified reads are not counted). With this
optional adjustment, our plots will show the fraction of total reads
(including unclassified).

input file: List of Bracken report files, SampleName

-i, input_file_specificiations:
    This is a .txt file with one entry per line in the order you want.
    The output plot will follow the order in this file.

    Example:
                Path/to/BrackenFilename, SampleA
                Path/to/BrackenFilename, SampleB
                Path/to/BrackenFilenmae, SampleC

optional input:

-k, kraken_file_specifications:
    This is a .txt file with one entry per line. The Sample names must
    match the input file sample neames and the Kraken files should 
    correspond to those use for the Bracken files. The order of this
    file does not matter but I find it is easiest if the order is the
    same as the input file.

    Example:
                Path/to/KrakenFilename, SampleA
                Path/to/KrakenFilename, SampleB
                Path/to/KrakenFilename, SampleC

-l, custom_legend_file:
    This is a .txt file with one entry per line of classification, color.
    The classification must match those used in the input file.
    The colors can be whatever you want. I prefer to use hex colors.
    Example: 
                Alteromonadales, #2ca25f
                Sphingobacteriales, #8856a7
                Rickettsiales, #43a2ca

output: stacked bar plots in vector pdf format. Publication ready.

Ref Links:
https://ccb.jhu.edu/software/kraken2
https://github.com/DerrickWood/kraken2
https://ccb.jhu.edu/software/bracken
https://github.com/jenniferlu717/Bracken