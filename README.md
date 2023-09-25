# CentIER
It is easy to use CentIER. For example, you can set the commond as python CentIER.py ColCEN.fasta 21 30 30. Run the commond and all analysis will be completed automatically.

Explanation of the commond line:

The first parameter refers to the pathway of the input genomic sequence file;

The second parameter refers to the size of the k-mer, such as 21 in the example;

The third parameter refers to the maximum distance between two intervals with abnormal k-mer signals that can be merged, such as 30 (the first) in the example.

The fourth parameter refers to the tolerance value for the left and right boundaries of the merged interval, such as 30 (the second) in the example.

The third and fourth parameters will be changed, according to the size of the input genome.

Generally, the parameters in the example are appropriate (these parameters have been tested on the genomes of Arabidopsis thaliana, rice and maize.) 

