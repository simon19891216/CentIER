# CentIER

Thank you for using CentIER

Before you get started with centIER, make sure you've installed LTR_retriever and gt on your own computer system and set up an environment variable pointing to them, giving that variable the name 'LTR_retriever' and 'gt', respectively. 

You can head over to https://github.com/oushujun/LTR_retriever and https://github.com/genometools/genometools/releases/tag/v1.6.5 for downloading the LTR_retriever and gt, respectively.

The Python version we used during the testing phase is Python 3.8.

The system we used during the testing phase is Ubuntu 20.04 (64-bit)

# download and install

You can get CentIER, at the ["Releases" ](https://github.com/simon19891216/CentIER/releases/tag/CentIERv2.0) module, and the ColCEN.fasta file is an example data (T2T genome of Arabidopsis).

After downloaded the CentIERv3.0:

tar -zxvf CentIERv3.0.tar.gz

cd CentIERv3.0

cp the_pathway_of_genome the_directory_containing_CentIERv3.0.py (Copy the genome to the same directory level as CentIERv3.0.py)

If you want to add the gff file (genome annotation file，the addition of GFF files can improve prediction accuracy.)：

cp the_pathway_of_gff the_directory_containing_CentIERv3.0.py (Copy the genome to the same directory level as CentIERv3.0.py)

Then, you can run centIER using the following command(the example data is the T2T genome of Arabidopsis thaliana. We strongly recommend that users run the program with example data to check if it is working properly before using their own data.):

# Options and Usage

The "-g" option is mandatory and is used for inputting the genome sequence. We recommend that users remove contigs from the genome sequence file and retain only the chromosome sequences.

python CentIERv3.0.py -g ColCEN1.fasta

The following options are optional.

The "-gff" option is used to input a GFF (General Feature Format) annotation file.

python CentIERv3.0.py -g ColCEN1.fasta -gff example.gff

The "-kmer" option is used to specify the k-mer size (default value is 21) in order to identify regions of sequence repetition.

python CentIERv3.0.py -g ColCEN1.fasta -gff example.gff -kmer 20

The "-c" option (default value is 15) is used to fine-tune the distance between two regions.

python CentIERv3.0.py -g ColCEN1.fasta -gff example.gff -kmer 20 -c 20

The "-step_len" option (default value is 10000) is used to motify the distance between two regions.

python CentIERv3.0.py -g ColCEN1.fasta -gff example.gff -kmer 20 -c 20 -step_len 20000

The "-mul_cents" option is used to identify centromeric regions of metapolycentric chromosomes.

python CentIERv3.0.py -g ColCEN1.fasta -gff example.gff -kmer 20 -c 20 -step_len 20000 -mul_cents

The distribution of Hi-C signals can also be utilized to predict centromeres. Users need to preprocess the Hi-C data using a separate script developed by our group in order to obtain the distribution of Hi-C signals.

The script can be downloaded from https://github.com/ruoyu1123/centromFind. 

The script will generate four files, which are then imported into CentIER using parameters “-matrix1”, “-matrix2”, “-bed1”, and “-bed2” respectively.

The remaining two parameters specify the size of the gap region ("-MINGAP", default value:2 (means 0.2 Mb)) and the threshold for the signal ("-SIGNAL_THRESHOLD", default value:0.7), respectively.

All results can be found in the folder of **CentIER_final_results**

Moreover, during the process of identifying the centromere range, some intermediate files will be generated and retained.

If you encounter any issues while running the program, you can contact me directly via email (xudongzhuanyong@163.com).
