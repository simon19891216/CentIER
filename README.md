# CenTIER

## Contents

1. [Introduction](#introduction)
2. [Install](#install)
3. [Usage](#usage)
4. [Output](#output)
5. [Cite](#cite)

## INTRODUCTION

CentIER is a bioinformatics tool for identifying and annotating centromere regions of the genome that can be accurately identified.Instead of identifying centromere by identifying tandem repeats, CentIER was designed with sequence features including tandem repeats, retrotransposons, k-mer frequency distributions in mind, and allows users to further define centromere regions by entering a genome structure annotation file and a Hi-C signal file. While recognizing monocentric Chromosomes, we also designed a program to recognize the centromere regions of Metapolycentric Chromosomes.  

By predicting the centromere of Arabidopsis, rice, corn, soybean, mulberry and other species, we found that CentIER maintained an accuracy of more than 90 percent. In maize, where centromere is difficult to predict, CentlER accurately predicted the centromere region of 9 out of 10 chromosomes.  

## INSTALL

The tool has the following software or package requirements.  

1. **genometools**: As its name, **gt** is a genome tools and the repository is https://github.com/genometools/genometools.git. You can install gt in `apt-get install genometools` when you are root user. Otherwise you only use `make -j4` to Compile and install it. The gt requires `cario` and `pango`, you can use conda to install them. If you don't need it, use the `cairo=no` make parameter. In that case, run `make cairo=no` cleanup before running `make cairo=no` to start with a clean build environment. We recommend using [gt_V1.6.5](https://github.com/genometools/genometools/releases/tag/v1.6.5).  
2. **LTR_retriever**: This is used to annotation genome and centremer regions. We recommend you build a **conda environment** to install it and use this enviroment to running **CentIER**. 
3. Python package：**numpy**, **pandas**, **pyfastx**, **scipy**. You can install them by `pip install`. Usually, You only need install **pyfastx** because **ltr_retriever** env have include other packages. You do not need to specify version number of packages. If some update case error, you can specify `ltr_retriever=3.0.0`, `pyfastx=2.1.0`, `pandas=2.2.2`, `numpy=2.0.1`, `scipy=1.14.0`.

The following is an example of the installation procedure  

``` sh
# Create centier enviroment and install ltr_retriever
conda create -n centier -c bioconda -c conda-forge ltr_retriever=3.0.0
conda activate centier
pip install pyfastx==2.1.0

# install gt
# if you have not download gt, please use this link to get it.
wget https://github.com/genometools/genometools/archive/refs/tags/v1.6.5.tar.gz

tar -zxf genometools_v1.6.5.tar.gz
cd genometools_v1.6.5
make -j4
export PATH=$PATH:your_PATH/genometools_v1.6.5/bin

# install CentIER
# if you have a dictoinary name software to install software. 
git clone  https://github.com/simon19891216/CentIER.git
# Add the CentIER to PATH. you can add this code to .bashrc.
export PATH=$PATH:your_PATH/CentIER
```

## USAGE

```shell
$ ./CentIERv3.0.py -h
usage: CentIERv3.0.py [-h] [-o OUTPUT] [--gff GFF] [-k KMER_SIZE]
                      [-c CENTER_TOLERANCE] [--step_len STEP_LEN]
                      [--mul_cents] [--matrix1 MATRIX1]
                      [--matrix2 MATRIX2] [--bed1 BED1] [--bed2 BED2]
                      [--MINGAP MINGAP]
                      [--SIGNAL_THRESHOLD SIGNAL_THRESHOLD]
                      genome

positional arguments:
  genome                Genome fasta file need to annotation

options:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output path
  --gff GFF             optional, annotation file (gff or gtf format)
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        the size of kmer
  -c CENTER_TOLERANCE, --center_tolerance CENTER_TOLERANCE
                        the fine-tuning size between two regions
  --step_len STEP_LEN   the size between two regions
  --mul_cents           This parameter is used to retain all potential
                        centromeric regions, especially when the
                        chromosomes of the species under study may not
                        have a single, defined centromere; it should be
                        included in such cases.
  --matrix1 MATRIX1     Path to matrix of 100000 file
  --matrix2 MATRIX2     Path to matrix of 200000 file
  --bed1 BED1           Path to bed file of matrix1
  --bed2 BED2           Path to bed file of matrix2
  --MINGAP MINGAP       Minimum gap value n*100000 (default: 2)
  --SIGNAL_THRESHOLD SIGNAL_THRESHOLD
                        Signal threshold value (default: 0.7)
```
**CentIER** have one positional argument **genome.fasta**. If you only have a T2T assembly (T2T level or Chromosome level genomes), you can running this tools as follow:  

``` shell
CentIER.py genome.fasta
```

If you have gff3 annotation file, you can running program as follow:

```shell
CentIER.py --gff genome.gff genome.fasta
```

If you have Hi-C data, you can use Hi-C matrix to identify regions of centremers more accuracy. The Hi-C matrix can be obtaion from **[Hi-Cpro](https://github.com/nservant/HiC-Pro.git)**. You can use the the Hi-C matrix as follow:  

```shell
CentIERv3.0.py --matrix1 sample.100000.matrix --matrix2 sample.20000.matrix --bed1 sample.100000_abs.bed --bed2 sample.20000_abs.bed genome.fasta
```
For more detailed use, or if you just want to locate the centromere region through Hi-C, you can refer to [CentromFind](https://github.com/ruoyu1123/centromFind.git) repository.

Depend on our test result, the **BUG** that ids will not be match if your Chr_ID of genome fasta does not follow the form of the following case.

```shell
# fasta sample. The sequence id must be "ChrN". if use Chr_1, the program will error.
>Chr1
AAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACC
>Chr2
CTGGTTGAAAATCATTGTGTATATAATGATAAT
>Chr3
GCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGT
>Chr4
ACCCTAAACCCTAAACCCTAAACCCTAA
>Chr5
CCTAAACCCTAAACCCTAAACCCTAAACCCTA
>Chr6
TCTTATTAATTGTTTGGACTGTTTATGTTTGGACATT
......
```

## OUTPUT

The **output_path** of **CentIER** can be specifed by `-o` and if you don't use `-o`, the default output path is `./CentIER_final_results`. The output file list as follow:  
   
```shell
$ ls
ColCEN.fasta_all_centromere_seq.txt
ColCEN.fasta_all_centromere_seq.txt_l1.txt
ColCEN.fasta_centromere_range.txt
ColCEN.fasta_draw_cen.svg
ColCEN.fastaLTR-hmmresult.txt
ColCEN.fasta_LTR_positionb.txt
ColCEN.fasta_ltr_position.txt
ColCEN.fasta_ltr_position.txt_seqt.txt
ColCEN.fasta_ltr_position.txt_seq.txt
ColCEN.fasta_LTR_statistics.txt
ColCEN.fasta_monomer_in_centromere.txt
ColCEN.fasta_monomer_seq.txt
```
The program will produce many tmp file with start with basename of your genome file in current path (./your_genome.fa.*). After running, you can remove them manually.
>**Wait for more information**  

## Cite  
if you use CentIER or [CentromFind](https://github.com/ruoyu1123/centromFind.git),
please cite this atricle:  
[CentIER: accurate centromere
identification for plant genome](https://doi.org/10.1016/j.xplc.2024.101046)  
Xu, D. et al. CentIER: accurate centromere identification for plant genome. Plant Communications 101046 (2024) doi:10.1016/j.xplc.2024.101046.

## Contact  

If you have any problem while running the program, you can contact us through the ways as follow:  

* Email: xudongzhuanyong@163.com  
* QQ Group: 1001056790  

You can submit issues in this repository too. Thanks for your using our program.  
**Author**: **Xu dong**, **Yang jinbao** 