# INTRODUCTION

Thank you for using CentIER. 


Before you get started with centIER, make sure you've installed LTR_retriever on your own and set up an environment variable pointing to it, giving that variable the name 'LTR_retriever'. 

Just head over to https://github.com/oushujun/LTR_retriever to download and get your hands on the program.

Furthermore, your system's Python version should be Python 3.0 or higher, and you should have the pyfastx and biopython Python packages installed.

After downloaded the CentIERv2.0:

tar -zxvf CentIERv2.0.tar.gz

cd CentIERv2.0

cp the_pathway_of_genome ./

if you want to add the gff file (genome annotation file，the addition of GFF files can improve prediction accuracy.)：

cp the_pathway_of_gff ./

Then, you can run centIER using the following command(the example data is the T2T genome of Arabidopsis thaliana. We strongly recommend that users run the program with example data to check if it is working properly before using their own data.):

python CentIERv2.0.py -g the_pathway_of_genome -gff the_pathway_of_annotation_file (-gff is optional)

all results can be found in the folder of CentIER_final_results

Moreover, several intermediate files will be generated, and you can choose to delete these temporary files by adding the -c option.

If you encounter any issues while running the program, you can contact me directly via email (xudongzhuanyong@163.com).