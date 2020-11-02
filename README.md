TumOnc
======

Installation
------------

### Prerequisites

Install python3 and pip3. Run (you may probably skip this step. But do
it if something goes wrong!):

    pip3 install -U --user numpy Cython scipy pandas matplotlib setuptools

Required versions:

> -   setuptools \>= 18.0.1

### Install and prepare

untar `TumOnc-VER.tar.gz` in some directory. Run, ex:

    tar xvzf TumOnc-VER.tar.gz
    cd TumOnc-VER
    python3 setup.py build_ext --inplace

where VER is the version of the program.

Usage
-----

### Step 0 \-- prepare genome files

Should be present in the current directory:

    human_g1k_v37.fasta
    Agil_RS_names.bed

`Agil_RS_names.bed` is the bed file with the coordinates of the exons in the cancer panel analysed, and is provided for convenience here.

To get the human genome you can run

	curl ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz | gunzip > human_g1k_v37.fasta


### Step 1 \-- prepare background model

Run:

    TumOnc.sh --mutfile final_SNP_EX.txt --savepmodel pmodel.pcl --covariates gene.covariates_bcc_expr_hic.txt

This uses the mutations from final\_SNP\_EX.txt (better have a large set
of mutations! I guess do not use CP for this, most probably this will
lead to bad results).

The covariates are read from `gene.covariates_bcc_expr_hic.txt` Note,
that the presence of all three covariates with proper names in the first
line is essential. All genes without proper covariates are discarded by
the current version of the code. As an example, the
`gene.covariates_bcc_expr_hic.txt` contains the new covariates. By
default the contents of `gene.covariates.txt` is used (now it contains
the Lawrence covariates).

#### Using other probability models:

To get bagels add the following option to the command line above:

    --pmodel_bagels True

Note \-- as I realized, this is *not* the exact bagels from MutSig \--it
calculates over 3-context instead of a (small) set of categories used by
Lawrence!

### Step 2 \-- test the mutation database against the background model

Run:

    TumOnc.sh --mutfile final_SNP_EX.txt --loadpmodel pmodel.pcl --outfilebase pm3_final

This will produce two files: `pm3_final_mis.csv` and `pm3_final_sil.csv`
containing the results of analysis of the mutations in the mutfile using
the background model from pmodel.pcl for missense and silent mutations.

Note \-- the mutfile can be different from the one used in Step 1. For
example, you can use Exome mutations in step one, and cancer panel in
step 2.

Additional modifications:

    --pval4 1

This switch makes a better (\"4-dim reduction\") calculation of the
p-values. Should be probably always used \-- but it is very slow! And
slows down a lot with the number of patients in the database.

    --mimnut N

Use minimal number of mutations per nucleotide \-- all nucleotides with
smaller number of mutations are just ignored. Default is 3. Can probably
safely put this to 4 to increase speed (You will loose some 3 mutated
nucleotides, but who cares?)
