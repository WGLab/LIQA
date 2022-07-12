# Instruction about how to use LIQA

The inputs of LIQA are aligned long-read RNA-seq data in BAM format and a reference isoform annotation file (Ensembl/Refseq). User needs to specify `-task` to perform in each step:
```
liqa -task <task>:

        refgene:   preprocess reference file
        
        quantify:   quantify isoform expression
        
        diff:   detect differential splicing gene/isoform
        
```

## Step 1: Transforming isoforms to compatible matrix based on reference annotation file
LIQA accepts two formats of reference annotation file. User can download the example reference file and data from the [example folder](https://github.com/WGLab/LIQA/tree/master/example):
### GTF format
For example `example.gtf`:
```
chr21   HAVANA  gene    25639258        25717562        .       +       .       gene_id "ENSG00000154721.15"; gene_type "protein_coding"; gene_name "JAM2"; level 2; hgnc_id "HGNC:14686"; tag "overlapping_locus"; havana_gene "OTTHUMG00000
078441.4";
chr21   HAVANA  transcript      25639258        25717562        .       +       .       gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; tran
script_name "JAM2-206"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078
441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21   HAVANA  exon    25639258        25639888        .       +       .       gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_n
ame "JAM2-206"; exon_number 1; exon_id "ENSE00003843052.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42
911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21   HAVANA  CDS     25639822        25639888        .       +       0       gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_n
ame "JAM2-206"; exon_number 1; exon_id "ENSE00003843052.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42
911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
```
### UCSC all fields
For example `example.ucsc`:
```
947     ENST00000397748.5       chr21   -       47556175        47575481        47556367        47575437        15      47556175,47556828,47557152,47558421,47558793,47565330,47565731,47566179,47570032,47570301,47571471,47571805,47572820,47574062,47575383,      47556408,47556967,47557248,47558560,47558837,47565492,47565861,47566241,47570164,47570439,47571651,47571894,47572949,47574246,47575481, 0       FTCD    cmpl    cmpl    1,0,0,2,0,0,2,0,0,0,0,1,1,0,0,
947     ENST00000397743.1       chr21   -       47556692        47575481        47557159        47575437        13      47556692,47557152,47558421,47565330,47565731,47566179,47570032,47570301,47571471,47571805,47572820,47574062,47575383,47556987,47557248,47558560,47565492,47565861,47566241,47570164,47570439,47571651,47571894,47572949,47574246,47575481,   0       FTCD    cmpl    cmpl    -1,1,0,0,2,0,0,0,0,1,1,0,0,
947     ENST00000397746.7       chr21   -       47556692        47575481        47556900        47575437        14      47556692,47557152,47558421,47558793,47565330,47565731,47566179,47570032,47570301,47571471,47571805,47572820,47574062,47575383,       47556987,47557248,47558560,47558837,47565492,47565861,47566241,47570164,47570439,47571651,47571894,47572949,47574246,47575481,  0       FTCD    cmpl    cmpl    0,0,2,0,0,2,0,0,0,0,1,1,0,0,
```
This reference file can be downloaded at [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?command=start) by selecting "all fields from selected table" in output format.


Then, we preprocess reference file by using `liqa -task refgene`. An example is given below.

- For gtf format
```
liqa -task refgene -ref example.gtf -format gtf -out example.refgene
```
- For ucsc all fields:
```
liqa -task refgene -ref example.ucsc -format ucsc -out example.refgene
```
The output isoform compatible matrix file `example.refgene` should be in following format:
```
RCAN3   chr1    +       24828840        24863510        NM_001251983,NM_001251978,NM_013441,NM_001251984,NM_001251979,NM_001251981,NM_001251977,NM_001251985,NM_001251982,NM_001251980,
RCAN3   chr1    +       24828840        24828953        0,0,0,1,1,0,0,0,0,0,
RCAN3   chr1    +       24829386        24829607        1,0,1,0,0,1,1,0,0,0,
RCAN3   chr1    +       24829608        24829640        1,0,1,0,0,0,0,0,0,0,
RCAN3   chr1    +       24834083        24834218        0,1,0,0,0,0,0,0,0,0,
RCAN3   chr1    +       24840803        24841057        0,1,1,0,1,1,1,1,1,1,
RCAN3   chr1    +       24857707        24857881        1,1,1,1,1,0,1,0,1,1,
RCAN3   chr1    +       24859572        24859601        1,1,1,1,1,1,1,0,0,0,
RCAN3   chr1    +       24859602        24859744        1,1,1,1,1,1,1,0,0,1,
RCAN3   chr1    +       24861582        24863510        1,1,1,1,1,1,1,1,1,1,
```
**Note**: user needs to specify correct reference file format (gtf or ucsc) in this step.

## Step 2: Quantifying isoform expression
In this step, user is first suggested to perform reads filtering using samtools:
```
samtools view bam_file -F 2308 -q 50 -O BAM -o bam_filtered
```

Then, user needs to sort and index the `bam_file`, and give `refgene_File`, `bam_file` to LIQA to estimate isoform expression using long-read RNA-seq data:
```
liqa -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file> -max_distance <max distance> -f_weight <weight of F function>
```
where
```
<bam>_file: A bam file.
<refgene_file>: A reference file obtained from step 1.
<max distance>: The maximum length of an alignment error at exon boundary. Recommend: 10.
<weight of F function>: The weight for bias correction in isoform usage estimation (between 0 and 10). Recommend: 1
```

For example:
```
liqa -task quantify -refgene example.refgene -bam example.bam -out isoform_expression_estimates -max_distance 20 -f_weight 1
```
## Step 3: Detecting differential splicing gene/isoform between conditions
In this step, user needs to give two lists of isoform expression estimates files for condition 1 and 2 to use LIQA for differential splicing gene/isoform detection.
```
liqa    -task diff
        -condition_1 <list_of_isoform_expression_estimation_file_for_condition1>
        -condition_2 <list_of_isoform_expression_estimation_file_for_condition2>
        -out <test_results_file>

```
Here are example format for input lists:

List_Condition_1
```
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample1
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample2
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample3
```
List_Condition_2
```
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition2_sample1
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition2_sample2
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition2_sample3
```

**Note**: user needs to provide the absolute path for each LIQA isoform estimates file

User can download the DAS detection example data from the [release page](https://github.com/WGLab/LIQA/releases/tag/1.0.0) and follow the command in [Demo and Examples](https://github.com/WGLab/LIQA/blob/master/doc/Examples.md).

Output:
```
	"geneDASResults" file: Gene-based test results file with 3 columns:
		• Column 1: gene name
		• Column 2: test P-Value
		
```
