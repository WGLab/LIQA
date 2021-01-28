# Instruction about how to use LIQA

The inputs of LIQA are aligned long-read RNA-seq data in BAM format and a reference isoform annotation file (Ensembl/Refseq). User needs to specify `-task` to perform in each step:
```
liqa -task <task>:

        refgene:   preprocess reference file
        
        quantify:   quantify isoform expression
        
        diff:   detect differential splicing gene/isoform
        
```

## Step 1: Transform isoforms to compatible matrix based on reference annotation file
LIQA accepts two formats of reference annotation file. User can download example reference file and data under the [example directory](https://github.com/WGLab/LIQA/tree/master/example):
### GTF format
For example:
```
chr1    ncbiRefSeq      exon    24828834        24828953        .       +       .       gene_id "RCAN3"; transcript_id "NM_001251979.2"; exon_number "1"; exon_id "NM_001251979.2.1"; gene_name "RCAN3";
chr1    ncbiRefSeq      exon    24840804        24841057        .       +       .       gene_id "RCAN3"; transcript_id "NM_001251979.2"; exon_number "2"; exon_id "NM_001251979.2.2"; gene_name "RCAN3";
```
### UCSC all fields
For example:
```
749     NM_001397       chr1    -       21543739        21616982        21546447        21616907        19      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,21616856,  21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616649,21616982,     0       ECE1cmpl     cmpl    0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,0,
93      NM_001113348    chr1    -       21543739        21672034        21546447        21671871        19      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,21671868,  21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616649,21672034,     0       ECE1cmpl     cmpl    0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,0,
749     NM_001113349    chr1    -       21543739        21616766        21546447        21616691        18      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,21616562,   21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21605825,21616766,      0       ECE1    cmpl    cmpl0,0,1,2,2,0,0,0,0,2,0,0,0,0,1,1,0,0,
749     NM_001113347    chr1    -       21543739        21606183        21546447        21605927        17      21543739,21548239,21551742,21553651,21554423,21560050,21562342,21563238,21564626,21571481,21573713,21582439,21584017,21585185
,21586763,21599191,21605683,    21546624,21548335,21551933,21553719,21554534,21560154,21562420,21563337,21564737,21571596,21573856,21582631,21584083,21585332,21586885,21599404,21606183,       0       ECE1    cmpl    cmpl    0,0,1,2,2,0,0
,0,0,2,0,0,0,0,1,1,0,
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

## Step 2: Quantify isoform expression
In this step, user needs to give  `refgene_File`, `bam_file` to LIQA to estimate isoform expression using long-read RNA-seq data:
```
liqa -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file> -max_distance <max distance> -f_weight <weight of F function>
```
where
```
<bam>_file: A bam file.
<refgene_file>: A reference file obtained from step 1.
<max distance>: The maximum length of an alignment error at exon boundary. Recommend: 20.
<weight of F function>: The weight for bias correction in isoform usage estimation. Recommend: 1
```

For example:
```
liqa -task quantify -refgene example.refgene -bam example.bam -out isoform_expression_estimates -max_distance 20 -f_weight 1
```
## Step 3: Detect differential splicing gene/isoform between conditions
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
		• Column 3: Hellinger distance which measures the splicing difference between two conditions for a gene in terms of isoform relative abundances
		
```
