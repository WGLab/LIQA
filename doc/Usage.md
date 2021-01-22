# Instruction about how to use LIQA

The inputs of LIQA are aligned long-read RNA-seq data in BAM format and a reference isoform annotation file (Ensembl/Refseq). User needs to specify `-task` to perform in each step:
```
liqa -task <task>:

        refgene:   preprocess reference file
        
        quantify:   quantify isoform expression
        
        diff:   detect differential splicing gene/isoform
        
```

## Step 1: Transform isoforms to compatible matrix based on reference annotation file
LIQA requires a reference annotation file `example.refFile` in following format:
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
Reference file in this format can be downloaded at [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?command=start) by selecting "all fields from selected table" in output format.

We preprocess `example.refFile` by using `liqa -task refgene`. An example is given below.
```
liqa -task refgene -ref example.refFile -out example.refgene
```

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

## Step 3: Detect differential splicing gene/isoform between conditions
In this step, user needs to give two lists of `isoform exression estimates` files for condition 1 and 2 to LIQA to detect differential splicing gene/isoform.
```
liqa    -task diff
        -condition_1 <list_of_isoform_expression_estimation_file_for_condition1>
        -condition_2 <list_of_isoform_expression_estimation_file_for_condition2>
        -out <test_results_file>

```
Here are example of input lists for condition 1 and 2:

List_Condition_1
```
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample1
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample2
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample3
```
List_Condition_2
```
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample1
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample2
/home/huy4/tmp/liqa_das_tmp/isoform_estimates_condition1_sample3
```

Output:
```
	"geneDASResults" file: Gene-based test results file with 3 columns:
		• Column 1: gene name
		• Column 2: test P-Value
		• Column 3: Hellinger distance which measures the splicing difference between two conditions for a gene in terms of isoform relative abundances
		
```
