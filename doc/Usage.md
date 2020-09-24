# Instruction about how to use LIQA

The inputs of LIQA are aligned long-read RNA-seq data in BAM format and a reference isoform annotation file (Ensembl/Refseq). User needs to specify `-task` to perform in each step:
```
LIQA.py -task <task>:

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

We preprocess `example.refFile` by using `python LIQA.py -task refgene`. An example is given below.
```
python LIQA.py -task refgene -ref example.refFile -out example.refgene
```

## Step 2: Quantify isoform expression
In this step, user needs to give  `refgene_File`, `bam_file` to LIQA to estimate isoform expression using long-read RNA-seq data:
```
python LIQA.py -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file> -max_distance <max distance>
```
where
```
<bam>_file: A bam file.
<refgene_file>: A reference file obtained from step 1.
<max distance>: The maximum length of an alignment error at exon boundary. Recommend: 10.
```

## Step 3: Detect differential splicing gene/isoform between conditions
In this step, user needs to give  `refgene_File`, `isoform exression estimates` to LIQA to detect differential splicing gene/isoform.
```
python LIQA.py -task diff -ref <reference_file> -est <isoformRelativeAbundances_estimations>
```
Outputs of `python LIQA.py -task das` are `das_*.sh` script files located at `./tmp/das_script`. User needs to run all of them to obtain differential alternative splicing even at exon group level across cell conditions.

```
python LIQA.py -task sum -gpinfo example.gpinfo
```
