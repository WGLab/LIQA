# Tutorial and FAQs
## Tutorial
This section contains examples isoform analysis using LIQA. It is assumed that LIQA has been successfully installed. If not, please install it and its required packages first according to the [Installation](https://github.com/WGLab/LIQA/blob/master/doc/Install.md).

### Step 1. Downloading the example data

To prepare to run LIQA, please download the [LIQA_example](https://github.com/WGLab/LIQA/releases/tag/1.0.0) from release page. Then, extract example data by:
```
$ tar xvf LIQA_example.tar.gz
```
Next, please make sure current working directory is `LIQA_example` by:
```
$ pwd
/home/huy4/bin/LIQA_example
```
### Step 2. Converting the reference annotation file to isoform compatible matrix
The example reference annotation file is called sample.gtf and is in gtf format. This is the first few lines of the file:
```
chr21	HAVANA	gene	25639258	25717562	.	+	.	gene_id "ENSG00000154721.15"; gene_type "protein_coding"; gene_name "JAM2"; level 2; hgnc_id "HGNC:14686"; tag "overlapping_locus"; havana_gene "OTTHUMG00000078441.4";
chr21	HAVANA	transcript	25639258	25717562	.	+	.	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	exon	25639258	25639888	.	+	.	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 1; exon_id "ENSE00003843052.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	CDS	25639822	25639888	.	+	0	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 1; exon_id "ENSE00003843052.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	start_codon	25639822	25639824	.	+	0	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 1; exon_id "ENSE00003843052.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	exon	25683883	25683948	.	+	.	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 2; exon_id "ENSE00001017308.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	CDS	25683883	25683948	.	+	2	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 2; exon_id "ENSE00001017308.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
chr21	HAVANA	exon	25689866	25689973	.	+	.	gene_id "ENSG00000154721.15"; transcript_id "ENST00000480456.6"; gene_type "protein_coding"; gene_name "JAM2"; transcript_type "protein_coding"; transcript_name "JAM2-206"; exon_number 3; exon_id "ENSE00001017301.1"; level 2; protein_id "ENSP00000420419.1"; transcript_support_level "1"; hgnc_id "HGNC:14686"; tag "basic"; tag "MANE_Select"; tag "appris_principal_3"; tag "CCDS"; ccdsid "CCDS42911.1"; havana_gene "OTTHUMG00000078441.4"; havana_transcript "OTTHUMT00000171347.2";
```

Run the following command to convert this to the LIQA-compatible refgene file. Note: This input must be a GTF file - see the FAQs section below for conversion from GFF to GTF - or a UCSC file (example/example.ucsc)
```
liqa -task refgene -format gtf -ref sample.gtf -out sample.refgene
```
The isoform compatible matrix will be saved to `sample.refgene`. Below is a preview of how it should look.
```
DSCAM-AS1       chr21   +       40383082        40385358        ENST00000422749.5,ENST00000455354.1,ENST00000427451.5,ENST00000444046.5,
DSCAM-AS1       chr21   +       40383082        40383408        1,1,1,1,
DSCAM-AS1       chr21   +       40383409        40383455        0,1,1,1,
DSCAM-AS1       chr21   +       40383456        40383504        0,1,0,0,
DSCAM-AS1       chr21   +       40383505        40383584        1,1,0,0,
DSCAM-AS1       chr21   +       40383651        40384138        0,0,1,0,
DSCAM-AS1       chr21   +       40384578        40385358        1,1,1,1,
SOD1    chr21   +       31659621        31668931        ENST00000270142.10,ENST00000389995.4,ENST00000470944.1,ENST00000476106.5,
SOD1    chr21   +       31659621        31659664        1,0,0,0,
SOD1    chr21   +       31659665        31659691        1,1,0,0,
SOD1    chr21   +       31659692        31659707        1,1,0,1,
```
Column interpretation
- Column 1: Gene Name
- Column 2: Chromosome Number
- Column 3: Strand
- Column 4-5: Positions
- Column 6: Isoform expression matrix

### Step 3. Quantifying isoform expression using LIQA

There are 10 samples to analyse. The next step is to run the quantify isoform expression step for each sample. Run the following 10 commands. Note that with your own data, if an error is outputted, then you should sort the BAM file (see Usage + FAQs). 
```
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_1 -max_distance 20 -f_weight 1 -bam data/simu1.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_2 -max_distance 20 -f_weight 1 -bam data/simu2.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_3 -max_distance 20 -f_weight 1 -bam data/simu3.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_4 -max_distance 20 -f_weight 1 -bam data/simu4.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_5 -max_distance 20 -f_weight 1 -bam data/simu5.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_6 -max_distance 20 -f_weight 1 -bam data/simu6.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_7 -max_distance 20 -f_weight 1 -bam data/simu7.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_8 -max_distance 20 -f_weight 1 -bam data/simu8.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_9 -max_distance 20 -f_weight 1 -bam data/simu9.bam
liqa -task quantify -refgene sample.refgene -out estimation/isoform_expression_estimates_10 -max_distance 20 -f_weight 1 -bam data/simu10.bam
```
The progress of this command will be outputted as well.

Then, the results will be under the directory of `estimation`. Below is a preview of isoform_expression_estimates_1:
```
GeneName        IsoformName     ReadPerGene_corrected   RelativeAbundance       infor_ratio
DSCAM-AS1       ENST00000422749.5       9.33340258530289        0.14814924738576016     0.8571428571428571
DSCAM-AS1       ENST00000455354.1       45.499747588441075      0.7222182156895409      0.8571428571428571
DSCAM-AS1       ENST00000427451.5       5.833381993667815       0.0925933649788542      0.8571428571428571
DSCAM-AS1       ENST00000444046.5       2.3334678325882163      0.037039171945844704    0.8571428571428571
SOD1    ENST00000270142.10      192.0   0.9846153846153847      0.3333333333333333
SOD1    ENST00000389995.4       0.0     0.0     0.3333333333333333
SOD1    ENST00000470944.1       3.0     0.015384615384615385    0.3333333333333333
SOD1    ENST00000476106.5       0.0     0.0     0.3333333333333333
```

Column interpretation
- Column 1: Gene Name
- Column 2: Isoform Name from ENST Ensembl Transcript ID
- Column 3: Corrected read counts per gene from quantify step
- Column 4: Relative abundance of transcript in genome
- Column 5: Isoform expression estimate in form of information ratio from matrix

### Step 4. Detecting DAS genes using LIQA
In this step, the user provides two lists of isoform expression estimates files for condition 1 and 2 to use LIQA for differential splicing gene/isoform detection.

In this example, we assume that sample 1 to 5 belong to group 1 and sample 6 to 10 belong to group 2. The following lists summarized isoform expression estimates file for each group:
```bash
$ more list1
estimation/isoform_expression_estimates_1
estimation/isoform_expression_estimates_2
estimation/isoform_expression_estimates_3
estimation/isoform_expression_estimates_4
estimation/isoform_expression_estimates_5
$ more list2
estimation/isoform_expression_estimates_6
estimation/isoform_expression_estimates_7
estimation/isoform_expression_estimates_8
estimation/isoform_expression_estimates_9
estimation/isoform_expression_estimates_10
```
Then, please using following command to detect DAS genes:
```
liqa -task diff -condition_1 list1 -condition_2 list2 -out das_detection_results
```
Then, the results will be saved to `das_detection_results`.

The output file is a gene-based test results file with 3 columns:
		• Column 1: gene name
		• Column 2: test P-Value
		• Column 3: Hellinger distance which measures the splicing difference between two conditions for a gene in terms of isoform relative abundances - this column does not appear with example data below

```
ATP5PF  0.529504372538969
BACE2   0.999992837256582
DOP1B   0.000895020175699646
DSCAM-AS1       0.376916750256904
FTCD    0.0570208071493823
GABPA   0.557111938189986
IFNAR1  0.999632675437593
IL10RB  0.258904461845951
JAM2    0.989221122845281
OLIG1   0.460240770280204
PDXK    1
RIPPLY3 1
S100B   1
SOD1    6.7783673930899e-05
```

## FAQs
Ensure you have the latest version of LIQA.

1. Can I use my GFF file as an input for the refgene task?
LIQA only takes a GTF or UCSC as an input. To convert your GFF file, you can use the [gffread](https://github.com/gpertea/gffread) command as follows:
```
gffread -i gffinput -o gtfoutput
```
2. Where is my BAM file/why is it truncated?
In the instance that you have an empty/truncated BAM file, make sure to sort the BAM file, which you can do with SAMtools.

3. Does LIQA have a novel isoform detection tool?
Yes, more information can be found in [Usage](https://github.com/WGLab/LIQA/blob/master/doc/Usage.md)
   
For all other questions, feel free to check past issues to see if any of your questions were previously answered.
