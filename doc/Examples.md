# Demo and Examples
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
```
liqa -task refgene -format gtf -ref sample.gtf -out sample.refgene
```
Then, isoform compatible matrix will be saved to `sample.refgene`.

### Step 3. Quantifying isoform expression using LIQA
There are 10 samples to analyse. Please use following command to quantify isoform expression for each sample:
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
Then, the results will be under the directory of `estimation`.

### Step 4. Detecting DAS genes using LIQA
It is assumed that sample 1 to 5 belong to group 1 and sample 6 to 10 belong to group 2. The following lists summarized isoform expression estimates file for each group:
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

