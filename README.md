# lfq_proceessing
Quantitative analysis of maxquant results

## Install required dependencies

To make sure all the required pachages are installed run the _install.R_ file in the parent folder of the cloned repo.

Assuming that git repository has benn cloned in the home directory and the application folder location is __/home/ubuntu/maxquant-quant-analysis-public__ then run:

```r
app_folder = "/home/ubuntu/maxquant-quant-analysis-public/"
source(file = file.path(app_folder, "install.R"))
```

___


## Set up enviroinment

### Input folder setup

The first step is to assign variables that point to the location of the maxquant output files and to the cloned git repository.

```r
mq_folder = "/home/ubuntu/data/mq-test/"
app_folder = "/home/peppe/maxquant-quant-analysis-public/"
```

The _mq_folder_ has to contain these files:

- _evidence.txt_
- _modificationSpecificPeptides.txt_
- _msms.txt_
- _parameters.txt_
- _peptides.txt_
- _proteinGroups.txt_
- _summary.txt_


### Experimental design info

Experimental dessign needs to be written in a tab separated file formatted like the example below:

| file_name | replicate | experiment | mqExperiment |
|------:|------:|------:|------:|
JD_06232014_sample1_A|1|sample1|1_A
JD_06232014_sample1_B|2|sample1|1_B
JD_06232014_sample1_C|3|sample1|1_C
JD_06232014_sample2_A|1|sample2|2_A
JD_06232014_sample2_B|2|sample2|2_B
JD_06232014_sample2_C|3|sample2|2_C
JD_06232014_sample3_A|1|sample3|3_A
JD_06232014_sample3_B|2|sample3|3_B
JD_06232014_sample3_C|3|sample3|3_C
JD_06232014_sample4_A|1|sample4|4_A
JD_06232014_sample4_B|2|sample4|4_B
JD_06232014_sample4_C|3|sample4|4_C

Where the column __file_name__ contains the names of the LC-MS/MS runs in the data, the column __experiment__ contains the experimental conditions, the column __replicate__ is a unique number for each of the biological or technical replicate in one experimental contition, the column __mqExperiment__ contains the values in the "Experiment" column in the Maxquant "Raw data" tab when setting up the Maxquant analysis (figure below).


![Maxquant experimental set up. Each "File" has to be associated to a unique value in the "Experiment" column](docs/MQ_Example_setup.png)
Maxquant experimental set up. Each "File" has to be associated to a unique value in the "Experiment" column


--- 
## Quantitative analysis

To load the main function into the enviroinment run:

```r
source(file = file.path(app_folder, "app/mq_transformer.R"))

```



Set up a variable with the name of the experiment design file, the function assumes the files is located in the same folder of the Maxquant output

```r
experimentalDesign_filename <- "experimentDesign.txt"

```



Running the quantitative analysis:

```r
run_analysis <- mq_transfomer(
    app_folder = app_folder, 
    mq_folder = mq_folder, 
    expdes_filename = experimentalDesign_filename
)
```


---

#### Imputation parameters 

Missing values imputation is done as described in the paper [here](www.massdynamics.com) 

Imputation parameters can be changed by altering the default values in the __mq_transfomer__ function:

```r
run_analysis <- mq_transfomer(
    app_folder = app_folder,
    mq_folder = mq_folder, 
    expdes_filename = experimentalDesign_filename,
    mputeStDev = 0.3,
    imputePosition = 1.8
)
```
