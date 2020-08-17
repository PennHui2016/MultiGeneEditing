## MultiGeneEditing
 
In this repository, a tool bulit with java swing and the python source codes are provided.

# Running the java swing tool
The following requirements should be installed at first:
1. java: jdk1.8.0_131 or later is recommended, adding the install path to the environment variable;
2. python 3: python 3.7 or later is recommended, also adding the install path to the environment variable;
3. python packages including scipy, numpy, pandas should be installed.

To use the tool, please following the guidence below:
1. download this repository, there will be two folders named "java" and "python source codes";
2. **download the required data from https://drive.google.com/file/d/1nyhlgcLlsbDYa4KbQ0S53LIvAAZXf3Cf/view?usp=sharing, the required data is named "data.zip"**;
3. extract the data.zip to the folder **"java"**;

the complete structure should be:
![example_folder](https://github.com/PennHui2016/images/blob/master/eg_folder1.png)
./MultiGeneEditing-master
├─**java
│  ├─build
│  │  └─classes
│  │      ├─images
│  │      └─multigeneediting
│  ├─code
│  │  └─__pycache__
│  ├─dist
│  ├─nbproject
│  │  └─private
│  ├─src
│  │   ├─images
│  │   └─multigeneediting
│  └─data
│       ├─multi_genes
│       ├─predict_genes_results
│       ├─single_editing_modify
│       ├─single_gene_MPCP1
│       └─two_gene_editing**
└─python source codes

4. double click the file "MultiGeneEditing.jar" under folder java, the tool will be open. Users can select one or more genes (if more than one genes need to be selected, just press the Ctrl and click those genes) and press "execute" button to submit the request. Finally, the results will be given.

# Running the python source codes
If users prefer to run the python source codes, then
python3 should be installed and the required packages such as scipy, numpy, pandas should be installed.

Then:
1. download this repository, there will be two folders named "java" and "python source codes";
2. **download the required data from https://drive.google.com/file/d/1nyhlgcLlsbDYa4KbQ0S53LIvAAZXf3Cf/view?usp=sharing, the required data is named "data.zip"**;
3. extract the data.zip to the folder **"python source codes"**;

the complete structure should be:
![example_folder](https://github.com/PennHui2016/images/blob/master/eg_folder2.png)
./MultiGeneEditing-master
├─java
│  ├─build
│  │  └─classes
│  │      ├─images
│  │      └─multigeneediting
│  ├─code
│  │  └─__pycache__
│  ├─dist
│  ├─nbproject
│  │  └─private
│  └─src
│      ├─images
│      └─multigeneediting
└─**python source codes
                    ├─data
                    │    ├─multi_genes
                    │    ├─predict_genes_results
                    │    ├─single_editing_modify
                    │    ├─single_gene_MPCP1
                    │    └─two_gene_editing
                    └─results**
Then, one can run the codes with the following command line:
  
   *python display_result.py input1 input2 input3*
   
   where input1 is the list of genes to be edited, e.g., "ENSG00000002834,ENSG00000005073,ENSG00000005339"
         input2 is a random number to make the name of result file.
         input3 is the parent folder of the folder "data", e.g., "./python source codes"
   
   an example command is 
   *python display_result.py ENSG00000002834,ENSG00000005073,ENSG00000005339 123 ./python source codes*
   
 after running the example command, the result file in the csv format can be found in the folder "./python source codes/results/multi_gene_spacer_123.csv"
 
 Any problem please contact Dr. Hui Peng, email: cdph2009@13.com
