# transposon_annotation_tools
A set of bioconda packages for transposon annotation and transposon feature annotation in nucleotide sequences. *transposon_annotation_tools* is part of [TransposonUltimate](https://github.com/DerKevinRiehl/TransposonUltimate).

## Installation
You can simply create a conda environment and install the tools you want as outlined in the following. We recommend an environment based on Python=2.7 but depending on the combination of packages you can try to install the packages as well in different manors. As complex dependencies cause long waiting times for environment resolving using conda, we recommend the use of mamba. 
* **Note1:** *Sinescan installation needs python=2.7. As conda is not taking the most recent version of the package, you need to specify the latest version number (1.1.2).* 
* **Note2:** *For some users the bioconda channel is reported to cause issues with genometools-genometools, therefore you might consider to download it from other channels, e.g. conda-forge: "conda install -y -c bioconda -c conda-forge genometools-genometools".* 
* **Note3:** *For some users sinescan was causing trouble with environment resolving of conda, therefore mamba worked better.*

**Installation using mamba (recommended)**
```
conda create -y --name transposon_annotation_tools_env python=2.7
conda activate transposon_annotation_tools_env
conda install -y mamba
mamba install -y -c bioconda genometools-genometools # for some users: mamba install -y -c bioconda -c conda-forge genometools-genometools
mamba install -y -c derkevinriehl transposon_annotation_tools_proteinncbicdd1000
conda install -y -c derkevinriehl transposon_annotation_tools_transposonpsicli # this one seems to not work with mamba
mamba install -y -c derkevinriehl transposon_annotation_tools_mitetracker
mamba install -y -c derkevinriehl transposon_annotation_tools_helitronscanner
mamba install -y -c derkevinriehl transposon_annotation_tools_mitefinderii
mamba install -y -c derkevinriehl transposon_annotation_tools_mustv2
mamba install -y -c derkevinriehl transposon_annotation_tools_sinefinder

# how to install sinescan
mamba install -y python=2.7 # if not done before of mentioned while creating the environment
mamba install -y -c derkevinriehl transposon_annotation_tools_sinescan=1.1.2

conda deactivate
```

**Installation using conda**
```
conda create -y --name transposon_annotation_tools_env python=2.7
conda activate transposon_annotation_tools_env
conda install -y -c bioconda genometools-genometools # for some users: conda install -y -c bioconda -c conda-forge genometools-genometools
conda install -y -c derkevinriehl transposon_annotation_tools_proteinncbicdd1000
conda install -y -c derkevinriehl transposon_annotation_tools_transposonpsicli
conda install -y -c derkevinriehl transposon_annotation_tools_mitetracker
conda install -y -c derkevinriehl transposon_annotation_tools_helitronscanner
conda install -y -c derkevinriehl transposon_annotation_tools_mitefinderii
conda install -y -c derkevinriehl transposon_annotation_tools_mustv2
conda install -y -c derkevinriehl transposon_annotation_tools_sinefinder

# how to install sinescan
conda install -y python=2.7 # if not done before of mentioned while creating the environment
conda install -y -c derkevinriehl transposon_annotation_tools_sinescan=1.1.2

conda deactivate
```

## What you will find in this repository
In the following list you will find name, publication, URL to conda package, URL to software and a short tutorial on how to run the package for each of the transposon annotation tools included into this package. Moreover, there is a *Software_Manual.pdf* in this git (containing the help-lines of all softwares mentioned below), the *demo.fasta* used in the examples below (not all softwares hit findings on this sample), as well as a folder *Manuals* with additional PDF files of the software authors.

## **Using MUSTv2** [CondaPackage](https://anaconda.org/DerKevinRiehl/transposon_annotation_tools_mustv2), [Publication](https://doi.org/10.1515/jib-2017-0029), [Code](http://www.healthinformaticslab.org/supp/resources.php)
In order to run "MUSTv2" which is a software for the detection of MITE transposons, please run following command:
```
mustv2 -help
mkdir temp
mustv2 demo.fasta result.txt temp
```


## **Using HelitronScanner** [CondaPackage](https://anaconda.org/derkevinriehl/transposon_annotation_tools_helitronscanner), [Publication](https://doi.org/10.1073/pnas.1410068111), [Code](https://sourceforge.net/projects/helitronscanner/files/)
In order to run "HelitronScanner" which is a software for the detection of HELITRON transposons, please run following command:
```
helitronscanner -help
helitronscanner scanHead -g demo.fasta -bs 0 -o scanHead.txt
helitronscanner scanTail -g demo.fasta -bs 0 -o scanTail.txt
helitronscanner pairends -hs scanHead.txt -ts scanTail.txt -o result.txt
```


## **Using SineFinder** [CondaPackage](https://anaconda.org/derkevinriehl/transposon_annotation_tools_sinefinder), [Publication](https://doi.org/10.1105/tpc.111.088682)
In order to run "SineFinder" which is a software for the detection of SINE transposons, please run following command:
```
sine_finder -help
sine_finder -V demo.fasta
```


## **Using MiteTracker** [CondaPackage](https://anaconda.org/derkevinriehl/transposon_annotation_tools_mitetracker), [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2376-y), [Code](https://github.com/INTABiotechMJ/MITE-Tracker)
In order to run "MiteTracker" which is a software for the detection of MITE transposons, please run following command:
```
mitetracker -help
mkdir results
mitetracker -g demo.fasta -j jobName -w 3
```


## **Using MiteFinderII** [CondaPackage](https://anaconda.org/derkevinriehl/transposon_annotation_tools_mitefinderii), [Publication](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-018-0418-y), [Code](https://github.com/jhu99/miteFinder)
In order to run "MiteFinderII" which is a software for the detection of MITE transposons, please run following command:
```
miteFinderII -help
miteFinderII -input demo.fasta -output result.txt
```


## **Using SineScan** [CondaPackage](https://anaconda.org/derkevinriehl/transposon_annotation_tools_sinescan), [Publication](https://doi.org/10.1093/bioinformatics/btw718), [Code](https://github.com/maohlzj/SINE_Scan)
In order to run "SineScan" which is a software for the detection of SINE transposons, please run following command:
```
sinescan -help
mkdir result
mkdir output
mkdir final
sinescan -s 123 -g demo.fasta -o output -d result -z final
```


## **Using TirVish** [CondaPackage](https://anaconda.org/bioconda/genometools), [Publication](https://ieeexplore.ieee.org/abstract/document/6529082), [Code](http://genometools.org/tools/gt_tirvish.htmln)
In order to run "TirVish" which is a software for the detection of TIR transposons, please run following command:
```
gt tirvish -help
gt suffixerator -db demo.fasta -indexname demo.index -tis -suf -lcp -des -ssp -sds -dna -mirrored
gt tirvish -index demo.index > result.txt
```


## **Using LtrHarvest** [CondaPackage](https://anaconda.org/bioconda/genometools), [Publication](https://link.springer.com/article/10.1186/1471-2105-9-18), [Code](https://www.zbh.uni-hamburg.de/forschung/gi/software/ltrharvest.html)
In order to run "LtrHarvest" which is a software for the detection of LTR transposons, please run following command:
```
gt ltrharvest -help
gt suffixerator -db demo.fasta -indexname demo.index -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index demo.index > result.txt
```


## **Using RepeatModeler** [CondaPackage](https://anaconda.org/bioconda/repeatmodeler), [Code](http://www.repeatmasker.org/RepeatModeler/)
In order to run "RepeatModeler" which is a software for the detection of all classes of transposons and repeats, please run following command:
```
RepeatModeler -help
BuildDatabase -name demo_index -engine ncbi demo.fasta
RepeatModeler -engine ncbi -pa 10 -database demo_index
```


## **Using RepeatMasker** [CondaPackage](https://anaconda.org/bioconda/repeatmasker), [Code](http://www.repeatmasker.org/)
In order to run "RepeatMasker" which is a software for the detection of all classes of transposons and repeats, please run following command:
```
RepeatMasker -help
RepeatMasker -pa 10 demo.fasta
```


## **Using TransposonPSI** [CondaPackage](https://anaconda.org/DerKevinRiehl/transposon_annotation_tools_transposonpsicli), [Code](http://transposonpsi.sourceforge.net/)
In order to run "TransposonPSI" which is a software for the detection of proteins characteristic for transposons, please run following command:
```
mkdir temp
mkdir result
transposonPSI -fastaFile demo.fasta -resultFolder result -tempFolder temp -mode nuc 
# modes: 'nuc' and 'prot'
```


## **Using TransposonProteinNCBICDD1000** [CondaPackage](https://anaconda.org/DerKevinRiehl/transposon_annotation_tools_proteinncbicdd1000), [Code](https://github.com/DerKevinRiehl/transposon_annotation_tools/)
In order to run "TransposonProteinNCBICDD1000" which is a software for the detection of proteins characteristic for transposons, please run following command:
```
mkdir result
proteinNCBICDD1000 -fastaFile demo.fasta -resultFolder result 
```
This software uses 1000 selected characteristic conserved domain models of proteins from the [NCBI CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml) that are found to occur frequently in transposons. The tool [RPSTBLASTN](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is used to annotate these proteins within a given fasta file.


## Background
During my masterthesis, I downloaded lots of these tools and I want to make it easier for the research community to install and run these softwares as bioconda packages from command line. Therefore, I packaged these softwares into conda packages. Please note: I am not the author of the softwares. My contribution lies in creating conda packages to allow a broader bioinformatic audience to use these tools. To find more specific information on the softwares, please refer to the URLs mentioned above. The only package created by myself within this GIT repository is *TransposonProteinNCBICDD1000*.


## Citations
Please cite our paper if you find TransposonUltimate useful:

Riehl, Kevin and Riccio, Cristian and Miska, Eric and Hemberg, Martin. TransposonUltimate: software for transposon classification, annotation and detection. bioRxiv doi: https://doi.org/10.1101/2021.04.30.442214 (submitted to GenomeBiology, available on BioRxiv https://www.biorxiv.org/content/10.1101/2021.04.30.442214v1 )

```
@article{riehl2021transposonultimate,
  title={TransposonUltimate: software for transposon classification, annotation and detection},
  author={Riehl, Kevin and Riccio, Cristian and Miska, Eric and Hemberg, Martin},
  journal={bioRxiv preprint https://www.biorxiv.org/content/10.1101/2021.04.30.442214v1},
  year={2021}
}
```

## Acknowledgements
We would like to thank Sarah Buddle, Simone Procaccia, Fu Xiang Quah and Alexandra Dallaire for their assistance with testing and debugging the software.
