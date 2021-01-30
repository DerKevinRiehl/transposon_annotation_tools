# transposon_annotation_tools
A set of bioconda packages for transposon annotations. During my masterthesis I downloaded lots of these tools and I want to make it easier for the research community to install and run these softwares as bioconda packages from command line. 


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
sinescan -s 123 -g demo.fasta -o result.txt -d result/
```


## **Using TirVish** [CondaPackage](https://anaconda.org/bioconda/genometools), [Publication](https://ieeexplore.ieee.org/abstract/document/6529082), [Code](http://genometools.org/tools/gt_tirvish.htmln)
In order to run "TirVish" which is a software for the detection of TIR transposons, please run following command:
```
gt suffixerator -db demo.fasta -indexname demo.index -tis -suf -lcp -des -ssp -sds -dna -mirrored
gt tirvish -index demo.index > result.txt
```


## **Using LtrHarvest** [CondaPackage](https://anaconda.org/bioconda/genometools), [Publication](https://link.springer.com/article/10.1186/1471-2105-9-18), [Code](https://www.zbh.uni-hamburg.de/forschung/gi/software/ltrharvest.html)
In order to run "LtrHarvest" which is a software for the detection of LTR transposons, please run following command:
```
gt suffixerator -db demo.fasta -indexname demo.index -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index demo.index > result.txt
```


## **Using RepeatModeler** [CondaPackage](https://anaconda.org/bioconda/repeatmodeler), [Code](http://www.repeatmasker.org/RepeatModeler/)
In order to run "RepeatModeler" which is a software for the detection of all classes of transposons and repeats, please run following command:
```
BuildDatabase -name demo_index -engine ncbi demo.fasta
RepeatModeler -engine ncbi -pa 10 -database demo_index
```



