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
In order to run "HelitronScanner" which is a software for the detection of helitron transposons, please run following command:
```
helitronscanner -help
helitronscanner scanHead -g demo.fasta -bs 0 -o scanHead.txt
helitronscanner scanTail -g demo.fasta -bs 0 -o scanTail.txt
helitronscanner pairends -hs scanHead.txt -ts scanTail.txt -o result.txt
```

