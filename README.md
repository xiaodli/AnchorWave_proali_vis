# AnchorWave `proali` subcommand's anchors result visualization

## Dotplot visualization

The length of the query and reference corresponding to each facet is equal to the actual chromosome length.

## Example

### maize(query) fai file

The following is part maize fai file and in order to determine whether the length of `chr10` has been modified to the actual chromosome length(152M->552M).

```text
chr1    308452471    6    80    81
chr2    243675191    312308139    80    81
chr3    238017767    559029276    80    81
chr4    250330460    800022272    80    81
chr5    226353449    1053481869    80    81
chr6    181357234    1282664743    80    81
chr7    185808916    1466288949    80    81
chr8    182411202    1654420483    80    81
chr9    163004744    1839111832    80    81
chr10    552435371    2004154143    80    81
```

### sorghum(reference) fai file

The following is part sorghum fai file.

```text
1    80884392    71    60    61
2    77742459    82232608    60    61
3    74386277    161270846    60    61
4    68658214    236896966    60    61
5    71854669    306699555    60    61
6    61277060    379751873    60    61
7    65505356    442050289    60    61
8    62686529    508647472    60    61
9    59416394    572378848    60    61
10    61233695    632785589    60    61
```

### Running R code

```R
source("dotplot.R")
plot_total("anchors", "query.fai", "ref.fai", "query", "ref", "anchors.png")
```

### Dotplot

<p align="center">
<img src="./anchors.png" alt= anchors.png width="800px" background-color="#ffffff" />
</p>

## line style visualization

You can vis multiple anchors file by raw_line_proali module

```bash
git clone git@github.com:xiaodli/AnchorWave_proali_vis.git
cd raw_line_proali
python main.py -h
python main.py line_proali -i 1.Kronos_Svevo.anno.3.anchormove,2.Svevo_XM001097.anno.3.anchormove,3.XM001097_NU00021.anno.3.anchormove,4.NU00021_IG77365.anno.3.anchormove,5.IG77365_IG99236.anno.3.anchormove,6.IG99236_PI294478.anno.3.anchormove,7.PI294478_NU01905.anno.3.anchormove,8.NU01905_NU01954.anno.3.anchormove,9.NU01954_Zavitan.anno.3.anchormove -o ten.line.png -l Kronos.length.txt,Svevo.length.txt,XM001097.length.txt,NU00021.length.txt,IG77365.length.txt,IG99236.length.txt,PI294478.length.txt,NU01905.length.txt,NU01954.length.txt,Zavitan.length.txt -n Kronos,Svevo,XM001097,NU00021,IG77365,IG99236,PI294478,NU01905,NU01954,Zavitan -rm "0" -cf 7 -sf 7 -it -sc "red,blue" -cs "four_colors" -al -gs "compact"
```

1. The following is length file (chr and length column is necessary, tab sep and no matter how many columns the file has, only the first two are used).

```text
chr  length  total_gene
1A  611319266  3772
1B  749393777  4158
2A  793380695  4950
2B  832257163  5398
3A  762016299  4530
3B  864951163  4947
4A  766591834  4214
4B  703924636  3578
5A  721915144  4785
5B  734861517  4869
6A  633522200  3588
6B  742827932  3996
7A  761831191  4622
7B  762926078  4011
```

2. prefix is species name.

3. remove_chromosome_prefix is chromosome prefix(comma separated).

### line plot

<p align="center">
<img src="./ten.line.png" alt= ./ten.line.png width="800px" background-color="#ffffff" />
</p>
