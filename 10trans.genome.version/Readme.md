# 不同版本编码基因名字转换
```
两个版本的gff文件的基因名字不能相同,
也就是假如版本1有一万个基因名字，版本2也有相同的一万个基因名字，如果你要使用这个流程进行不同版本基因名字的对应。
你需要将版本2的基因名字加上一个"_"或者其他的符号，与版本1区分
也就是gff文件中，gene行的ID=后的基因名字后加一个下划线或者其他符号，还有mRNA行的  Parent=  后的基因名字加同样的符号。
```
## 安装quota_Anchor并下载一些可有可无的脚本

```bash
conda create -n quota_Anchor bioconda::quota_anchor -y
git clone git@github.com:baoxingsong/quota_Anchor.git
cp -r quota_Anchor/scripts .
```

## 下载ncbi和ensembl的玉米基因组

```bash
mkdir raw_data
cd raw_data
mwget -n 16 https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_902167145.1/download\?include_annotation_type\=GENOME_FASTA\&include_annotation_type\=GENOME_GFF
mwget -n 16 https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-61/plants/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.61.chr.gff3.gz
mwget -n 16 https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-61/plants/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
```

把基因组和gff文件放在raw_data文件夹下, 并且修改名字如下所示

```bash
tree raw_data
```
```
├── Zea.mays.fa
├── Zea.mays.gff3
├── zm.ncbi.fa
└── zm.ncbi.gff3
```

## 提取最长转录本  1）不支持GTF文件，2）将基因组和注释文件名字分别改为 gff3 fa后缀

```bash
conda activate quota_Anchor
python ./scripts/longest_pipeline.py -i raw_data -o output_dir --overwrite
```

报错信息如下，这是由于NCBI有些行中strand位置没有写 +  -链信息，NCBI玉米注释版本太早了（通常是非染色体片段才有这种问号）

```text
[2025/07/03 18:16:01 INFO] Error parsing strand (?) from GFF line:
[2025/07/03 18:16:01 INFO] NC_007982.1	RefSeq	mRNA	50490	267232	.	?	.	ID=rna-ZeamMp186;Parent=gene-ZeamMp186;Dbxref=GeneID:4055939;exception=trans-splicing%2C RNA editing;gbkey=mRNA;gene=nad1;locus_tag=ZeamMp186
```

所以最好提前把非染色体的注释给删掉，下边我提前把非染色体的注释行给删掉

```bash
awk '{print $1}' raw_data/zm.ncbi.gff3|grep -v "#"|uniq|less -SN
```

确定第一个非染色体的名字叫什么，这里我只确定了报错的非染色体片段的名字确定叫做 NC_007982.1, 在Gff文件中确定了NC_007982.1位于990391行，删除GFF文件中该行及之后的行

```bash
sed -i '990391,$d' raw_data/zm.ncbi.gff3
```

## 再次提取最长转录本

```bash
python ./scripts/longest_pipeline.py -i raw_data -o output_dir --overwrite
```

## 统计染色体长度和编码基因数目信息

```bash
find ./raw_data/*fai |awk '{printf "%s,", $1}'
find ./raw_data/*gff3 |awk '{printf "%s,", $1}'
find ./raw_data/*gff3 |awk '{printf "%s,", $1}'|sed s/gff3/length\.txt/g  
```

上边三个命令可以帮助你生成参数，不用自己手动输入, -s 参数含义使用quota_Anchor get_chr_length -h 查看

```bash
quota_Anchor get_chr_length -f ./raw_data/Zea.mays.fa.fai,./raw_data/zm.ncbi.fa.fai -g ./raw_data/Zea.mays.gff3,./raw_data/zm.ncbi.gff3 -s 0-9,CHR,chr,Chr:0-9,CHR,chr,Chr,NC -o ./raw_data/Zea.mays.length.txt,./raw_data/zm.ncbi.length.txt --overwrite
```

## blast 或者diamond比对, 并且结合了gff文件的信息生成了一个table文件, ncbi玉米作为query，ensembl玉米作为reference

```bash
quota_Anchor pre_col -a diamond -qs output_dir/01longest/zm.ncbi.longest.pep -rs output_dir/01longest/Zea.mays.longest.pep -db output_dir/01longest/Zea.mays.longest.pep.diamond     -mts 20 -e 1e-5 -b zea.mays_ncbi.blast -qg raw_data/zm.ncbi.gff3 -rg raw_data/Zea.mays.gff3 -o zea.mays_ncbi.table -bs 0 -al 0 -ql raw_data/zm.ncbi.length.txt -rl raw_data/Zea.mays.length.txt --overwrite
```

## 共线性分析，下边这个算法和JCVI中的quota_align不同，JCVI中的quota_align结果中存在一些基因对不符合你设置的R Q值（需要你写脚本自己去处理）， 下边这个是严格符合R，Q的，下边的结果严格符合1:1的

```bash
quota_Anchor col -i zea.mays_ncbi.table -o zea.mays_ncbi.table.collinearity -s 0 -W 0 -E -0.005 -D 25 -r 1 -q 1 --overwrite -m -1
```

## 之前的步骤已经完成了不同版本基因组的id对应，下边是可视化

```bash
quota_Anchor dotplot -i zea.mays_ncbi.table.collinearity  -o zea.mays_ncbi.table.collinearity.png -r raw_data/Zea.mays.length.txt -q raw_data/zm.ncbi.length.txt -t order -r_label "Ensembl" -q_label "Ncbi" -w 1500 -e 1200 --overwrite                      
```

## Note
```
你也可以使用当前目录下的shell脚本代替quota_Anchor `pre_col` `col`两个步骤,从而使用你自定义参数生成的blast文件而不是quota_Anchor所使用的blast参数（不论是blastp还是diamond都存在极个别序列相同但是比对不上的情况，因此如果你想或者更准确的结果，建议增加max_target_seq参数的大小，比如设置为50. 现在很多人可能会将该参数设置为10或者20左右.对于diamond你也可以使用--sensitive模式）.
```
