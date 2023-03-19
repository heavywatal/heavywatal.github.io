+++
title = 'biomaRt'
subtitle = "プログラム的にデータ取得"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -33
+++

<https://www.bioconductor.org/packages/release/bioc/html/biomaRt.html>

<https://qiita.com/yuifu/items/a757629506c1cd98156b>

BioMartからプログラマチックにデータを取得するための [Bioconductor]({{< relref "bioconductor.md" >}}) 拡張。
ウェブアプリ [MartView](http://www.biomart.org/biomart/martview) のGUIでひととおり慣れておくと良い。

## クラス

`Mart`
:   `@biomart` --- `ensembl`

    `@host` --- `http://www.biomart.org:80/biomart/martservice`

    `@dataset` --- `scerevisiae_gene_ensembl`

    `@filters`
    :   `$name` --- `chromosome_name`, `biotype`, ...\
        `$options` --- `[ncRNA,protein_coding,pseudogene,rRNA,snoRNA,snRNA,tRNA]`, ...\
        `$operation` --- `=`, `>=`, `only,excluded`, ...\
        `$description`\
        `$fullDescription`\
        ...

    `@attributes`
    :   `$name` --- `ensembl_gene_id`, `entrezgene`, ...\
        `$page` --- `feature_page`, `structure`, `homologs`, `snp`, `sequences`\
        `$description`\
        `$fullDescription`\
        ...

## 関数

`listMarts(mart, host="www.biomart.org", ..., archive=FALSE, ...)`
:   利用可能なマート列挙。 e.g. `ensembl`, `fungi_mart_21`, `unimart`

`listDatasets(mart, verbose=FALSE)`
:   マートで利用可能なデータセット列挙。 e.g. `scerevisiae_gene_ensembl`

`attributePages(mart)`
:   マートで利用可能なデータのカテゴリ分け(ページ)を列挙。
    `unique(mart@attributes$page)` と同じ。
    データを取得するときはこれらの間を跨がないようにアトリビュートを選ぶ。

`listAttributes(mart, page, what=c("name", "desciption"))`
:   マートで利用可能なアトリビュート列挙。
    `mart@attributes` へのアクセサ。

`listFilters(mart, what=c("name", "description"))`
:   マートで利用可能なフィルター列挙。
    `mart@filters` へのアクセサ。

`useMart(martname, dataset, host="www.biomart.org", ..., archive=FALSE, ...)`
:   利用するマート(とデータセット)を指定して `Mart` オブジェクトを作る。
    データセットを決めずにも作れる。 e.g. `mart = useMart("ensembl", "scerevisiae_gene_ensembl")`

`useDataset(dataset, mart, verbose=FALSE)`
:   データセットを決めた `Mart` オブジェクトを作る。
    でもこれって `useMart()` にもできるので不要...？ e.g. `mart = useDataset("scerevisiae_gene_ensembl", mart)`

`getBM(attributes, filters="", values="", mart, curl=NULL, checkFilters=TRUE, verbose=FALSE, uniqueRows=TRUE, bmHeader=FALSE)`
:   このパッケージのメイン関数。
    フィルターは名前付きリストで渡す。
    `filter` に渡すと `AND` 結合:

        getBM(c("ensembl_gene_id", "sgd_gene"),
              filters=list(chromosome_name="IV", biotype="tRNA"),
              mart=ensembl)
        # 第四染色体上のtRNA遺伝子

    `values` に渡すと `OR` 結合:

        getBM(c("ensembl_gene_id", "sgd_gene"),
              values=list(chromosome_name="IV", biotype="tRNA"),
              mart=ensembl)
        # 第四染色体上の遺伝子とtRNA遺伝子

    このへんとか `useDataset()` らへんとか、
    インターフェイスがあまり洗練されてない印象だなぁ...。

`getSequence(chromosome, start, end, id, type, seqType, upstream, downstream, mart)`
:   Sequencesページのデータをダウンロードすることに特化した `getBM()` ラッパー。 結果は2列(配列とID)の data.frame。\
    **seqType**: `gene_exon`, `transcript_exon`, `transcript_exon_intron`, `gene_exon_intron`, `cdna`, `coding`, `coding_transcript_flank`, `coding_gene_flank`, `transcript_flank`, `gene_flank`, `peptide`, `3utr`, `5utr`

`getGene(id, type, mart)`
:   IDを指定して遺伝子情報をダウンロードすることに特化した `getBM()` ラッパー。
    得られるdata.frameは9列:
    指定したID, `external_gene_id`, `description`, `chromosome_name`,
    `band`, `strand`, `start_position`, `end_position`, `ensembl_gene_id`

`select(mart, keys, columns, keytype)`
:   `getBM()` を生のSQLっぽくした感じ。
    keytype(フィルター)は1つしか使えない。
    `dplyr::select()` と名前が衝突する。

`columns(mart)`
:   `mart@attributes$name` と同じ

`keytypes(mart)`
:   `mart@filters$name` と同じ

`keys(mart, keytype)`
:   `subset(mart@filters, name==keytype)$options`
    をちゃんと要素ごとに切った文字列ベクタで。

`exportFASTA(sequences, file)`

## 使用例

### Ensembl

マートとデータセットの決定

```r
library(biomaRt)
ensembl = useMart("ensembl", "scerevisiae_gene_ensembl")

## あるいは
listMarts()
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("scerevisiae_gene_ensembl", ensembl)
```

どんなデータやフィルタが利用可能か調べる

```r
attributePages(ensembl)
listAttributes(ensembl, "feature_page")
subset(ensembl@filters, select=c(name, description, type, operation))
```

フィルタの選択肢を調べる

```r
> subset(ensembl@filters, name=="biotype")$options
[1] "[ncRNA,protein_coding,pseudogene,rRNA,snoRNA,snRNA,tRNA]"
> keys(ensembl, "biotype")
[1] "ncRNA" "protein_coding" "pseudogene" "rRNA" "snoRNA" "snRNA" "tRNA"
> keys(ensembl, "go_evidence_code")
[1] "IBA" "IC"  "IDA" "IEA" "IEP" "IGI" "IMP" "IPI" "ISA" "ISM" "ISS" "NAS" "ND"  "TAS"
```

近いミラーを使う

```r
> listMarts(host="asia.ensembl.org")
               biomart               version
1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 75
2     ENSEMBL_MART_SNP  Ensembl Variation 75
3 ENSEMBL_MART_FUNCGEN Ensembl Regulation 75
4    ENSEMBL_MART_VEGA               Vega 55
5                pride        PRIDE (EBI UK)
> ensembl = useMart("ENSEMBL_MART_SNP", "scerevisiae_snp", host="asia.ensembl.org")
```

### UniProt

<https://www.uniprot.org/>

フィルタ列挙

```r
> unimart = useMart("unimart", "uniprot")
> subset(unimart@filters, select=c(name, description, type, operation))
               name       description    type operation
1  superregnum_name  Superregnum name    list         =
2     proteome_name Complete proteome    list         =
3         accession         Accession    text         =
4      protein_name           Protein    text         =
5    length_greater          Length >    text         >
6    length_smaller          Length <    text         <
7  protein_evidence Protein existence    list         =
8           embl_id          EMBL IDs id_list      =,in
9   arrayexpress_id  ArrayExpress IDs id_list      =,in
10       ensembl_id       Ensembl IDs id_list      =,in
11        pdbsum_id        PDBSum IDs id_list      =,in
12        intact_id        IntAct IDs id_list      =,in
13      interpro_id      InterPro IDs id_list      =,in
14            go_id Gene Ontology IDs id_list      =,in
15        gene_name         Gene name    text         =
16       entry_type        Entry type    list         =
17        organelle         organelle    list         =
18        plasmid_f           Plasmid    text         =
```

フィルタの選択肢

```r
> subset(unimart@filters, 3 < nchar(options) & nchar(options) < 120, select=c(name, options))
               name                                                                                                            options
1  superregnum_name                                                                               [Eukaryota,Bacteria,Archaea,Viruses]
7  protein_evidence [1: Evidence at protein level,2: Evidence at transcript level,3: Inferred from homology,4: Predicted,5: Uncertain]
16       entry_type                                                                                                [Swiss-Prot,TrEMBL]
```

"Complete proteome" の選択肢(すげえ長い)を抜き出す

```r
proteome_name = biomaRt::keys(unimart, "proteome_name")
grep("Sac.* cer.*", proteome_name, value=TRUE)
```

アトリビュート列挙
(ただし `embl_id` 以降の項目はほとんど使えない)

```r
> listAttributes(unimart)
                      name       description
1                accession         Accession
2                     name        Entry name
3             protein_name      Protein name
4                gene_name         Gene name
5                 organism          Organism
6         protein_evidence Protein existence
7               entry_type            Status
8                    go_id             GO ID
9                  go_name           GO name
10  db2go_p__dm_primary_id          GO ID(p)
11 db2go_p__dm_description           GO name
12 db2go_f__dm_description       GO name (F)
13  db2go_f__dm_primary_id         GO ID (F)
14  db2go_c__dm_primary_id         GO ID (C)
15 db2go_c__dm_description       GO name (C)
16                 embl_id          EMBL IDs
17              ensembl_id       Ensembl IDs
18             interpro_id      InterPro IDs
19               pdbsum_id        PDBSum IDs
20                  pdb_id           PDB IDs
21            arrayexpress  ArrayExpress IDs
22                pride_id         PRIDE IDs
23             interact_id        IntAct IDs
24                comments          Comments
25               ec_number         Ec number
26                 keyword           Keyword
27            plasmid_name      Plasmid name
28          organelle_name    organelle name
```

------------------------------------------------------------------------

<a href="https://www.amazon.co.jp/dp/4621062506//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=c58e6e9dc365558cc336d9ea0a2c8a12" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621062506&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4621062506" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320057082//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=b0a0cf1dfe769f34f7544db70a0f6711" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320057082&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320057082" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320123700//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=f7d9a9e5ba94b0fa7da4a582b70da85a" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320123700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320123700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
