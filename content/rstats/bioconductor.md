+++
title = 'Bioconductor'
subtitle = "Genomicデータ解析ツール群"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -40
+++

<a href="https://www.amazon.co.jp/dp/4621062506//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=c58e6e9dc365558cc336d9ea0a2c8a12" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621062506&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4621062506" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320057082//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=b0a0cf1dfe769f34f7544db70a0f6711" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320057082&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320057082" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320123700//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=f7d9a9e5ba94b0fa7da4a582b70da85a" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320123700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320123700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />

- <https://www.bioconductor.org>
- <http://blog.hackingisbelieving.org/2010/09/r-package_06.html>
- <https://qiita.com/tags/bioconductor>
- <http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual>
- <http://search.bioconductor.jp/>

## 基本操作

### インストール

<https://www.bioconductor.org/install/>

Bioconductor 3.8 から
[BiocManager](https://cran.r-project.org/package=BiocManager)
を使う方法に変わった。
`source()` や `BiocLite()` はもう使わない。

1.  [R本体をインストール]({{< relref "intro.md" >}})
1.  Rの中でコマンドを実行:

    ```r
    install.packages("BiocManager")
    BiocManager::install(c("Biostrings", "GenomicRanges", "rtracklayer"))
    ```

### パッケージ管理

標準の関数ではなく
[BiocManager](https://cran.r-project.org/package=BiocManager)
を使う:

```r
# バージョンなどを確認
BiocManager::valid()

# インストール済みのものを更新
BiocManager::install()

# 使いたいパッケージを入れる
BiocManager::install("edgeR")
BiocManager::install("VariantAnnotation")
```

インストールしたものを使うときには普通と同じように読み込む:

```r
library(edgeR)
```

一覧
:   <https://www.bioconductor.org/packages/release/bioc/>\
    <https://www.bioconductor.org/packages/release/data/annotation/>\
    <https://www.bioconductor.org/packages/release/data/experiment/>

検索
:   <https://www.bioconductor.org/packages/release/BiocViews.html>

使い方を調べる

```r
help(package='Biostrings')
browseVignettes(package='Biostrings')
```

## `Biostrings`

<https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html>

### クラス

配列
:   `BString`\
    `DNAString`\
    `RNAString`\
    `AAString`

配列セット e.g. ある遺伝子の複数のcds
:   `BStringSet`\
    `DNAStringSet`\
    `RNAStringSet`\
    `AAStringSet`

配列セットリスト e.g. 遺伝子ごとの配列セット `cdsBy(txdb, by='gene')`
:   `BStringSetList`\
    `DNAStringSetList`\
    `RNAStringSetList`\
    `AAStringSetList`

多重整列
:   `DNAMultipleAlignment`\
    `RNAMultipleAlignment`\
    `AAMultipleAlignment`

ペアワイズ整列
:   `PairwiseAlignments`
:   このコンストラクタを直接呼び出してはいけない。
    `pairwiseAlignment()` 関数を使うべし。

`Views` `(x, start=NULL, end=NULL, width=NULL, names=NULL)`
:   「ある配列のこの領域とこの領域」を表現したいとき、
    配列を切り出すのではなく、その切り取り方のみを保持する。
    `IRanges` みたいな感じ。
    `as.character()` でその部分配列を文字列として取り出せる。

`PDict` `(x, max.mismatch=NA, ...)`
:   正規表現 `TAA|TAG|TGA` のように複数の配列で検索する場合にまとめておく。
    `matchPDict()` とかで使う。

デモデータ

```r
.file = system.file("extdata", "someORF.fa", package="Biostrings")
bss = readDNAStringSet(.file)
```

### 関数

配列変換
:   `reverse(bs)`\
    `complement(bs)`\
    `reverseComplement(bs)`

    `transcribe(dna)` --- *deprecated*\
       3'-鋳型鎖=アンチセンス鎖-5' を引数として相補的な 5'-mRNA-3' を返す。非推奨。

    `translate(bs, genetic.code=GENETIC_CODE, if.fuzzy.codon='error')`\
       翻訳して `AAString` に変換。 引数はDNAでもRNAでもよい。 `genetic.code` は名前付き文字列ベクタで、 組み込みで使えるやつは `getGeneticCode('SGC2')` のように指定できる。 何がどの名前で使えるかは `GENETIC_CODE_TABLE` で確認。

    `codons(bs)`\
       3塩基ずつ切って見る `Views` オブジェクトを返す。

    `subseq(bs, start, end)<-`\
       部分置換

    `xscat(...)`\
       `Bstring` 版 `paste0()`

頻度を数える
:   `alphabetFrequency(x, as.prob=FALSE, ...)`\
    `letterFrequency(x, letters, OR='|', as.prob=FALSE, ...)`\
    `dinucleotideFrequency(x, ...)`\
    `trinucleotideFrequency(x, ...)`\
    `oligonucleotideFrequency(x, width, ...)`\
    `nucleotideFrequencyAt(set_or_views, at)`

コンセンサス
:   `consensusMatrix(set_or_views)`\
    `consensusString(set_or_views, ambiguityMap, threshold, ...)`

配列ファイルを読む、書く (FASTA, FASTQ)
:   `readBStringSet(file, format='fasta', ...)`\
    `readDNAStringSet(file, format='fasta', ...)`\
    `readRNAStringSet(file, format='fasta', ...)`\
    `readAAStringSet(file, format='fasta', ...)`\
    `writeXStringSet(x, filepath, append=FALSE, compress=FALSE, format='fasta', ...)`

アラインメントを読む (FASTA, stockholm, clustal)
:   `readDNAMultipleAlignment()`\
    `readRNAMultipleAlignment()`\
    `readAAMultipleAlignment()`

配列を読み込まずに情報だけを見る
:   `fasta.info(file, ...)`\
    `fastq.geometry(file, ...)`

検索
:   `matchPattern(pattern, subject, ...)`\
    `matchPDict(PDict, subject, ...)`\
    `matchPWM(pwm, subject, min.score='80%', widh.score=FALSE, ...)`\
    対象は `Bstring` だけでなく `Views` でもよい。 例えば `codons()` の結果を対象とすれば読み枠限定の検索となる。 結果の返し方の違う `vmatchXXX()` と `countXXX()` もある。

## `GenomicRanges`

<https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html>

<https://qiita.com/yuifu/items/556af5d4d086c96ec783>

### クラス

`GRanges`
:   `@seqnames`\
    `@ranges`\
    `@strand`\
    `@elementMetadata`\
    `@seqinfo`\
    `@metadata`

`makeGRangesFromDataFrame(df, keep.extra.columns=FALSE, ignore.strand=FALSE, seqinfo=NULL, ...)`
:   対応する列名 `seqnames.field`, `start.field`, `end.field`, `strand.field`
    を指定することもできるが、大概いい感じに拾ってくれる。

`GRangesList` `(gr1, gr2, ...)`

`Seqinfo` `(seqnames, seqlengths, isCircular, genome)`
:   `@genome`\
    `@is_circular`\
    `@seqlengths`\
    `@seqnames`

### 関数

配列の情報を get, set, modify
:   `seqinfo(x)`\
    `seqnames(x)`\
    `seqlevels(x)`\
       特定の染色体だけ対象にしたい場合は\
       `seqlevels(x, force=TRUE) = paste0('chr', c(1:22, 'X', 'Y'))`\
       のように削り `seqlevels0(x)` で戻せる\
    `seqlengths(x)`\
    `isCircular(x)`\
    `genome(x)`

個々の区間のプロパティを参照
:   `start(gr)`, `end(gr)`, `width(gr)`\
    `strand(gr)`\
    `names(gr)`\
    `metadata(gr)`

区間の集合を操作 (デフォルトではstrand毎に)
:   <https://qiita.com/wakuteka/items/9634e5ed96db3536756f>

    `range(gr, ..., ignore.strand=FALSE, na.rm=FALSE)`\
       端から端まで1つの区間として返す

    `reduce(gr, ...)`\
       重なってる区間をつなげる

    `gaps(gr, start=1, end=seqlengths(gr))`\
       含まれていない部分だけ抽出 (= 全区間 - `reduce(gr)`)

    `disjoin(gr, ...)`\
       重ならない部分だけ抽出

    `isDisjoint(gr, ...)`\
    `disjointBins(gr, ...)`\
    `coverage(gr, shift=0, width=NULL, wight=1, method=c('auto', 'sort', 'hash'))`

個々の区間を操作
:   `shift(gr, shift=0, ...)`\
    `narrow(gr, start=NA, end=NA, width=NA, ...)`\
    `flank(gr, width, start=TRUE, both=FALSE, ...)`\
    `promoters(gr, upstream=2000, downstream=200, ...)`\
    `resize(gr, width, fix='start', ...)`\
    `restrict(gr, start=NA, end=NA, ...)`\
    `trim(gr, ...)`

重なり <https://qiita.com/wakuteka/items/10027edccc6c2e244cd2>
:   `findOverlaps(query, subject, ...)`\
    `countOverlaps(query, subject, ...)`\
    `overlapsAny(query, subject, ...)`\
    `subsetByOverlaps(query, subject, ...)`

## `GenomicFeatures`

<https://www.bioconductor.org/packages/release/bioc/html/GenomicFeatures.html>

<https://qiita.com/yuifu/items/4bab5f713aa75bd18a84>

`TranscriptDB` からいろんな条件で絞り込み、
該当する区間を `GRanges` または `GRangesList` で返す。

絞り込んで `GRanges` を返す
:   `transcripts(txdb, vals=NULL, columns)`\
       `vals=list(tx_chrom=c('chrI', 'chrV'))` のように指定する。 取りうる値: gene\_id, tx\_id, tx\_name, tx\_chrom, tx\_strand, exon\_id, exon\_name, exon\_chrom, exon\_strand, cds\_id, cds\_name, cds\_chrom, cds\_strand, exon\_rank\
    `exons(txdb, ...)`\
    `cds(txdb, ...)`\
    `genes(txdb, ...)`\
    `promoters(txdb, upstream=2000, downstream=200, ...)`\
    `disjointExons(txdb, aggregateGenes=FALSE, includeTranscripts=TRUE, ...)`\
    `microRNAs(txdb)`\
    `tRNAs(txdb)`\
       `library(FDb.UCSC.tRNAs)` の情報を使ってtRNAコーディング部分を抽出

グループ毎に絞り込んで `GRangesList` を返す
:   `transcriptsBy(txdb, by, use.names=FALSE, ...)`\
       `by` は `'gene'`, `'exon'`, `'cds'`, `'tx'` のどれか。\
    `exonsBy(txdb, by, ...)`\
    `cdsBy(txdb, by, ...)`\
    `intronsByTranscript(txdb, ...)`\
    `fiveUTRsByTranscript(txdb, ...)`\
    `threeUTRsByTranscript(txdb, ...)`

`extractTranscriptSeqs(s, transcripts)`
:   `DNAString` または `BSgenome` から配列を抜き出す。
    `transcripts` 引数は `GRanges`, `GRangesList`, `TranscriptDb` のどれか。

`TranscriptDB` の形式変換
:   `makeTranscriptDbFromBiomart(biomart='ensembl', dataset='hsapiens_gene_ensembl', ...)`\
    `makeTranscriptDbFromUCSC(genome='hg18', table='knownGene', ...)`\
    `makeTranscriptDBFromGFF(file, format=c('gff3', 'gtf'), ...)`\
    `asBED(txdb)`\
    `asGFF(txdb)`

### `TranscriptDB`

データのインストールと読み込み

```r
BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")

library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
txdb = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
```

## `BSgenome`

<https://www.bioconductor.org/packages/release/bioc/html/BSgenome.html>

インストール、利用

```r
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")

library("BSgenome.Scerevisiae.UCSC.sacCer3")
bsg = BSgenome.Scerevisiae.UCSC.sacCer3
```

### クラス

`BSgenome`
:   `@organism`\
    `@species`\
    `@provider`\
    `@provider_version`\
    `@release_date`\
    `@release_name`\
    `@source_url`\
    `@seqinfo`\
    `@user_seqnames`\
    `@masks`\
    `@single_sequences`\
    `@multiple_sequences`\
    `@pkgname` ほか

染色体(`DNAString`)にアクセス

```r
bsg$chrI
bsg[['chrI']]
bsg[[1]]
```

------------------------------------------------------------------------

`BSParams`
:   `@X`\
    `@FUN`\
    `@exclude`\
    `@simplify`\
    `@maskList`\
    `@motifList`\
    `@userMask`\
    `@invertUserMask`

`new()()` で作って `bsapply()()` に渡すらしい

------------------------------------------------------------------------

`GenomeData`

`GenomeDescription`

### 関数

アクセサー
:   `organism(bsg)`\
    `species(bsg)`\
    `provider(bsg)`\
    `providerVersion(bsg)`\
    `releaseDate(bsg)`\
    `releaseName(bsg)`\
    `bsgenomeName(bsg)`\
    `seqlengths()` などseqinfo系関数も使える

利用可能なデータを調べる
:   `available.genomes(splitNameParts=FALSE, type=getOption('pkgType'))`\
    `installed.genomes(splitNameParts=FALSE)`

`getBSgenome(name, masked=FALSE)`
:   インストール済みデータから読み込み。
    `BSgenome.Scerevisiae.UCSC.sacCer3` あるいは `sacCer3` のような名前で。

`getSeq(bsg, names, start=NA, end=NA, width=NA, starnd='+', as.character=FALSE)`
:   `BSgenome` オブジェクトから配列を抜き出す。
    `names` は配列名の文字列ベクタか `GRanges` か `GRangesList`

### データパッケージを作る

<https://www.bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf>

1.  染色体ごとのFASTAファイルを用意する e.g. Ensemblから `*.chromosome.*.fa.gz` をダウンロードして展開

    染色体名と拡張子だけを残したファイル名にする必要がある。
    しかもgzipのままでは読めない残念な仕様なので展開しておく。
    e.g. `IV.fa`

1.  既存の `BSgenome` のやつを参考に適当にseedファイルを作る:

        Package: BSgenome.Scerevisiae.EF4.74.ensembl
        Title: BSgenome.Scerevisiae.EF4.74.ensembl
        Description: BSgenome.Scerevisiae.EF4.74.ensembl
        Version: 0.74.1
        organism: Saccharomyces_cerevisiae
        species: Scerevisiae
        provider: ENSEMBL
        provider_version: release-74
        release_date: 2013-11-25T21:47:17
        release_name: EF4.74
        source_url: ftp://ftp.ensembl.org/pub/release-74/fasta/Saccharomyces_cerevisiae/dna/
        organism_biocview: Saccharomyces_cerevisiae
        BSgenomeObjname: BSgenome.Scerevisiae.EF4.74.ensembl
        seqnames: c("I", "II", "III", "IV", "IX", "Mito", "V", "VI", "VII", "VIII", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI")
        circ_seqs: c("Mito")
        SrcDataFiles1: Saccharomyces_cerevisiae.EF4.74.dna.chromosome.I.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.II.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.III.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.IV.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.IX.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.Mito.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.V.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.VI.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.VII.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.VIII.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.X.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XI.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XII.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XIII.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XIV.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XV.fa.gz, Saccharomyces_cerevisiae.EF4.74.dna.chromosome.XVI.fa.gz
                from ftp://ftp.ensembl.org/pub/release-74/fasta/Saccharomyces_cerevisiae/dna/
        PkgExamples: genome$I  # same as genome[["I"]]
        seqs_srcdir: /Users/watal/db/ensembl/release-74/fasta/saccharomyces_cerevisiae/dna

    パッケージ名に使える文字は結構限られてるので注意

1.  R で `forgeBSgenomeDataPkg(seedfile)` を実行
1.  できあがったパッケージディレクトリをビルドしてインストール:

        R CMD build BSgenome.Scerevisiae.EF4.74.ensembl
        R CMD INSTALL BSgenome.Scerevisiae.EF4.74.ensembl

## `VariantAnnotation`

<https://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html>

### クラス

`VCF`, `CollapsedVCF`, `ExpandedVCF`
:   `@assays`\
    `@colData`\
    `@exptData`\
    `@fixed` `DataFrame`\
       `$REF`: 参照配列の塩基 `ref()`\
       `$ALT`: 変異配列の塩基 `alt()`\
       `$QUAL`\
       `$FILTER`\
    `@info` `DataFrame`\
       `$TSA`: SNV, deletion, insertion\
       `$VE`: `*_incl_consequences.vcf.gz` の追加情報\
    `@rowData`\
       変異の位置情報 `GRanges`。 `rowData()` でアクセスすると `@fixed` の情報込みで表示される。

### 関数

`readVcf(file, genome, param)`
:   予めターミナルで `tabix -h some_file.vcf.gz` を実行して
    indexファイル `some_file.vcf.gz.tbi` を作っておく。
    file には生のファイル名ではなく `Rsamtools::TabixFile` を渡す。
    genome は `sacCer3` みたいな名前か `Seqinfo` を指定。
    param に `GRanges` などを入れると範囲限定で読み込む。
    染色体の名前(`seqlevels`)が合ってないと怒られるので修正する。

`writeVcf(obj, filename, index=FALSE)`
:   `index=TRUE` とすると `bgzip` 圧縮と `tabix` 生成を試みるが失敗する。

`predictCoding(query, subject, seqSource, varAllele, ..., ignore.strand=FALSE)`
:   それぞれの変異がCDSやprotein上のどの位置でどういう変化を引き起こすか というのを計算して `GRanges` のmetadataに書き出す。 (引数の型判別がバグってるっぽいので `vcf@rowData, txdb, bsgenome, alt(vcf)` とする)\
    query : `VCF` あるいは `rowData(vcf)`\
    subject : `TxDb` オブジェクト\
    seqSource : `BSgenome` オブジェクト\
    varAllele : 省略あるいは `alt(vcf)`

## データ読み込み、取得

[rtracklayer]({{< relref "rtracklayer.md" >}})

[biomart]({{< relref "biomart.md" >}})

[GEOquery](https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html)

## モチーフ検索

<http://blog.hackingisbelieving.org/2012/02/dna-bioconductor.html>

[/bio/motif]({{< relref "motif.md" >}})

## 作図

<http://blog.hackingisbelieving.org/2012/02/rbioconductor.html>

<http://blog.hackingisbelieving.org/2012/02/r.html>

<https://qiita.com/wakuteka/items/a99d5fb9f24367f55461>

[ggbio]({{< relref "ggbio.md" >}})
