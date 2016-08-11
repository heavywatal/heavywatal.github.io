+++
title = 'rtracklayer'
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -35
+++

<https://www.bioconductor.org/packages/devel/bioc/html/rtracklayer.html>

## 読み書き

<https://genome.ucsc.edu/FAQ/FAQformat.html>

GFF, BED, bedGraph, BED15, WIG, BigWig, 2bit &lt;---&gt; `GRanges`, `GRangesList`

`import(filename, format, text, ...)`
:   ファイルあるいは文字列を読み込んで `GRanges` に変換する。
    `format` は拡張子から判断してくれる。
    `import.gff()`, `import.bed()` など形式ごとの関数も定義されていて、
    引数などの詳しいヘルプもそっちから見られる。

    GFFは `feature.type='exon'` などとして興味のある行だけ読んだほうが良さそう。
    結果として、余分な列も減る。

`export(gr, filename, format, ...)`
:   `GRanges` をGFFなどの形式でファイルに書き出す。

## UCSC Table Browserにアクセス

<https://genome.ucsc.edu/cgi-bin/hgTables>

ゲノム、トラック、テーブルの指定方法が独特

```r
session = browserSession()
genome(session) = 'hg19'
query = ucscTableQuery(session, 'knownGene')
tableName(query) = 'kgXref'
kgXref = getTable(query)  # data.frame

query = ucscTableQuery(session, 'knownGene', 'hg19', 'kgXref')
```

テーブルを指定しなければ先頭が使われる

```r
query = ucscTableQuery(session, 'cytoBand', 'hg19')
cytoBand = track(query)  # GRanges
```

トラックやテーブルの名前を閲覧

```r
trackNames(session)
tableNames(query)
```
