+++
title = 'ggbio'
subtitle = "Bioconductor データを ggplot2 的に描画"
tags = ["r", "bioconductor", "graph"]
[menu.main]
  parent = "rstats"
  weight = -31
+++

<http://www.tengfei.name/ggbio/>

Yin, Cook, and Lawrence (2012)
ggbio: an R package for extending the grammar of graphics for genomic data.
*Genome Biology* 13:R77
<http://www.ncbi.nlm.nih.gov/pubmed/22937822>

[Bioconductor]({{< relref "bioconductor.md" >}}) を読み込んでインストール

```r
source("http://bioconductor.org/biocLite.R")
biocLite('ggbio')
library(ggbio)
```

## 使い方

<http://www.tengfei.name/ggbio/docs/man/>

`autoplot(gr, ...)`
:   `xlab, ylab, main`\
    `trancate.gaps=FALSE`\
    `truncate.fun=NULL`\
    `ratio=0.0025`\
    `space.skip=0.1`\
    `legend=TRUE`\
    `geom=NULL`\
    `stat=NULL`\
    `coord=c('default', 'genome', 'truncate_gaps')`\
    `layout=c('linear', 'karyogram', 'circle')`

`plotSingleChrom()`, `plotIdeogram()`, `Ideogram()`

`plotStackedOverview()`, `plotKaryogram()`

`plotGrandLinear()`

------------------------------------------------------------------------

`geom_alignment()`

`geon_arch()`

`geom_arrow()`

`geom_arrowrect()`

`geom_bar()`

`geom_chevron()`

`geom_rect()`

`geom_segment()`

------------------------------------------------------------------------

`layout_karyogram(data, ...)`

`layout_circle()`

------------------------------------------------------------------------

`stat_aggregate`
`stat_bin`
`stat_coverage`
`stat_gene`
`stat_identity`
`stat_mismatch`
`stat_reduce`
`stat_slice`
`stat_stepping`
`stat_table`

------------------------------------------------------------------------

`theme_alignment`
`theme_clear`
`theme_genome`
`theme_noexpand`
`theme_null`
`theme_pack_panels`
`theme_tracks_sunset`

------------------------------------------------------------------------

`tracks(...)`
:   複数のトラックをひとつの図にまとめる\
    `heights`, `xlim`,\
    `xlab`, `main`, `title`, `theme`,\
    `track.plot.color`, `track.bg.color`,\
    `main.height=unit(1.5, 'lines')`,\
    `scale.height=unit(1, 'lines')`,\
    `xlab.height=unit(1.5, 'lines')`,\
    `padding=unit(-1, 'lines')`,\
    `label.bg.color='white'`, `label.bg.fill='gray80'`,\
    `label.text.color='black'`, `label.text.cex=1`,\
    `label.width=unit(2.5, 'lines')`

`zoom()`, `zoom_in()`, `zoom_out()`

`nextView()` , `prevView()`

## `biovizBase`

<http://master.bioconductor.org/packages/release/bioc/html/biovizBase.html>

------------------------------------------------------------------------

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4621062506/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41aQBFtkgBL._SX160_.jpg" alt="RとBioconductorを用いたバイオインフォマティクス" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320057082/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51yBjAPptKL._SX160_.jpg" alt="Rによるバイオインフォマティクスデータ解析 第2版 －Bioconductorを用いたゲノムスケールのデータマイニング－" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320123700/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41aoEmhUR0L._SX160_.jpg" alt="トランスクリプトーム解析 (シリーズ Useful R 7)" /></a>
