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
<https://www.ncbi.nlm.nih.gov/pubmed/22937822>

[Bioconductor]({{< relref "bioconductor.md" >}}) を読み込んでインストール

```r
source("https://bioconductor.org/biocLite.R")
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

<https://master.bioconductor.org/packages/release/bioc/html/biovizBase.html>

------------------------------------------------------------------------

<a href="https://www.amazon.co.jp/dp/4621062506//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=c58e6e9dc365558cc336d9ea0a2c8a12" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621062506&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4621062506" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320057082//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=b0a0cf1dfe769f34f7544db70a0f6711" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320057082&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320057082" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4320123700//ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=f7d9a9e5ba94b0fa7da4a582b70da85a" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4320123700&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4320123700" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
