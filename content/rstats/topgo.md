+++
title = 'topGO'
subtitle = "Bioconductor でenrichment解析"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -30
+++

<https://www.bioconductor.org/packages/release/bioc/html/topGO.html>

Rの中からインストール

```r
BiocManager::install("topGO")
```

## 使い方

### `topGOdata` を作る

```r
tg_data = new('topGOdata',
    ontology='BP',
    allGenes=setNames(score, name),
    geneSelectionFun=function(x) {x > threshold},
    nodeSize=10,
    annotationFun=annFUN.org,
    mapping='org.Hs.eg.db',
    ID='entrez')
```

`ontology`
:   `BP`, `CC`, `MF` のどれか

`description` (省略可)
:   説明string

`allGenes`
:   遺伝子名を名前とするnamed vector。
    値はP-valueとかなんとか、好きなスコア。

`geneSelectionFun`
:   `allGenes` の数値を引数として、
    今回興味のある遺伝子に `TRUE` を返すような関数。

`nodeSize`
:   1以上の整数。
    これ以下の遺伝子数しかないGO termを結果から除外。

`annotationFun`
:   遺伝子IDとGO termを結びつける関数。
    一般的なデータベースのIDなら `annFUN.org` で足りるはず。

`...`
:   以降は `annotationFun()` に渡す引数。

------------------------------------------------------------------------

`annFUN.org(whichOnto, feasibleGenes, mapping, ID='entrez')`

`mapping`
:   [Bioconductor]({{< relref "bioconductor.md" >}}) のマッピングパッケージ。
    例えばヒトなら `org.Hs.eg.db` 。

`ID`
:   `allGenes` に与えた名前の種類。
    `entrez`, `genbank`, `alias`, `ensembl`,
    `symbol`, `genename`, `unigene`

### 解析

```r
whichTests()
whichAlgorithms()
tg_result = runTest(tg_data, algorithm='classic', statistic='fisher')
```

`algorithm`
:   `classic`, `elim`, `weight`,
    `weight01`, `lea`, `parentChild`

`statistic`
:   `fisher`, `ks`, `t`, `globaltest`, `sum`

`scoreOrder`
:   デフォルトはP値を扱うように `increasing` 。
    興味のある遺伝子で値が高くなるスコアの場合は `decreasing` を指定。

結果は `topGOresult` オブジェクト

```r
score(tg_result)
geneData(tg_result)
```

### 解釈・描画

GO termのOver-representationランキング

```r
num_significant = geneData(tg_result)['Significant']
GenTable(tg_data, classic_fisher=tg_result, topNodes=num_significant)
```

`runTest()` の結果は好きな名前で複数並べることができる。

------------------------------------------------------------------------

`Rgraphviz` を使ってDAG描画

```r
showSigOfNodes(tg_data, score(tg_result), firstSigNodes=10, useInfo='all')
printGraph(tg_data, tg_result, 20, fn.prefix='go_fisher', pdfSW=TRUE)
```

前者はプロットだけ、後者はPDFに書き出し。
significant nodeが四角で、赤いほど低いP値。

## GO terms

```r
BPterms = ls(GOBPTerm)
MFterms = ls(GOMFTerm)
CCterms = ls(GOCCTerm)
```
