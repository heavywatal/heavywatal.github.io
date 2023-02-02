+++
title = 'topGO'
subtitle = "Bioconductor でenrichment解析"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -30
+++

[Gene Ontology]({{< relref "gene_ontology.md" >}}) を利用して、
ある遺伝子セットにどんな機能が多めに含まれているかを解析する。

<https://www.bioconductor.org/packages/release/bioc/html/topGO.html>

Rの中からインストール

```r
BiocManager::install("topGO")
```

## 使い方

ヒト遺伝子にランダムなスコアをつけた架空データ

```r
library(tidyverse)
library(org.Hs.eg.db)
entrez_ids = mappedkeys(org.Hs.egGO)
scores = runif(length(entrez_ids), 0, 1)  # p-value-like
# scores = rnorm(length(entrez_ids), 0, 0.4)  # log2FC-like
named_scores = setNames(scores, entrez_ids)
```

### `topGOdata` を作る

```r
tg_data = new("topGOdata",
    ontology="BP",
    allGenes=named_scores,
    geneSelectionFun=function(x) {x < 0.01},
    # geneSelectionFun=function(x) {abs(x) > 1},
    nodeSize=10,
    annotationFun=annFUN.org,
    mapping="org.Hs.eg.db",
    ID="entrez")
```

`ontology`
:   `BP`, `CC`, `MF` のどれか

`description` (省略可)
:   説明string

`allGenes`
:   遺伝子名を名前とするnamed vector。
    値はp-valueとかlog2FCとかなんとか、好きなスコア。
:   名前は"all"だけど常に全遺伝子である必要はなくて、
    むしろ目的に応じて適切な"background"を選んで渡すべき。

`geneSelectionFun`
:   `allGenes` の数値を引数として今回興味のある遺伝子に `TRUE` を返すような関数。
    スコアが閾値より大きい・小さいとか、上位100までとか。
    `function(p) {p < 0.01}`,
    `function(log2FC) {abs(log2FC) > 1}`,
    `function(score) {rank(score) <= 100L}`
:   必要なのは `statistic="fisher"` のときだけで、なおかつ実際に適用されるのは
    `runTest()` 実行時なのにここで入力させる、という筋の悪いデザイン。
    スコアをそのまま使うはずの `statistic="ks"` などでも要求され、
    Significant genes という謎の結果が計算されてしまう。
    fisherを使わない場合は結果の誤解を防ぐためにも一律 `FALSE` を返す関数にしておくのが安全。
    `function(x) {logical(length(x))}`

`nodeSize`
:   紐付けられた遺伝子の数がこれより少ないGO termを除外。1以上の整数。

[`annotationFun`](https://www.rdocumentation.org/packages/topGO/topics/annFUN)
:   遺伝子IDとGO termを結びつける関数。
:   後述のように `org.**.**.db` パッケージがあるようなメジャー種の遺伝子IDなら `annFUN.org` 。
    チップが登録されているマイクロアレイなら `annFUN.db` 。
    自作マップで頑張るなら `annFUN.gene2GO`, `annFUN.GO2genes`, `annFUN.file` 。

`...`
:   以降は `annotationFun()` に渡す引数。
:   e.g., `annFUN.org(whichOnto, feasibleGenes, mapping, ID="entrez")`

    `mapping`
    :   IDマッピング用のBioConductorパッケージ。
        [BioConductor AnnotationData Packages](https://bioconductor.org/packages/release/data/annotation/)
        から探す。例えばヒトなら `org.Hs.eg.db` 。

    `ID`
    :   `allGenes` に与えた名前の種類。
    :   `entrez`, `genbank`, `alias`, `ensembl`,
        `symbol`, `genename`, `unigene`.
        (`annFUN.org()` の中に書いてある。case-insensitive)


### 検定

```r
whichAlgorithms()
whichTests()
so="increasing"
resClassicFisher = runTest(tg_data, algorithm="classic", statistic="fisher", sortOrder=so)
resElimFisher = runTest(tg_data, algorithm="elim", statistic="fisher", sortOrder=so)
resWeightFisher = runTest(tg_data, algorithm="weight", statistic="fisher", sortOrder=so)
resClassicKS = runTest(tg_data, algorithm="classic", statistic="ks", sortOrder=so)
resElimKS = runTest(tg_data, algorithm="elim", statistic="ks", sortOrder=so)
resWeightKS = runTest(tg_data, algorithm="weight01", statistic="ks", sortOrder=so)
```

`algorithm`
:   `classic`: GO termをそのまま使って計算。偽陽性多め。
:   `elim` [(Alexa et al. 2006)](https://doi.org/10.1093/bioinformatics/btl140):
    DAG上での隣接関係を考慮して補正。
    下位termから検定を始め、有意なものが見つかったらそこに含まれる遺伝子を祖先ノードから除外する。
    classicに比べて上位termの偽陽性が減ってconservative。
    `cutOff` オプションでこの判定基準を 0.01 から変更可能。
:   `weight` (Alexa et al. 2006):
    elimを一般化して少しマイルドにしたような感じ。
    有意な子ノードを多く持つ親ノードは生き残る。
    classicよりconservativeだがelimより取りこぼさない。
    計算が複雑すぎるせいか検定は `fisher` しかサポートされていない。
:   `weight01`: "mixture between the `elim` and the `weight` algorithms"
    とのことだが詳細は不明。
    topGOデフォルトに据えるくらい自信あるんだろうけど。検定は全て可能。
:   `lea`: ドキュメントでも論文でも言及無し。お蔵入りしたプロトタイプ？
:   `parentChild` [(Grossmann et al. 2007)](https://doi.org/10.1093/bioinformatics/btm440):
    classicよりも下位termでの偽陽性が少ない。
    親を複数持つ場合の扱い2つ(union or intercection)のうちどちらを採用してるかは不明。
    検定は `fisher` のみ。

`statistic`
:   `fisher`: 遺伝子の数に基づいて検定。セットに含まれているか否かの二値。
:   `ks` [(Ackermann and Strimmer 2009)](https://doi.org/10.1186/1471-2105-10-47):
    遺伝子のスコアや順位に基づいて検定。
    閾値で遺伝子セットを区切らずに済む。
    (普通のKSだったら上位への偏りだけを見ているとは限らないけどそのあたりは調整済み？)
:   `globaltest` [(Goeman and Bühlmann 2007)](https://doi.org/10.1093/bioinformatics/btm051):
    統計量を挟まず生データに基づいて検定。
    (`topGOdata` にどうやってデータ渡すんだろう？)
:   `t`, `sum`, `ks.ties`: 不明。

`scoreOrder`
:   `allGenes` に与えた値が小さいほど良いP値とかなら `increasing` (デフォルト)。
    興味のある遺伝子で値が高くなるlog2FCのようなスコアなら `decreasing` を指定。
    ドキュメントにあんまりちゃんと載ってない？

結局どれを使うか？
:   `algorithm` はとりあえず最もconservativeで検定も自由な `elim` 。
    もしくは作者Alexaさんを信じて `weight01` 。デフォルト設定は論文にも書きやすい。
:   `statistic` は、遺伝子セットが既に区切ってあるなら `fisher` 、
    DEG解析やら何やらで遺伝子がp値やスコアを持ってるなら `ks` 。
    `globaltest` を使いたい場合は
    [本家 `globaltest` パッケージ](http://bioconductor.org/packages/release/bioc/html/globaltest.html)
    を参照。

結果は `topGOresult` オブジェクト。
計算されるp値は多重検定の補正をしていない生の値。
"Significant genes" は `geneSelectionFun(allGenes)` で `TRUE` になったものの数で、
無関係なはずのKSでも計算されてしまって気持ち悪い。

```r
str(resElimKS)
geneData(resElimKS)
score(resElimKS) |> head()
```

### 解釈・描画

`GenTable()` を使って GO term Over-representation ランキングを表示できる。
複数の `runTest()` 産物を好きな名前で複数並べたりすることもできる。
が、P値が桁の小さい文字列型になってたり勝手に行が削られたりして怖いので使わない。
`GenTable()` の実装を参考に似た形式のより良いテーブルを自分で作る。

```r
# num_nodes=length(tg_data@graph@nodes)  # GO terms >nodeSize
# tg_table = GenTable(tg_data,
#   classicFisher=resClassicFisher,
#   elimFisher=resElimFisher,
#   classicKS=resClassicKS,
#   elimKS=resElimKS,
#   topNodes=num_nodes) |>
#   tibble::as_tibble() |>
#   print()

annoStat = termStat(tg_data, sort(tg_data@graph@nodes)) |>
  tibble::rownames_to_column(var = "GO.ID") |>
  tibble::as_tibble() |>
  dplyr::mutate(Term=topGO:::.getTermsDefinition(terms, ontology(tg_data), 65535L)) |>
  dplyr::relocate(Term, .after=GO.ID) |>
  print()

tg_table = annoStat |> dplyr::mutate(
    classicFisher=score(resClassicFisher, whichGO=GO.ID),
    elimFisher=score(resElimFisher, whichGO=GO.ID),
    weightFisher=score(resWeightFisher, whichGO=GO.ID),
    classicKS=score(resClassicKS, whichGO=GO.ID),
    elimKS=score(resElimKS, whichGO=GO.ID),
    weightKS=score(resWeightKS, whichGO=GO.ID),
  ) |>
  print()
```

あとは dplyr, tidyr, ggplot2 などを使って自由に整形、可視化する。

```r
tg_table |> dplyr::arrange(elimKS)

tg_table |>
  dplyr::select(matches("Fisher$|KS$")) |>
  pairs()

p0 = ggplot(tg_table) +
  geom_point(shape = 16, alpha = 0.5) +
  coord_fixed() +
  theme_minimal()

p0 + aes(classicFisher, classicKS)
p0 + aes(elimFisher, elimKS)
p0 + aes(classicKS, elimKS)
p0 + aes(classicKS, weightKS)
p0 + aes(elimKS, weightKS)
```

---

`Rgraphviz` を使ってDAG描画

```r
# BiocManager::install("Rgraphviz")
showSigOfNodes(tg_data, score(resElimKS), firstSigNodes=6, useInfo="all")
printGraph(tg_data, resElimKS, firstSigNodes=6, fn.prefix="filename", pdfSW=TRUE)
```

前者はプロットだけ、後者はPDFに書き出し。
significant nodeが四角で、赤いほど低いP値。

## GO terms

```r
BPterms = ls(GOBPTerm)
MFterms = ls(GOMFTerm)
CCterms = ls(GOCCTerm)
```
