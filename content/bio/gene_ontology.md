+++
title = 'Gene Ontology'
[menu.main]
  parent = "bio"
+++

生物種や分野によらない共通の語彙で遺伝子産物の機能を記述するための用語体系。

機能に応じて遺伝子にたくさんのタグ(**GO term**)を付けましょうってこと。

e.g. Human RB1
<http://amigo.geneontology.org/amigo/gene_product/UniProtKB:P06400>

## GO term

<http://geneontology.org/page/ontology-structure>

e.g. <http://amigo.geneontology.org/amigo/term/GO:0043065#display-graphics-tab>

小さくて専門的な下位termから、大きくて一般的な上位termに向かう
directed acyclic graph (DAG) を構成している。
かなり細かい多層構造。
分岐したあと上位でまた合流するので、木構造での表示には限界がある。

最上位のtermは3つ:

[Biological Process](http://geneontology.org/page/biological-process-ontology-guidelines)
:   生物学的な機能。
    下位には例えば、分化、細胞分裂、細胞死など。
    `GO:0008150`

[Cellular Component](http://geneontology.org/page/cellular-component-ontology-guidelines)
:   細胞内での局在。
    下位には例えば、核内、小胞体など。
    `GO:0005575`

[Molecular Function](http://geneontology.org/page/molecular-function-ontology-guidelines)
:   化学的な機能。
    下位には例えば、加水分解酵素、DNA結合など。
    `GO:0003674`

ほかのアノテーションと同じようにEvidenceのレベルもいろいろある。

### GO relation

<http://geneontology.org/page/ontology-relations>

DAGのエッジのことをrelationと呼ぶ。

-   *is a*
-   *part of*
-   *regulates*: *positively* or *negatively*

フィードフォワード的なショートカットも推定される。。

## GO解析

何らかの解析で遺伝子のサブセットが抽出されたとして、
それらが持つGO termの組成を全ゲノムでの組成と比較してみたときに、
特定の機能群におけるenrichmentが見られるかどうか？

e.g. 対照区と処理区で発現が変化した遺伝子には、XXXというtermが多く含まれる

超幾何分布

## ツール

### Amigo

<http://amigo.geneontology.org>

### BioConductor

[GO.db](http://www.bioconductor.org/packages/release/data/annotation/html/GO.db.html)
:   GO term/relationの構造データ。遺伝子のデータではない。

[GOstats](http://www.bioconductor.org/packages/release/bioc/html/GOstats.html)
:   enrichment解析。あまり人気がないっぽい。

[topGO](http://www.bioconductor.org/packages/release/bioc/html/topGO.html)
:   enrichment解析

[GOexpress](http://www.bioconductor.org/packages/release/bioc/html/GOexpress.html)
:   発現パターンで遺伝子をクラスタリングしてGO解析

[GOSemSim](http://www.bioconductor.org/packages/release/bioc/html/GOSemSim.html)
:   2つの遺伝子群の機能的類似度をGOベースで

[goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html)
:   RNA-seqからGOまで直行？
