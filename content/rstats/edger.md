+++
title = 'edgeR'
subtitle = "リードカウントから発現変動遺伝子を検出"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -30
+++

<https://bioconductor.org/packages/release/bioc/html/edgeR.html>

Robinson MD, McCarthy DJ, Smyth GK (2010)
"edgeR: a Bioconductor package for differential expression analysis of digital gene expression data."
*Bioinformatics* **26** (1):139--140
<https://www.ncbi.nlm.nih.gov/pubmed/19910308>

[Bioconductor]({{< relref "bioconductor.md" >}}) パッケージとしてインストール

```r
BiocManager::install("edgeR")
```

ユーザーガイドPDFを開く

```r
edgeRUsersGuide()
```

とても良くできたドキュメントなので必読

## 使い方

1.  `HTSeq` などで遺伝子ごとのリードカウントを用意する
1.  Rで読み込む

    ```r
    library(edgeR)

    targets = data.frame(group=c("control", "case"),
                         files=c("control.txt", "case.txt"))
    dge = readDGE(targets, header = FALSE)
    ```

1.  低カウント過ぎる遺伝子を除去

    ```r
    ok_c = (dge$counts > 5) |> rowSums() |> as.logical()
    ok_cpm = (cpm(dge) > 1) |> rowSums() |> as.logical()
    dge = dge[ok_c & ok_cpm, , keep.lib.sizes = FALSE]
    ```

1.  正規化係数を計算

    ```r
    dge = calcNormFactors(dge)
    dge$samples
    ```

1.  モデル

    ```r
    design = model.matrix(~ 0 + group, data = targets)
    ```

1.  common dispersion, trended dispersion, tagwise dispersionを推定

    ```r
    dge = estimateDisp(dge, design)
    ```

    ただしこれにはbiological replicatesが必要。
    1グループ1サンプルずつしか無い場合は4つの選択肢がある。

    1.  検定を諦めてdescriptive discussion (推奨)
    1.  経験的な値を適当に入れとく

        ```r
        bcv = 0.4  # for human
        bcv = 0.1  # for genetically identical model organisms
        bcv = 0.01 # for technical replicates
        dge$common.dispersion = bcv ^ 2
        ```

    1.  説明変数を減らしたモデルでdispersionを推定し、フルモデルで使う

        ```r
        reduced = estimateGLMCommonDisp(dge, reduced.design, method="deviance", robust=TRUE, subset=NULL)
        dge$common.dispersion = reduced$common.dispersion
        ```

        合計カウントが多くてDEGが比較的少ないときに有効。
        当然、1変数2群比較では使えない。

    1.  群間で差がないと考えられるhousekeeping genesから推定。
        100個以上の遺伝子を使うのが望ましい。
        例えばヒトなら

        ```r
        system("wget https://www.tau.ac.il/~elieis/HKG/HK_genes.txt")
        hk_genes = read_tsv("HK_genes.txt", col_names=c("gene_symbol", "refseq"))
        tmp = dge
        tmp$samples$group = 1
        housekeeping = rownames(tmp$counts) %in% hk_genes$gene_symbol
        hk = estimateCommonDisp(tmp[housekeeping,])
        dge$common.dispersion = hk$common.dispersion
        ```

1.  1変数2群ならシンプルに検定

    ```r
    et = exactTest(dge)
    topTags(et, 20, adjust.method="BH", sort.by="PValue")
    ```

    多群ならGLMで尤度比検定

    ```r
    fit = glmFit(dge, design)
    lrt = glmLRT(fit, coef=2:3)

    lrt_1vs2 = glmLRT(fit, coef=2)
    lrt_1vs3 = glmLRT(fit, coef=3)
    lrt_2vs3 = glmLRT(fit, contrast=c(0, -1, 1))
    ```

    検定結果からDEGを抜き出す

    ```r
    min_lfc = 1
    de = decideTestsDGE(et, adjust.method="BH", p.value=0.01, lfc=min_lfc)
    de_tags = rownames(dge)[as.logical(de)]
    ```

    P値とlogFCの両方で切ったほうがいいらしい。

    `topTags()` や `decideTestsDGE()`
    で多重検定の補正をするということは、
    `et` や `lrt` の `PValue` は補正前。

1.  プロット `logFC ~ mean(logCPM)` してみる

    ```r
    plotSmear(et, de.tags=de_tags)
    abline(h=min_lfc * c(-1, 1), col="blue")
    ```

1.  Gene Ontology解析

    ```r
    go = goana(lrt, species="Hs")
    topGO(go)
    ```

    NCBI RefSeqアノテーションを使うのでIDはEntrez系で

## Appendix

レプリケート *i* における遺伝子 *g* の観察リード数を $y _{gi}$、
知りたい真の発現fractionを $\pi _{gi}$ とする。

<div>\[\begin{aligned}
\sum _g \pi _{gi} &= 1 \\
\sqrt {\phi _g} &\equiv \text{CV}[\pi _{gi}]_i
                = \frac {\operatorname{sd}[\pi _{gi}]_i} {\operatorname{mean}[\pi _{gi}]_i} \\
N_i \pi _{gi} &\sim \text{Gamma}(\phi _g^{-1}, \mu _{gi} \phi _g) \\
\operatorname{E}[N_i \pi _{gi}] &= k\theta = \mu _{gi} \\
\operatorname{var}[N_i \pi _{gi}] &= k\theta^2 = \mu _{gi}^2 \phi _g \\
y_{gi} &\sim \text{Poisson}(N_i \pi _{gi}) \\
       &= \int _0^\infty \text{Poisson}(N_i \pi _{gi})~
                         \text{Gamma}(\phi _g^{-1}, \mu _{gi} \phi _g)~ \mathrm d N_i \pi _{gi} \\
       &= \text{NB}(\phi _g^{-1}, \frac {\mu _{gi} \phi _g} {1 + \mu _{gi} \phi _g}) \\
\operatorname{E}[y_{gi}] &= \mu _{gi} \\
\operatorname{var}[y_{gi}] &= \operatorname{E} \left[\operatorname{var}[y_{gi} \mid \pi _{gi}]_i \right] _\pi
                      + \operatorname{var} \left[\operatorname{E}[y_{gi}\mid \pi _{gi}]_i \right] _\pi \\
                     &= \operatorname{E}[N _i \pi _{gi}] _\pi + \operatorname{var}[N _i \pi _{gi}] _\pi \\
                     &= \mu _{gi} + \mu _{gi}^2 \phi _g \\
\text{CV}^2[y_{gi}] &= 1 / \mu _{gi} + \phi _g \\
                   &= \text{CV}^2[y_{gi} \mid \pi _{gi}] + \text{CV}^2[\pi _{gi}] \\
                   &= \text{Technical~CV}^2 + \text{Biological~CV}^2 \\
\end{aligned}\]</div>

dispersion $\phi _g$
:   普通は $D = \sigma^2 / \mu$ と定義されるけど、
    ここでは $\text{CV}^2$

BCV $\sqrt{\phi _g}$
:   biological coefficient of variation.
    こっちは普通と同じように $\text{CV} = \sigma / \mu$ 。
    sequencing depthがどんなに大きくても残るばらつき。

DGE
:   digital gene expression

CPM
:   counts per million

logFC
:   log2-fold-change
