+++
title = 'edgeR'
subtitle = "リードカウントから発現変動遺伝子を検出"
[menu.main]
  parent = "rstats"
+++

<https://bioconductor.org/packages/release/bioc/html/edgeR.html>

Robinson MD, McCarthy DJ, Smyth GK (2010)
"edgeR: a Bioconductor package for differential expression analysis of digital gene expression data."
*Bioinformatics* **26** (1):139--140
<http://www.ncbi.nlm.nih.gov/pubmed/19910308>

[Bioconductor]({{< relref "bioconductor.md" >}}) パッケージとしてインストール

```r
source('http://bioconductor.org/biocLite.R')
biocLite('edgeR')
```

ユーザーガイドPDFを開く

```r
edgeRUsersGuide()
```

とても良くできたドキュメントなので必読

## 使い方

1.  `HTSeq` などで遺伝子ごとのリードカウントを用意する
2.  Rで読み込む

    ```r
    library(edgeR)

    targets = data.frame(group=c('control', 'case'),
                         files=c('control.txt', 'case.txt'))
    dge = readDGE(targets, header=FALSE)
    ```

3.  低カウント過ぎる遺伝子を除去

    ```r
    ok_c = (dge$counts > 5) %>% rowSums %>% {. > 0}
    ok_cpm = dge %>% cpm %>% {. > 1} %>% rowSums %>% {. > 0}
    dge = dge[ok_c & ok_cpm, , keep.lib.sizes=FALSE]
    ```

4.  正規化係数を計算

    ```r
    dge = calcNormFactors(dge)
    dge$samples
    ```

5.  モデル

    ```r
    design = model.matrix(~ 0 + group, data=targets)
    ```

6.  common dispersion, trended dispersion, tagwise dispersionを推定

    ```r
    dge = estimateDisp(dge, design)
    ```

    ただしこれにはbiological replicatesが必要。
    1グループ1サンプルずつしか無い場合は4つの選択肢がある。

    1.  検定を諦めてdescriptive discussion (推奨)
    2.  経験的な値を適当に入れとく

        ```r
        bcv = 0.4  # for human
        bcv = 0.1  # for genetically identical model organisms
        bcv = 0.01 # for technical replicates
        dge$common.dispersion = bcv ^ 2
        ```

    3.  説明変数を減らしたモデルでdispersionを推定し、フルモデルで使う

        ```r
        reduced = estimateGLMCommonDisp(dge, reduced.design, method='deviance', robust=TRUE, subset=NULL)
        dge$common.dispersion = reduced$common.dispersion
        ```

        合計カウントが多くてDEGが比較的少ないときに有効。
        当然、1変数2群比較では使えない。

    4.  群間で差がないと考えられるhousekeeping genesから推定。
        100個以上の遺伝子を使うのが望ましい。
        例えばヒトなら

        ```r
        system('wget http://www.tau.ac.il/~elieis/HKG/HK_genes.txt')
        hk_genes = read_tsv('HK_genes.txt', col_names=c('gene_symbol', 'refseq'))
        tmp = dge
        tmp$samples$group = 1
        housekeeping = rownames(tmp$counts) %in% hk_genes$gene_symbol
        hk = estimateCommonDisp(tmp[housekeeping,])
        dge$common.dispersion = hk$common.dispersion
        ```

7.  1変数2群ならシンプルに検定

    ```r
    et = exactTest(dge)
    topTags(et, 20, adjust.method='BH', sort.by='PValue')
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
    de = decideTestsDGE(et, adjust.method='BH', p.value=0.01, lfc=min_lfc)
    de_tags = rownames(dge)[as.logical(de)]
    ```

    P値とlogFCの両方で切ったほうがいいらしい。

    {{%div class="note"%}}
`topTags()` や `decideTestsDGE()`
で多重検定の補正をするということは、
`et` や `lrt` の `PValue` は補正前。
    {{%/div%}}

8.  プロット `logFC ~ mean(logCPM)` してみる

    ```r
    plotSmear(et, de.tags=de_tags)
    abline(h=min_lfc * c(-1, 1), col="blue")
    ```

9.  Gene Ontology解析

    ```r
    go = goana(lrt, species='Hs')
    topGO(go)
    ```

    NCBI RefSeqアノテーションを使うのでIDはEntrez系で

## Appendix

レプリケート *i* における遺伝子 *g* の観察リード数を `$y _{gi}$` 、
知りたい真の発現fractionを `$\pi _{gi}$` とする。

<div>$$\begin{split}
\sum _g \pi _{gi} &= 1 \\
\sqrt {\phi _g} &\equiv \mathrm{CV}[\pi _{gi}]_i
                = \frac {\mathrm {sd}[\pi _{gi}]_i} {\mathrm{mean}[\pi _{gi}]_i} \\
N_i \pi _{gi} &\sim \mathrm{Gamma}(\phi _g^{-1}, \mu _{gi} \phi _g) \\
\mathrm E[N_i \pi _{gi}] &= k\theta = \mu _{gi} \\
\mathrm{var}[N_i \pi _{gi}] &= k\theta^2 = \mu _{gi}^2 \phi _g \\
y_{gi} &\sim \mathrm{Poisson}(N_i \pi _{gi}) \\
       &= \int _0^\infty \mathrm{Poisson}(N_i \pi _{gi})~
                         \mathrm{Gamma}(\phi _g^{-1}, \mu _{gi} \phi _g)~ \mathrm d N_i \pi _{gi} \\
       &= \mathrm{NB}(\phi _g^{-1}, \frac {\mu _{gi} \phi _g} {1 + \mu _{gi} \phi _g}) \\
\mathrm E[y_{gi}] &= \mu _{gi} \\
\mathrm{var}[y_{gi}] &= \mathrm E \big[\mathrm{var}[y_{gi} \mid \pi _{gi}]_i \big]_\pi
                      + \mathrm{var} \big[\mathrm E[y_{gi}\mid \pi _{gi}]_i \big]_\pi \\
                     &= \mathrm E[N_i \pi _{gi}]_\pi + \mathrm{var}[N_i \pi _{gi}]_\pi \\
                     &= \mu _{gi} + \mu _{gi}^2 \phi _g \\
\mathrm{CV}^2[y_{gi}] &= 1 / \mu _{gi} + \phi _g \\
                   &= \mathrm{CV}^2[y_{gi} \mid \pi _{gi}] + \mathrm{CV}^2[\pi _{gi}] \\
                   &= \mathrm{Technical~CV}^2 + \mathrm{Biological~CV}^2 \\
\end{split}$$</div>

dispersion `$\phi _g$`
:   普通は $D = \sigma^2 / \mu$ と定義されるけど、
    ここでは $\mathrm{CV}^2$

BCV `$\sqrt{\phi _g}$`
:   biological coefficient of variation.
    こっちは普通と同じように $\mathrm{CV} = \sigma / \mu$ 。
    sequencing depthがどんなに大きくても残るばらつき。

DGE
:   digital gene expression

CPM
:   counts per million

logFC
:   log2-fold-change

------------------------------------------------------------------------

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4621062506/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41aQBFtkgBL._SX160_.jpg" alt="RとBioconductorを用いたバイオインフォマティクス" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320057082/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51yBjAPptKL._SX160_.jpg" alt="Rによるバイオインフォマティクスデータ解析 第2版 －Bioconductorを用いたゲノムスケールのデータマイニング－" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320123700/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41aoEmhUR0L._SX160_.jpg" alt="トランスクリプトーム解析 (シリーズ Useful R 7)" /></a>
