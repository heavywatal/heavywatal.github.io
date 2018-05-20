+++
title = 'R自学自習の基礎知識'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -99
+++

## 初学者向け講義資料

<ol start="0">
<li><a href="/slides/nagoya2018/0-why-r.html">どうしてRを使うの？</a> (12分)
<li><a href="/slides/nagoya2018/1-basic-r.html">Rの基本</a> (15分)
<li><a href="/slides/nagoya2018/2-ggplot.html">R + ggplot2 — きれいなグラフを簡単に合理的に</a> (25分)
<li><a href="/slides/nagoya2018/3-tidy-data.html">R + tidyverse — 使える形にデータを整える</a> (30分)
</ol>

## インターネットで調べる・尋ねる

- <https://www.rdocumentation.org/>
- <http://r4ds.had.co.nz/> --- R for Data Science (体系的に学びたい人はぜひ通読を)
- [Stack Overflow](https://stackoverflow.com/questions/tagged/r)
- [r-wakalang](https://github.com/tokyor/r-wakalang) --- Slack上の日本語コミュニティ


## RStudioを使う

Rの実行方法はいくつかあるが、これから始めるひとは
[RStudio](https://www.rstudio.com/) を使うとよい。


## スクリプトを保存

R のコンソールにコマンドを打ち込むと、即座に結果が返ってくる。
一度きりで済む短い処理ならそのように対話的な方法が分かりやすいかもしれないが、
別のデータにも使い回すとか、タイプミスでやり直すとか、
同じような処理を何度か繰り返す場合には(大抵はそうなる)、
一連の処理をまとめてスクリプト(テキストファイル)に書き出しておくとよい。
ファイルの拡張子は `.txt` でも何でもいいが `.R` にすることが多い。

[R Markdown](https://rmarkdown.rstudio.com/)
を使えばコードと解析結果の図表を同時に保存・プレゼンすることができる。


## パッケージ

便利な関数やサンプルデータなどをひとまとめにしたもの。
R開発チームが公式に作ってるものから、
ユーザーが自分用に作ってアップロードしたものまで、さまざまある。
自分でやろうとしてることは既に誰かがやってくれてる可能性が高いので、
車輪の再発明をする前に、まずは既存のパッケージを調べるべし。
プログラミングは、イチから「書く」のではなく、
既存のパーツを利用して「組む」感覚のほうが近い。


### Standard Packages

https://stat.ethz.ch/R-manual/R-devel/doc/html/packages.html

R の標準機能。
何もしなくても使用可能な状態になっているので、
パッケージであることはあまり意識しなくてもいい。

`base`
: `c()`, `data.frame()`, `sum()` などホントに基本的なもの

`graphics`, `grDevices`, `grid`
: `plot()` などグラフ描画関連

`stats`
: `anova()`, `glm()`, `t.test()` など統計解析関連

`utils`
: `help()`, `install.packages()`, `read.table()` など

ほかに
`compiler`, `datasets`, `methods`, `parallel`,
`splines`, `stats4`, `tcltk`, `tools`

### Recommended Packages

R と一緒にインストールされるが、
使用する前に `library()` で呼び出しておく必要があるパッケージ。

例えば `MASS` に入ってる `stepAIC()` を使うには

```r
> stepAIC(model)
Error: could not find function "stepAIC"
> library(MASS)
> stepAIC(model)  # OK
```

ほかに
`boot`, `class`, `cluster`, `codetools`,
`foreign`, `KernSmooth`, `lattice`, `Matrix`,
`mgcv`, `nlme`, `nnet`, `rpart`, `spatial`, `survival`

### Contributed Packages

<a href="https://tidyverse.tidyverse.org/">
<img src="https://tidyverse.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

数千ものパッケージが有志により開発され、CRANにまとめて公開されている。

例えば [rstudio.com/products/rpackages](https://www.rstudio.com/products/rpackages/)
で紹介されているもの、特にHadley Wickhamらによる
[tidyverse](https://www.tidyverse.org/) パッケージ群
([ggplot2]({{< relref "ggplot2.md" >}}),
[dplyr]({{< relref "dplyr.md" >}}),
[purrr]({{< relref "purrr.md" >}}),
[tidyr]({{< relref "tidyr.md" >}}),
[readr]({{< relref "readr.md" >}}),
[stringr]({{< relref "stringr.md" >}})など)
はどんな解析にも有用で、標準になってもいいくらい便利。
というより、tidyverseなき裸のRは使う気が起きない。
Rの中から下記のようなコマンドで一括インストール・読み込みできる。

```r
> install.packages('tidyverse')
> library(tidyverse)
```

そのほか
<https://cran.r-project.org/web/views/>
で用途別に紹介されている。

パッケージを作るには [devtools]({{< relref "devtools.md" >}}) を使う。


## 作業ディレクトリ

R ではどこかのディレクトリ (= フォルダ) に身をおいて作業する。
[ファイルを開く]({{< relref "readr.md" >}})ときなどに
R がファイルを探すのはこの作業ディレクトリである。
`No such file or directory` と怒られる場合は、
作業ディレクトリとファイルの場所が合っていないかったり、
ファイル名のタイプミスだったりすることが多い。

`getwd()`, `setwd()`
:   現在地の working directory を get/set する関数

    ```r
    > getwd()
    [1] /Users/watal
    > setwd("~/Desktop")
    > getwd()
    [1] /Users/watal/Desktop
    ```

`list.files(path = ".", ...)` または `dir(path = ".", ...)`
:   ディレクトリ内のファイルを列挙する関数。
    path を省略するとワーキングディレクトリが対象となる

    ```r
    > read.table("mydata.txt", header=TRUE)
    Warning in file(file, "rt") :
      cannot open file 'mydata.txt': No such file or directory
    Error in file(file, "rt") : cannot open the connection
    > list.files()
    [1] mydate.txt
    # ファイル名が微妙に違う!!
    ```

## R で調べる

`help(topic, package=NULL, ...)`
:   ヘルプを表示する。
    知りたいオブジェクトの先頭にクエスチョンマークを付けるという方法もある。
    以下の2つは等価

    ```r
    > help(sum)
    > ?sum
    ```

`help.start()`
:   組み込みヘルプをウェブブラウザで開く。
    インストール済みパッケージのヘルプも Packages のリンクから見られる。

`attributes(obj)`, `str(obj)`
:   オブジェクトの属性や構造を調べる。
    `$names` に列挙されてる属性はダラー `$` を挟んで取り出せる

    ```r
    > x = rnorm(100)
    > y = rnorm(100)
    > mymodel = lm(y ~ x)
    > attributes(mymodel)
    $names
     [1] "coefficients"  "residuals"     "effects"       "rank"
     [5] "fitted.values" "assign"        "qr"            "df.residual"
     [9] "xlevels"       "call"          "terms"         "model"

    $class
    [1] "lm"
    >
    > mymodel$coefficients
    (Intercept)           x
     -0.1494812   0.1218096
    >
    > str(mymodel)
    List of 12
     $ coefficients : Named num [1:2] -0.171 -0.198
      ..- attr(*, "names")= chr [1:2] "(Intercept)" "x"
     $ residuals    : Named num [1:100] -0.5143 -1.3148 0.1954 1.6039 -0.0875 ...
      ..- attr(*, "names")= chr [1:100] "1" "2" "3" "4" ...
     $ effects      : Named num [1:100] 1.2027 2.0571 0.1179 1.8384 -0.0946 ...
      ..- attr(*, "names")= chr [1:100] "(Intercept)" "x" "" "" ...
    ```

## データ読み込み

See [readr]({{< relref "readr.md" >}})

## データ処理・整形

See [dplyr]({{< relref "dplyr.md" >}}), [purrr]({{< relref "purrr.md" >}}), and [tidyr]({{< relref "tidyr.md" >}})

## グラフ作図

2D: See [ggplot2]({{< relref "ggplot2.md" >}})

3D: See [rgl]({{< relref "rgl.md" >}})

## 関連書籍

Rのモダンな使い方、考え方。
まずは[公開オンライン版](http://r4ds.had.co.nz/)を読んでみて。

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=7512fc64786c3444187031221a043ee8" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=76fdb18e30b26131bc43d81752150e58" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />

古典的な仮説検定から一般化線形モデル・階層ベイズモデルくらいまでならこのへんで:

<a href="https://www.amazon.co.jp/%E3%83%87%E3%83%BC%E3%82%BF%E8%A7%A3%E6%9E%90%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AE%E7%B5%B1%E8%A8%88%E3%83%A2%E3%83%87%E3%83%AA%E3%83%B3%E3%82%B0%E5%85%A5%E9%96%80__%E4%B8%80%E8%88%AC%E5%8C%96%E7%B7%9A%E5%BD%A2%E3%83%A2%E3%83%87%E3%83%AB%E3%83%BB%E9%9A%8E%E5%B1%A4%E3%83%99%E3%82%A4%E3%82%BA%E3%83%A2%E3%83%87%E3%83%AB%E3%83%BBMCMC-%E7%A2%BA%E7%8E%87%E3%81%A8%E6%83%85%E5%A0%B1%E3%81%AE%E7%A7%91%E5%AD%A6-%E4%B9%85%E4%BF%9D-%E6%8B%93%E5%BC%A5/dp/400006973X/ref=as_li_ss_il?ie=UTF8&qid=1486144002&sr=8-1&keywords=%E7%B5%B1%E8%A8%88%E3%83%A2%E3%83%87%E3%83%AB&linkCode=li3&tag=heavywatal-22&linkId=bde06ff31718be83e792333c4e10a110" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=400006973X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=400006973X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E5%AD%A6%E3%81%B6%E7%B5%B1%E8%A8%88%E5%AD%A6%E5%85%A5%E9%96%80-%E5%B6%8B%E7%94%B0-%E6%AD%A3%E5%92%8C/dp/4807908596/ref=as_li_ss_il?ie=UTF8&qid=1486144186&sr=8-26&keywords=%E6%95%B0%E7%90%86%E7%B5%B1%E8%A8%88%E5%AD%A6&linkCode=li3&tag=heavywatal-22&linkId=5a14aebbec604ea6e9a6c104dfb50a0d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4807908596&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4807908596" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
