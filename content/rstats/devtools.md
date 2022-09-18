+++
title = 'devtools'
subtitle = "Rパッケージ作成支援"
tags = ["r", "tidyverse", "package"]
[menu.main]
  parent = "rstats"
  weight = -50
+++

自分で書いた関数が多くなってきたら、まとめてパッケージを作るとよい。
少しだけ面倒だが、以下のようなメリットがある。

-   普通に変数や関数を定義するとワークスペースが名前でいっぱいになってしまうが、
    パッケージ内で定義されている変数や関数は `ls()` で出てこないのでスッキリ
-   既存のオブジェクトと名前が衝突するような場合でも、
    `mypackage::func1` のように名前空間を明示的に指定して呼び出せる

CRANに上げる程ではないとしても、
GitHubに公開しておけば誰でも使えるようになるので、
共同研究者と解析環境を共有したり、
ひとつの論文のワークフローを置いておいたり、いろいろ使い道はある。

## Rパッケージ

- <https://r-pkgs.org/>
- [How to develop good R packages (for open science)](https://masalmon.eu/2017/12/11/goodrpackages/)

### 最低限の作成手順

1.  開発支援パッケージをインストールする:
    `install.packages(c("devtools", "usethis"))`

1.  [usethis](https://usethis.r-lib.org/) の関数をいくつか使って骨組みを作る:

    ```r
    usethis::create_package("hello")
    usethis::use_mit_license()
    usethis::use_roxygen_md()
    usethis::use_package_doc()
    ```

1.  `devtools::check()` で様子を見てみる。
    LaTeX 関連で怒られたら足りないパッケージを入れる:
    `sudo tlmgr install inconsolata helvetic`

1.  ローカルgitリポジトリを作って最初のコミットをする:
    `usethis::use_git()`

1.  GitHubかどこかに空のリポジトリを作る。
    パッケージと同名でなくてもよい。

1.  そのリモートリポジトリにプッシュ:

    ```sh
    git remote add origin https::/github.com/heavywatal/rhello.git
    git push -u origin master
    ```

1.  とりあえず誰でもインストール可能なパッケージができたはず:
    ```r
    remotes::install_github("heavywatal/rhello")
    ```


### ソース構造

https://r-pkgs.org/package-structure-state.html

```sh
DESCRIPTION  # 一番大事
NAMESPACE    # 見せるオブジェクトを列挙
README.md    # 全体の説明を簡単に
R/           # Rソースコード
data/        # サンプルデータなど
inst/        # CITATION
man/         # ヘルプファイル.Rd
src/         # C++ソースコード
tests/
vignettes/
```

CRANから落としてくる `.tar.gz` ソースコード (bundle package) とは違う。


### [`DESCRIPTION`](https://r-pkgs.org/description.html)

-   どうでも良さそうなファイル名とは裏腹に、ちゃんと書かないと動かない。
-   始めはusethisとかに生成してもらい、他のパッケージを参考にしつつ修正していく。
-   `Imports` に列挙したものは依存パッケージとして同時にインストールされる。
    しかし `library()` 時に同時に読み込まれるわけではない。
-   `Title` はピリオドを含まない一文でタイトルケース。
    `Description` はピリオドを含む一段落。
-   ライセンスを別ファイルにする場合は `License: file LICENSE` と書く
-   `Authors@R` のとこは後でRで評価されるので変な形。

### [`NAMESPACE`](https://r-pkgs.org/namespace.html)

-   [後述のroxygen2](#roxygen2)がソースコードから自動生成するので**直接触らない**。
-   ここで `export()` された関数だけがユーザーから見える。
-   外部パッケージから `importFrom(package, function)`
    された関数はattachされてパッケージ内で利用可能になる。
    `import(package)` で全ての関数をまとめて処理できるけど名前の衝突が怖い。

### [`R/` ソースコード](https://r-pkgs.org/r.html)

-   `NAMESPACE` や `man/*.Rd` を自動生成してもらえるように
    [後述のroxygen](#roxygen2)形式でコメントを書く
-   ファイルの数や名前は何でもいいので、開発者が分かりやすいようにしとく。
    例えば、似たような関数をひとつのファイルにまとめ、
    同名の`@rdname`を指定して`man/`と同じ構造にするとか。
-   `library()` や `require()` を書かない。
    必要なパッケージは `DESCRIPTION` の `Imports` に書き、
    `名前空間::関数()` のようにして使う。
    どうしても名前空間を省略したいものだけ `@importFrom` を指定する。
-   パッケージ読み込み時などに実行される特殊関数
    `.onLoad(...)`, `.onAttach(...)`, `.onUnload(...)`
    は慣習的に`zzz.R`というファイルに記述する。

### [`src/` C++ソースコード](https://r-pkgs.org/src.html)

[Rcppページの"Rパッケージで使う"セクション]({{< relref "rcpp.md#Rパッケージで使う" >}})を参照


### [`vignettes/`](https://r-pkgs.org/vignettes.html)

個々の関数の使用例はRソースファイルの `@examples` に書くとして、
複数の関数を組み合わせてどう使うかとか、
パッケージ全体の使い方とかを説明するのが`vignettes/`の役割。
Rmarkdown形式で書いてHTMLやPDFで表示できるので表現力豊か。

`usethis::use_vignette("hello")` で雛形を作ってもらうのが楽。

`pandoc` と `pandoc-citeproc` が必要なので
[Homebrew]({{< relref "homebrew.md" >}}) とかでインストールしておく。

`check()`がデフォルトで`vignette=TRUE`かつ処理がやや重いので、
毎回その重さを受け入れるか、わざわざ`FALSE`指定するかというのは悩みどころ。

[pkgdown](https://pkgdown.r-lib.org)でウェブサイトを構築すると、
ここに置いてあるvignettesはArticlesという位置づけになる。


### [`inst/`](https://r-pkgs.org/inst.html)

ここに入ってるものはインストール先でトップディレクトリに移される謎仕様。

論文で引用されることを想定している場合は `inst/CITATION` を作る。
`citation("ggplot2")` のように参照できるようになる。


### [`tests/`](https://r-pkgs.org/tests.html)

<a href="https://testthat.r-lib.org/">
<img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/SVG/testthat.svg" align="right" width="120" height="139">
</a>

[testthat](https://testthat.r-lib.org)パッケージを使うのがデファクトスタンダード。
`use_testthat()` で初期設定して
`use_test("somefunction")` のようにテストファイルを追加する。
`tests/testthat/` 以下のファイル構成は `R/` と同じにしておくとわかりやすい。

さらに[covr](https://covr.r-lib.org/)パッケージを使って、
ソースコードのうちどれくらいがテストでカバーされてるかを可視化すると良い。

### `data/`

- <https://r-pkgs.org/data.html>
- <https://usethis.r-lib.org/reference/use_data.html>

`usethis::use_data_raw()` で `data-raw/<dataset>.R` をセットアップし、
その中で `usethis::use_data()` を呼んで
`data/<dataset>.rda` や `R/sysdata.rda` に配置する。

`data/`
: R から `save()`, `load()` で読むようなバイナリ形式のファイルを置く。
  1オブジェクト1ファイルで同名にして `.rda`, `.RData` 拡張子を付ける。
: 勝手にexportされるので `R/` 内のソースにroxygenドキュメントを書く。

`R/sysdata.rda`
: ユーザーに公開せずパッケージ内部で使うためのデータ。

`data-raw/`
: `data/` のファイルを作るためのソースコードやテキストデータを置いておく。
: bundled package には入れない
  (ように `usethis::use_data_raw()` が `.Rbuildignore` を設定する)。

`inst/extdata/`
: データ読み書きを例示するためのデータファイルを置いておく。
: `system.file("extdata", "mtcars.csv", package = "readr")`
  のようにアクセスできる。


### [その他](https://r-pkgs.org/misc.html)

`demo/`
:   vignettesに取って代わられた古い機能。
    ソースコード`*.R`を置いておくと`demo()`関数で呼び出せるというだけ。
:   `check()`でソースコードの中身は実行されないが、
    `demo/00Index` というファイルとの整合性が取れてないと警告される。
    「デモの名前 + スペース3つ以上かタブ + 適当な説明」という形式。
    ファイル名から拡張子を取り除いたものがデモの名前になる。

    ```
    mydemo1    Description of demo 1
    mydemo2    Description of demo 2
    ```

`exec/`
: 実行可能スクリプトの置き場所。
  インストール先はそのままパッケージ内の `exec/`
  (つまり最初からパスが通ってるような場所ではない)。
  パーミッションはインストール時に設定してもらえるのでソースツリーでは `644` でOK。
  `Rscript` をshebangで指定することも可能。
  例えばRStudioを使いこなせないターミナル勢としては
  [`knitr::knit()` するだけのコマンド](https://github.com/heavywatal/rwtl/blob/master/exec/knit.R)
  とか作っておくと便利。

`tools/`
: たぶん何を入れても一切チェックされずインストールもされない自由なディレクトリ。
  `tests/` や `vignettes/` に入れる前のガラクタコードを一時的に置いとくとか。

`.Rbuildignore`
: Rパッケージらしからぬ変なファイルやディレクトリがあると怒られるので、
  そういうやつをこのファイルに列挙して無視してもらう。


## `devtools`

<a href="https://devtools.r-lib.org/">
<img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/SVG/devtools.svg" align="right" width="120" height="139">
</a>

骨組みを作るとこからCRANにデプロイするとこまでお世話してくれる。
いろんな専門パッケージの集合体。

### 主な関数

`document(pkg = ".", roclets = NULL, quiet = FALSE)`
:   `roxygen2` を呼び出してソースコードから
    `NAMESPACE` や `man/*.Rd` を自動生成する。
:   `src/` 内のドキュメント変更はこれ1回の実行では反映されない。

`check(pkg = ".", document = NA, ...)`
:   パッケージとしての整合性を確認。
    ついでに `document()` は実行できるけど `spell_check()` はできないので手動で。

`test(pkg = ".", filter = NULL, ...)`
:   `testthat` を呼び出して `tests/` 以下のテストコードを実行する

`build(pkg = ".", path = NULL, ...)`
:   `R CMD install` でインストール可能な tar ball (bundle package) を作る。
    Rcppコードをコンパイルするという意味ではない。

`install(pkg = ".", reload = TRUE, quick = FALSE, build = !quick, ...)`
:   ローカルにあるソースからインストール。
    `build = TRUE` のとき(デフォルト)、わざわざ bundle package を
    `tempdir()` に作ってからそいつでインストールする。

`install_github(repo, ref = "master", subdir = NULL, ...)`
:   GitHubリポジトリからインストール。

`unload(pkg = ".", quiet = FALSE)`
:   `datach("package:XXX")` とか `unloadNamespace(XXX)`
    よりもちゃんとまっさらにパッケージを外す。


`load_all(pkg = ".", reset = TRUE, recompile = FALSE, export_all = TRUE, ...)`
:   `install()` せずファイルから直接 `library()` する。
    ロード済みでもまず `unload()` が呼ばれるので安心。

`clean_dll(pkg = ".")`
:   `src/` 以下に生成される `.o`, `.so` を消す。
    普段は触る必要ないが、たまにこれが必要な不具合に出くわす。


## `roxygen2`

<a href="https://github.com/klutometis/roxygen">
<img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/SVG/roxygen2.svg" align="right" width="120" height="139">
</a>

Rソースコードのコメントから`NAMESPACE`とヘルプ(`man/*.Rd`)を自動生成する。

- <https://cran.r-project.org/web/packages/roxygen2/>
- <https://github.com/klutometis/roxygen>
- <https://r-pkgs.org/man.html>
- <https://kbroman.org/pkg_primer/pages/docs.html>
- <https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html>

`roxygen2::roxygenise(package.dir=".", ..., clean=FALSE)`
を直接呼んでもよいが、
基本的には `devtools::document()` を使って間接的に利用する。

### 使い方

```r
#' Title of the simple function to add 1
#'
#' The second paragraph is recognized as the description.
#' It is not recommended to use @@title and @@description explicitly.
#' @param x A numeric vector
#' @export
#' @examples
#' increment(42)
increment = function(x) {x + 1}
```

-   `#' ` から始まる行がroxygenコメントとして扱われる。
-   タグは `@` で始まる。
    `@` そのものを入力したいときは重ねて `@@` とする。
-   1行目にタイトル。1行あけて2段落目に説明文を書く。
    明示的に `@title`, `@description` タグを使うことも可能だが非推奨らしい。
    タイトルをコピペして全部同じにすると怒られる。
    2段落目を省略するとタイトルが流用される。
-   空行だけでは切れ目として扱われないので `NULL` などを置いたりする。
-   [Rd形式の代わりにMarkdown形式で記述できる](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html)。
    `usethis::use_roxygen_md()` を実行して
    `DESCRIPTION` に `Roxygen: list(markdown = TRUE)` を書き加える全体設定か、
    使いたいブロックにいちいち `@md` を書く個別設定か。
-   `"_PACKAGE"` という文字列の上に書かれたブロックは、
    パッケージそのものに対応するヘルプ `*-package.Rd` になる。
    `@useDynLib` など全体に関わる設定はここでやるのが良いのでは。

    ```r
    #' Example package to say hello
    #' @aliases NULL hello-package
    #' @useDynLib hello, .registration = TRUE
    #' @importFrom Rcpp sourceCpp
    #' @importFrom magrittr %>%
    #' @importFrom rlang .data
    #' @keywords internal
    "_PACKAGE"
    ```

-   dplyrなどで列名を直に指定すると
    `undefined global variables`
    という警告が出るので
    `dplyr::filter(diamonds, .data$cut == "Ideal")`
    のように pronoun を使って抑える。
    そのためにどこかに `#' @importFrom rlang .data` を書いておく。
    ただし `$` によるアクセスがやや遅いので、
    `group_by()` などで多数のグループを処理するときなど気になる場合は
    `!!as.name("carat")` のようにする。


### タグ

使用可能なタグ一覧は[準備中？](https://github.com/klutometis/roxygen/issues/792)

`@import pkg1, pkg2, ...`
:   `NAMESPACE` で `import()` するパッケージを指定。
    名前の衝突が怖いので基本的に使わない。

`@importFrom pkg func`
:   `NAMESPACE` で `importFrom()` するパッケージと関数を指定。
    パッケージ内でよほど何回も登場する関数や演算子を登録する。
    e.g., `@importFrom magrittr %>%` 。
    重複して何度も書いちゃっても大丈夫。

`@export`
:   `NAMESPACE` で `export()` する関数を指定。
    一般ユーザーからはこれが付いてる関数だけ見える。

`@param arg1 description...`
:   関数の引数。型や役割の説明を書く。

`@inherit package::function [fields...]`
:   別の関数から継承する。
    部分的に継承したければ関数名の後に指定。
:   params return title description details seealso sections references examples author source note
:   `@inheritParams` は `@param` の不足を補うショートカット。
:   `@noRd` と組み合わせて使えればいろいろ楽できそうだけどダメっぽい。

`@template template-name`
:   `man-roxygen/template-name.R` の内容を利用する。
    `@inheritParams` と違って、その関数で使用しないparamまで展開されちゃうのが欠点。

`@eval fun()`
:   関数を評価して出てきた文字列ベクタをroxygenコメントとして処理する。
    各要素の先頭は `#' ` ではなく@タグで。
    関数には引数も渡せるのでいろいろできる。

`@return description`
:   関数の返り値。

`@examples code...`
:   例となるコードを記述する。
    exportしないやつには書いてはいけない。
    `\dontrun{}` に入れるとチェックから除外される。
    単数形の `@example` は外部ファイルのパスを受け取る。

`@rdname basename`
:   `man/`に書き出すRdファイルの名前。
    複数の関数で同じものを指定すればヘルプをひとつにまとめられる。
    このときタイトルや `@param` などは共有される。
    似たような引数を持つ関数や、関連する関数をまとめるのによく使う。

`@include other-file.R`
:   指定したファイルを先に読み込む。
    メソッドの定義などで順序が重要になる場合に使う。

`@seealso ...`
:   `[mean()]`, `[ggplot2::diamonds]`, `<https://...>` のようにしてリンクを張れる。

`@family`
:   これを共通して持つ関数同士に `@seealso` が自動生成される。
    ただし `@rdname` とかで関数をまとめたりしてるとうまくいかないっぽい。

`@section Some Title:`
:   新しいセクションを作る。前には空行、行末にコロン、後には普通の文章。

`@docType`, `@name`
:   パッケージやデータを記述するのに必要だったが今では不要っぽい。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/1491910399/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=2ac4a1dea2487eff5b1034acbbad7ddb&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873117593/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=68d0ffbb02a546eca75795f147bd54e2&language=ja_JP" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117593&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=4873117593" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
