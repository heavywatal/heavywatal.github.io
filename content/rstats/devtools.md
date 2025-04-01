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

<a href="https://usethis.r-lib.org/">
<img src="/_img/hex-stickers/usethis.webp" align="right" width="120" height="139">
</a>

- <https://r-pkgs.org/>
- [How to develop good R packages (for open science)](https://masalmon.eu/2017/12/11/goodrpackages/)

### 最低限の作成手順

1.  開発支援パッケージをインストールする:
    `install.packages(c("devtools", "usethis"))`

1.  [usethis](https://usethis.r-lib.org/) の関数をいくつか使って骨組みを作る:

    ```r
    usethis::create_package("hello")
    usethis::use_mit_license()
    usethis::use_package_doc()
    ```

1.  `devtools::check()` で様子を見てみる。

1.  ローカルgitリポジトリを作って最初のコミットをする:
    `usethis::use_git()`

1.  GitHubに同名のリポジトリを作ってプッシュ:
    `usethis::use_github()`

    （リポジトリの名前をパッケージ名と違うものにしたい場合などは手動で）

1.  とりあえず誰でもインストール可能なパッケージができたはず:
    ```r
    remotes::install_github("heavywatal/rhello")
    ```


### ソース構造

<https://r-pkgs.org/structure.html>

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


### `DESCRIPTION`

<https://r-pkgs.org/description.html>

-   どうでも良さそうなファイル名とは裏腹に、ちゃんと書かないと動かない。
-   始めはusethisとかに生成してもらい、他のパッケージを参考にしつつ修正していく。
-   `Imports` に列挙したものは依存パッケージとして同時にインストールされる。
    しかし `library()` 時に同時に読み込まれるわけではない。
-   `Title` はピリオドを含まない一文でタイトルケース。
    `Description` はピリオドを含む一段落。
-   ライセンスを別ファイルにする場合は `License: file LICENSE` と書く
-   `Authors@R` のとこは後でRで評価されるので変な形。

### `NAMESPACE`

<https://r-pkgs.org/namespace.html>

-   [後述のroxygen2](#roxygen2)がソースコードから自動生成するので**直接触らない**。
-   ここで `export()` された関数だけがユーザーから見える。
-   外部パッケージから `importFrom(package, function)`
    された関数はattachされてパッケージ内で利用可能になる。
    `import(package)` で全ての関数をまとめて処理できるけど名前の衝突が怖い。

### `R/` ソースコード

<https://r-pkgs.org/code.html>

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

### `src/` C++ソースコード

<https://r-pkgs.org/src.html>

[Rcppページの"Rパッケージで使う"セクション]({{< relref "rcpp.md#Rパッケージで使う" >}})を参照


### `vignettes/`

<https://r-pkgs.org/vignettes.html>

<a href="https://pkgdown.r-lib.org/">
<img src="/_img/hex-stickers/pkgdown.webp" align="right" width="120" height="139">
</a>

個々の関数の使用例はRソースファイルの `@examples` に書くとして、
複数の関数を組み合わせてどう使うかを説明するのが`vignettes/`の役割。
R Markdown形式で書いてHTMLやPDFで表示できるので表現力豊か。

`usethis::use_vignette("hello")` で雛形を作ってもらうのが楽。

`pandoc` と `pandoc-citeproc` が必要なので
[Homebrew]({{< relref "homebrew.md" >}}) とかでインストールしておく。

`check()`がデフォルトで`vignette=TRUE`かつ処理がやや重いので、
毎回その重さを受け入れるか、わざわざ`FALSE`指定するかというのは悩みどころ。

[pkgdown](https://pkgdown.r-lib.org)でウェブサイトを構築すると、
ここに置いてある文書は
[Articles](https://pkgdown.r-lib.org/articles/pkgdown.html#articles)
という位置づけで出力される。

ただしパッケージと同名で `vignettes/<package-name>.Rmd` というファイルを作ると、
Articles一覧の中ではなくReferenceの隣に "Get started" としてリンクされる。
ここでパッケージの使い方を説明することが期待されている。
`index.md`/`README.md` から生成されるホームページも似たような役割じゃないかと思うが、
そちらはパッケージを使うかどうか判断するための情報を提供するところとして書き、
「実際に使いたくなったら "Get started" を見よ」
というリンクを末尾に置いておくのが良さそう。
[pkgdown#2372](https://github.com/r-lib/pkgdown/issues/2372)


### `tests/`

<https://r-pkgs.org/tests.html>

<a href="https://testthat.r-lib.org/">
<img src="/_img/hex-stickers/testthat.webp" align="right" width="120" height="139">
</a>

<a href="https://covr.r-lib.org/">
<img src="/_img/hex-stickers/covr.webp" align="right" width="120" height="139">
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


### その他

<https://r-pkgs.org/misc.html>

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

`inst/`
: ここに入ってるものはインストール先でトップディレクトリに移される謎仕様。
: 論文で引用されることを想定している場合は `inst/CITATION` を作る。
  `citation("ggplot2")` のように参照できるようになる。

`tools/`
: たぶん何を入れても一切チェックされずインストールもされない自由なディレクトリ。
  `tests/` や `vignettes/` に入れる前のガラクタコードを一時的に置いとくとか。

`.Rbuildignore`
: Rパッケージらしからぬ変なファイルやディレクトリがあると怒られるので、
  そういうやつをこのファイルに列挙して無視してもらう。

`README.md`
: GitHubでいい感じに見えるようにMarkdown形式でパッケージの概要を書く。
  Rコードを含む `README.Rmd` から `devtools::build_readme()` で生成することも可能。
: pkgdownでもホームページの材料として使われるが、`index.md` を置けばそちらが優先される。
  開発者向けと一般向けを使い分けるとか。
  `index.md` をルートに置きたくない場合は `pkgdown/index.md` でもいい。


## `devtools`

<a href="https://devtools.r-lib.org/">
<img src="/_img/hex-stickers/devtools.webp" align="right" width="120" height="139">
</a>

骨組みを作るとこからCRANにデプロイするとこまでお世話してくれる。
いろんな専門パッケージの集合体。

### 主な関数

`document(pkg = ".", roclets = NULL, quiet = FALSE)`
:   [`roxygen2`](#roxygen2) を呼び出してソースコードから
    `NAMESPACE` や `man/*.Rd` を自動生成する。

`check(pkg = ".", document = NA, ...)`
:   パッケージとしての整合性を確認。
    ついでに `document()` は実行できるけど `spell_check()` はできないので手動で。

`test(pkg = ".", filter = NULL, ...)`
:   `testthat` を呼び出して [`tests/`](#tests) 以下のテストコードを実行する

`build(pkg = ".", path = NULL, ...)`
:   `R CMD install` でインストール可能な tar ball (bundle package) を作る。
    `src/` のコードをコンパイルするという意味ではない。

`install(pkg = ".", reload = TRUE, quick = FALSE, build = !quick, ...)`
:   ローカルにあるソースからインストール。
    `build = TRUE` のとき(デフォルト)、わざわざ bundle package を
    `tempdir()` に作ってからそいつでインストールする。

`install_github(repo, ref = "HEAD", subdir = NULL, ...)`
:   GitHubリポジトリからインストール。

`unload(pkg = ".", quiet = FALSE)`
:   `detach("package:XXX")` とか `unloadNamespace(XXX)`
    よりもちゃんとまっさらにパッケージを外す。


`load_all(pkg = ".", reset = TRUE, recompile = FALSE, export_all = TRUE, ...)`
:   `install()` せずファイルから直接 `library()` する。
    ロード済みでもまず `unload()` が呼ばれるので安心。
:   `load_all` 状態のパッケージに対して `system.file()` を呼び出すと、
    [`pkgload::system.file()`](https://pkgload.r-lib.org/reference/system.file.html)
    が間に割り込み、ソースのトップと `inst/` を起点にして探索してくれる。
    `configure` で生成するファイルを見つけてもらうには
    `$R_PACKAGE_DIR` に直接送り込まず一旦 `inst/` などに置く必要がある。

`clean_dll(pkg = ".")`
:   `src/` 以下に生成される `.o`, `.so` を消す。
    普段は触る必要ないが、たまにこれが必要な不具合に出くわす。


## `roxygen2`

<a href="https://github.com/klutometis/roxygen">
<img src="/_img/hex-stickers/roxygen2.webp" align="right" width="120" height="139">
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

<https://roxygen2.r-lib.org/articles/rd.html>

```r
#' Title of the simple function to add 1 (without explicit @@title tag)
#'
#' The second paragraph is recognized as the description.
#' Explicit @@description is unnecessary unless you want to include
#' an empty line to express multiple paragraphs or bullet lists.
#'
#' The third and subsequent paragraphs are recognized as the details.
#' @param x A numeric vector.
#' @returns A numeric vector with 1 added to each element.
#' @export
#' @examples
#' increment(42)
increment = function(x) {x + 1}
```

-   `#' ` から始まる行がroxygenコメントとして扱われる。
-   タグは `@` で始まる。
    `@` そのものを入力したいときは重ねて `@@` とする。
-   1行目にタイトル。1行あけて2段落目に説明文を書く。
    明示的な `@title`, `@description`, `@details` タグは大概不要。
    箇条書きや複数段落を表現したい場合にだけ適宜使う。
    タイトルをコピペして全部同じにすると怒られる。
    2段落目を省略するとタイトルが流用される。
-   空行だけでは切れ目として扱われないので `NULL` などを置いたりする。
-   [Rd形式の代わりにMarkdown形式で記述できる](https://roxygen2.r-lib.org/articles/markdown.html)。
    `usethis::create_package()` がデフォルトで
    `Roxygen: list(markdown = TRUE)` を `DESCRIPTION` に書いてくれる。
    その全体設定をせず使いたいブロックにいちいち `@md` を書く個別設定も可能。
-   `"_PACKAGE"` という文字列の上に書かれたブロックは、
    パッケージそのものに対応するヘルプ `*-package.Rd` になる。
    `@useDynLib` など全体に関わる設定はここでやるのが良いのでは。

    ```r
    #' Example package to say hello
    #' @useDynLib hello, .registration = TRUE
    #' @importFrom rlang .data :=
    #' @keywords internal
    "_PACKAGE"
    ```

-   dplyrなどで列名を直に指定すると
    `undefined global variables`
    という警告が出るので
    `dplyr::filter(diamonds, .data$cut == "Ideal")`
    のように pronoun を使って抑える。
    そのためにどこかに `#' @importFrom rlang .data` を書いておく。
    ただし `group_by()` で多数のグループを処理するときなど、
    `.data$` の遅さが気になる場合は `!!as.name("carat")` のようにする。


### タグ

使用可能なタグ一覧は[準備中？](https://github.com/klutometis/roxygen/issues/792)

`@import pkg1, pkg2, ...`
:   `NAMESPACE` で `import()` するパッケージを指定。
    名前の衝突が怖いので基本的に使わない。

`@importFrom pkg func`
:   `NAMESPACE` で `importFrom()` するパッケージと関数を指定。
    パッケージ内でよほど何回も登場する関数や
    `名前空間::` の形で使いにくい演算子を登録する。
    e.g., `@importFrom rlang .data :=` 。
    重複して書いても大丈夫。

`@export`
:   `NAMESPACE` で `export()` する関数を指定。
    一般ユーザーからはこれが付いてる関数だけ見える。
    基本的にはこれをつける関数にだけドキュメントを書く。

`@param arg1 description...`
:   関数の引数。型や役割の説明を文として書く。つまり大文字で始まりピリオドで終わる。

`@returns ...`
:   返り値の型などを簡潔に、文として書く。
    "Value" というセクションに出力される。
    `@description` で済むような場合でもCRANに求められるらしい。
    `@return` も同じだが tidyverse 界隈では `@returns` のほうが好まれている。

`@examples code...`
:   例となるコードを記述する。
    最後に空行を入れるとそれもしっかり含まれてしまう。
:   親環境に影響を与えないように注意する。
    ファイルの書き出しや `options()` など副作用を伴う例を書く場合は手動で現状復帰する。
    `tempdir()` は使えるが `on.exit()` と [`withr`](https://withr.r-lib.org/)
    は使えないらしい。
:   使い方は見せるけど実行しない、結果を見せない場合は `\dontrun{}` に入れる。
    `@examplesIf` で条件付きにもできる。
:   単数形の `@example` は外部ファイルのパスを受け取る。

`@rdname basename`
:   `man/`に書き出すRdファイルの名前。
    複数の関数で同じものを指定すればヘルプをひとつにまとめられる。
    このときタイトルや `@param` などは共有される。
    似たような引数を持つ関数や、関連する関数をまとめるのによく使う。

`@include other-file.R`
:   指定したファイルを先に読み込むように `DESCRIPTION` の `Collate` を自動生成する。
    ひとつでも使うと全ファイルを列挙する形になってしまうので、
    メソッドの定義などで順序が重要になる場合にのみ使う。

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

`@seealso ...`
:   `[mean()]`, `[ggplot2::diamonds]`, `<https://...>` のようにしてリンクを張れる。

`@family`
:   これを共通して持つ関数同士に `@seealso` が自動生成される。
    ただし `@rdname` とかで関数をまとめたりしてるとうまくいかないっぽい。

`@section Some Title:`
:   `@description` や `@details` にセクションを作るための古いタグ。
    今はMarkdownで簡単に書ける: `# Section`, `## Subsection`

`@docType`, `@name`
:   パッケージやデータを記述するのに必要だったが今では不要っぽい。
