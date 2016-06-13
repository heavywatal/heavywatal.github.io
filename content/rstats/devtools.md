+++
title = 'devtools'
subtitle = "Rパッケージ作成支援"
tags = ["r", "hadley"]
[menu.main]
  parent = "rstats"
  weight = -50
+++

<http://r-pkgs.had.co.nz/>

自分で書いた関数が多くなってきたら、まとめてパッケージを作るとよい。
少しだけ面倒だが、以下のようなメリットがある。

-   普通に変数や関数を定義するとワークスペースが名前でいっぱいになってしまうが、
    パッケージ内で定義されている変数や関数は `ls()` で出てこないのでスッキリ
-   既存のオブジェクトと名前が衝突するような場合でも、
    `mypackage::func1` のように名前空間を明示的に指定して呼び出せる
-   バイトコンパイルされる分だけ関数がちょっと速くなるかも

CRANに上げる程ではないとしても、
GitHubに公開しておけば誰でも使えるようになるので、
共同研究者と解析環境を共有したり、
ひとつの論文のワークフローを置いておいたり、いろいろ使い道はある。

## Rパッケージ

### 最低限の作成手順

1.  Rを起動して `devtools` をインストールする

    ```r
    install.packages('devtools')
    library(devtools)
    ```

2.  `devtools::create()` で骨組みを作る

    ```r
    setwd('~/tmp/')
    library(devtools)
    pkgname = 'namaespace'
    title_desc = 'Utility to create dummy R packages for namespace'
    github_repo = sprintf('https://github.com/heavywatal/%s', pkgname)
    devtools::create(pkgname, description=list(
        Package=pkgname,
        Title=title_desc,
        Description=paste0(title_desc, '.'),
        `Authors@R`="person('Watal M.', 'Iwasaki', email='user@example.com', role=c('aut', 'cre'))",
        License='MIT',
        Suggests='readr, tidyr, dplyr, ggplot2',
        Imports='devtools, stringr',
        URL=github_repo,
        BugReports=paste0(github_repo, '/issues'))
    ```

3.  `devtools::check(pkgname)` で様子を見てみる

    {{%div class="note"%}}
LaTeX 関連で怒られたら足りないパッケージを入れる:
`sudo tlmgr install inconsolata helvetic`
    {{%/div%}}

4.  GitHubに空のリポジトリを作る
    -   パッケージと同名でなくてもよい
    -   `README` や `LICENSE` は後から作る

5.  コミットしてプッシュ

    ```sh
    % git init
    % git add --all
    % git commit -m "first commit"
    % git remote add origin git@github.com:heavywatal/namaespace.git
    % git push -u origin master
    ```

6.  とりあえず誰でもインストール可能なパッケージができたはず

    ```r
    install_github('heavywatal/namaespace')
    ```

7.  あとは `R/` にソースコードを置いたり `README.md` を書いたり

### 構造

```sh
DESCRIPTION  # 一番大事
NAMESPACE    # 見せるオブジェクトを列挙
README.md    # 全体の説明を簡単に
R/           # Rソースコード
data/        # サンプルデータなど
man/         # ヘルプファイル.Rd
src/         # C++
tests/
vignettes/
```

### `DESCRIPTION`

-   どうでも良さそうなファイル名とは裏腹に、ちゃんと書かないと動かない
-   始めは `devtools::create()` で自動生成し、それから修正していく
-   `Imports` に列挙したパッケージは同時にインストールされる。
    しかし `library()` 時に同時に読み込まれるわけではない。
-   `Title` はピリオドを含まず、Descriptionはピリオドを含むようにする
-   ライセンスを別ファイルにする場合は `License: file LICENSE` と書く
-   `Authors@R` のとこは後でRで評価されるので変な形

### `NAMESPACE`

-   `roxygen2` (下記) がソースコードから自動生成するので**直接触らない**
-   ここで `export()` された関数だけがユーザーから見える
-   ここで `import()` された関数はパッケージ内でattachされた状態になるが、
    そうしないで毎回 `名前空間::` 越しにアクセスしたほうがよい。
    `magrittr::%>%` とか `pipeR::%>>%` のような演算子は仕方がないので
    `importFrom()` で個別指定する。

### Rソースコード

-   `R/` 以下に配置
-   ファイルの数や名前は何でもいいので、開発者が分かりやすいようにしとく。
    例えば、似たような関数をひとつのファイルにまとめ、同名の`@rdname`を指定し、
    `man/`と同じ構造にするのがちょうどいいのではなかろうか。
-   `library()` や `require()` を書かない。
    必要なパッケージは `DESCRIPTION` の `Imports` に書き、
    `名前空間::関数()` のようにして使う。
    どうしても名前空間を省略したいものだけ `@importFrom` を指定する。
-   `NAMESPACE` や `man/*.Rd` を自動生成してもらえるように
    roxygen形式でコメントを書く

    ```r
    #' A simple function to add 1
    #' @export
    #' @param x A numeric vector
    #' @return A numeric vector
    #' @examples
    #' increment(42)
    increment = function(x) {x + 1}
    ```

    -   <http://r-pkgs.had.co.nz/man.html>
    -   <http://kbroman.org/pkg_primer/pages/docs.html>
    -   <http://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html>

## `devtools`

<https://github.com/hadley/devtools>

骨組みを作るとこからCRANにデプロイするとこまでお世話してくれる。

### 関数

`create(path, description=getOption('devtools.desc'), check=FALSE, rstudio=TRUE)`
:   まっさらな状態から骨組みを作る。
    `path` が既に存在している場合は先に進めないので、
    やり直すときは `system('rm -rf path')` などとして一旦消す必要がある。

`load_all(pkg='.', reset=TRUE, recompile=FALSE, export_all=TRUE, quiet=FALSE)`
:   リロード

`check(pkg='.', document=TRUE, cleanup=TRUE, cran=TRUE, check_version=FALSE, ...)`
:   パッケージとしての整合性を確認

`document(pkg='.', ...)`
:   `roxygen2` を呼び出してソースコードから
    `NAMESPACE` や `man/*.Rd` を自動生成する

`test(pkg='.', filter=NULL, ...)`
:   `testthat` を呼び出して `test/` 以下のテストコードを実行する

`install(pkg='.', reload=TRUE, quick=FALSE, local=TRUE, ...)`
:   ローカルにあるディレクトリからインストール

`install_github(repo, username=NULL, ref='master', ...)`
:   GitHubリポジトリからインストール

### 設定

項目の説明を読む

```r
?devtools
```

例えば `.Rprofile` に

```r
options(devtools.desc.author='Watal M. Iwasaki <user@example.com> [aut, cre]')
options(devtools.desc.license='MIT')
```

## `roxygen2`

Rソースコードのコメントからヘルプを自動生成する。

<http://cran.r-project.org/web/packages/roxygen2/>

<https://github.com/klutometis/roxygen>

`roxygen2::roxygenise(package.dir='.', ..., clean=FALSE)`
を直接呼んでもよいが、
基本的には `devtools::document()` を使って間接的に利用する。

### 使い方

-   `#'` から始まる行がroxygenコメントとして扱われる。
-   タグは `@` で始まる。
    `@` そのものを入力したいときは重ねて `@@` とする。
-   1行目には必ずタイトルを書く。
    タグ `@title` は省略可能。
    とりあえずコピペして全部同じタイトルにしようとすると、
    関数ごとにユニークなタイトルをつけろと怒られる。
-   2段落目は `@description` 扱いされる。
    省略するとタイトルが流用される。
-   空行だけでは切れ目として扱ってくれないので、
    データやパッケージのドキュメントを書く場合は `NULL` を入れる。

### タグ

`@import pkg1, pkg2, ...`
:   `NAMESPACE` で `import()` するパッケージを指定。
    名前の衝突が怖いので基本的に使わない。

`@importFrom pkg func`
:   `NAMESPACE` で `importFrom()` するパッケージを指定。
    e.g., `@importFrom pipeR %>>%`

`@export`
:   `NAMESPACE` で `export()` するパッケージを指定。

`@param arg1 description...`
:   関数の引数。型や役割の説明を書く。

`@return description`
:   関数の返り値。

`@examples code...`
:   例となるコードを記述する。
    exportしないやつには書いてはいけない。
    単数形の `@example` は外部ファイルのパスを受け取る

`@rdname basename`
:   `man/`に書き出すRdファイルの名前。
    複数の関数で同じものを指定すればひとつのヘルプにまとめられる。
    このとき`@param`などは共有されるので、
    同じ名前のものはどこかで1度だけ記述する。
    逆に言えば、中身が違うものに同じ名前をつけてはいけないし、
    引数や機能がほとんど重ならない関数をまとめると分かりにくくなる。

`@docType`
:   関数やクラスには不要だが `data` か `package` の場合はここで指定。

`@name name`
:   パッケージやデータはこのタグで明示的に設定する必要がある。
