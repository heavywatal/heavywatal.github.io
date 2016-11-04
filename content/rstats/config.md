+++
title = 'R環境設定'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -95
+++

https://cran.r-project.org/manuals.html

## インストール

https://cran.r-project.org/doc/manuals/R-admin.html

### Mac

- [R本体](https://cran.rstudio.com/bin/macosx/)
- [Rstudio (任意)](https://www.rstudio.com/products/rstudio/download/)

### Ubuntu

ターミナルからリポジトリを追加して

```sh
sudo sh -c 'echo "deb https://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -cs)/" > /etc/apt/sources.list.d/cran-mirror.list'
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get --quiet update
sudo apt-get install r-base
```

## 起動オプション

ターミナルから `R --help` を見るか
https://cran.r-project.org/doc/manuals/R-intro.html#Invoking-R

例えば以下のようなエイリアスを `.zshrc` で設定する。

```sh
alias R='R --quiet --no-save --no-restore-data'
```

## 環境変数

https://cran.r-project.org/doc/manuals/R-admin.html#Environment-variable-index

https://stat.ethz.ch/R-manual/R-patched/library/base/html/EnvVar.html

```r
Sys.getenv()
```


## .Renviron

<https://stat.ethz.ch/R-manual/R-patched/library/base/html/Startup.html>

R起動時に読み込まれ、環境変数を設定するファイル。
Rスクリプトではなく、シェルスクリプトっぽい代入式で書く。
例 (<https://github.com/heavywatal/rwtl/blob/master/.R/.Renviron>):

    R_USER=${HOME}/.R
    R_LIBS_USER=${R_USER}/library
    R_ENVIRON_USER=${R_USER}/.Renviron
    R_PROFILE_USER=${R_USER}/.Rprofile
    R_HISTFILE=${R_USER}/.Rhistory
    R_HISTSIZE=65535
    LANG=en_US.UTF-8
    LANGUAGE=en_US.UTF-8
    LC_ALL=en_US.UTF-8

読み込まれる順序はだいたい以下のとおりなので、ホームディレクトリにシムリンクを張っておけば読み込まれる。

1.  `$R_ENVIRON`
2.  `$R_HOME/etc/Renviron.site`
3.  `$R_ENVIRON_USER`
4.  `./.Renviron`
5.  `~/.Renviron`

```sh
% cd
% ln -s .R/.Renviron
```

## .Rprofile

https://cran.r-project.org/doc/manuals/R-intro.html#Customizing-the-environment

<https://stat.ethz.ch/R-manual/R-patched/library/base/html/options.html>

R起動時に読み込まれるファイル。
中身はRスクリプトなので、パッケージの読み込みや関数の定義など、Rでできることは何でもできるはず。
例: <https://github.com/heavywatal/rwtl/blob/master/.R/.Rprofile>

読み込まれる順序はだいたい以下のとおり。
`.Renviron` のほうが先に読み込まれるので、
上記のように `R_PROFILE_USER` を定義しておいて、そこに置いとけば読み込まれる。

1.  `$R_PROFILE`
2.  `$R_HOME/etc/Rprofile.site`
3.  `$R_PROFILE_USER`
4.  `./.Rprofile`
5.  `~/.Rprofile`

`.First()` と `.Last()` はそれぞれ起動時と終了時に実行される関数。
これらが原因で `install.packages()` がエラーを引き起こすこともあるので注意。

### `options()`

`?options` で項目の一覧を見られる。

スクリプトの結果が設定依存で変わっては困るので、
そういう本質的なものはいじらずに、
表示関連の項目だけ変えるに留めたほうがよい。
例えば `options(stringsAsFactors=FALSE)` とやってしまいたくなるが、
それを設定した環境とそうでない環境で `read.table()` の結果が変わるのは危険。
関数の引数として毎回指定するか、
[readr]({{< relref "readr.md" >}})の関数を使うべし。

`warn=1`
: 警告レベルの設定。
  デフォルト(`warn=0`)では、警告があっても計算は滞り無く進行し、
  最後に "There were 50 or more warnings (use warnings() to see the first 50)"
  などと軽く表示されるだけなので、
  見落としたりして後々大変なバグ取り作業に発展する恐れがある。
  警告が発生するごとに警告文を表示する(`warn=1`)か、
  エラー扱いにして計算をストップするようにしておく(`warn=2`)ことでそれを回避すべし。
  ちなみに負数だと警告無視。

`showWarnCalls=TRUE`, `showErrorCalls=TRUE`
: 警告やエラーの出処をたどって表示する。
  Pythonほどわかりやすくないが、ちょっとはマシになる。


## ライブラリの管理

https://cran.r-project.org/doc/manuals/R-admin.html#Add_002don-packages
