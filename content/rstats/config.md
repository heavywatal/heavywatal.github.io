+++
title = 'R環境設定'
tags = ["r"]
[menu.main]
  parent = "rstats"
  weight = -95
+++

## インストール

### Mac

<http://cran.rstudio.com/bin/macosx/>

<http://www.rstudio.com/products/rstudio/download/>

### Ubuntu

ターミナルからリポジトリを追加して

```sh
sudo sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -cs)/" > /etc/apt/sources.list.d/cran-mirror.list'
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get --quiet update
sudo apt-get install r-base
```

## 起動オプション

<http://cran.r-project.org/doc/manuals/r-release/R-intro.html#Invoking-R>

あるいはターミナルから

```sh
% R --help
```

## 環境変数

<http://stat.ethz.ch/R-manual/R-patched/library/base/html/EnvVar.html>

## .Renviron

<http://stat.ethz.ch/R-manual/R-patched/library/base/html/Startup.html>

R起動時に読み込まれ、環境変数を設定するファイル。
Rスクリプトではなく、シェルスクリプトっぽい代入式で書く:

    R_USER=${HOME}/rstats
    R_LIBS_USER=${R_USER}/library
    R_ENVIRON_USER=${R_USER}/.Renviron
    R_PROFILE_USER=${R_USER}/.Rprofile
    R_HISTFILE=${R_USER}/.Rhistory
    R_HISTSIZE=1024
    LANG=en_US.UTF-8
    LANGUAGE=en_US.UTF-8
    LC_ALL=en_US.UTF-8

読み込まれる順序はだいたい以下のとおりなので、ホームディレクトリにシムリンクを張っておけば読み込まれる。

1.  `$R_ENVIRON`
2.  `$R_HOME/etc/Renviron.site`
3.  `$R_ENVIRON_USER`
4.  `./.Renviron`
5.  `$HOME/.Renviron`

```sh
% cd
% ln -s local/lib/R/.Renviron
```

## .Rprofile

<http://cran.r-project.org/doc/manuals/r-release/R-intro.html#Customizing-the-environment>

<http://stat.ethz.ch/R-manual/R-patched/library/base/html/options.html>

R起動時に読み込まれるファイル。
中身はRスクリプトなので、パッケージの読み込みや関数の定義など、Rでできることは何でもできるはず。
読み込まれる順序はだいたい以下のとおり。
`.Renviron` のほうが先に読み込まれるので、
上記のように `R_PROFILE_USER` を定義しておいて、そこに置いとけば読み込まれる。

1.  `$R_PROFILE`
2.  `$R_HOME/etc/Rprofile.site`
3.  `$R_PROFILE_USER`
4.  `./.Rprofile`
5.  `$HOME/.Rprofile`

`.First()` と `.Last()` はそれぞれ起動時と終了時に実行される関数。
これらが原因で `install.packages()` がエラーを引き起こすこともあるので注意。

## ライブラリの管理

<http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Add_002don-packages>
