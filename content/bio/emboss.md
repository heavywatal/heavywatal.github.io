+++
title = 'EMBOSS'
[menu.main]
  parent = "bio"
+++

The European Molecular Biology Open Software Suite

<http://emboss.sourceforge.net/>

## Installation

1.  <ftp://emboss.open-bio.org/pub/EMBOSS/>
    からソースコードをダウンロードしてきて展開し、
    そのディレクトリに `cd` する:

        % wget -O- ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-latest.tar.gz | tar xz
        % cd EMBOSS-6.6.0/

2.  コンパイルしてインストール:

        % ./configure --disable-debug --disable-dependency-tracking --enable-64 --with-thread --without-x --prefix=/usr/local/emboss
        % make && make check
        % sudo make install

3.  パスを通す。
    例えば `.zshenv` や `.bashrc` に以下のように記述する:

        export PATH=${PATH}:/usr/local/emboss/bin

{{%div class="note"%}}
Mac なら [Homebrew]({{< relref "mac/homebrew.md" >}}) でも入れられる

```sh
% brew tap homebrew/science
% brew install emboss
```
{{%/div%}}

## EMBASSY

<http://emboss.sourceforge.net/embassy/>

### PHYLIP

PHYLogeny Inference Package

<http://evolution.genetics.washington.edu/phylip.html>

本家PHYLIPに含まれる protdist や dnadist
はコマンドライン引数が取れないので、そのままでは大変使いにくい。
そのへんを改良したPHYLIPNEWというのがEMBASSYパッケージとして提供されているのでこっちを使う。
しかもこっちなら [BioPython]({{< relref "python/biopython.md" >}}) から使うこともできるっぽい。
元のPHYLIPと区別するためか、プログラム名の先頭にすべて `f` が付く。

1.  EMBOSSをインストール
2.  <ftp://emboss.open-bio.org/pub/EMBOSS/> からソースコードをダウンロードして展開し、
    そこに `cd` する:

        % wget -O- ftp://emboss.open-bio.org/pub/EMBOSS/PHYLIPNEW-3.69.tar.gz | tar xz
        % cd PHYLIPNEW-3.69/

3.  EMBOSSと同じ `prefix` を指定して `configure` 、コンパイル、インストール:

        % ./configure --prefix=/usr/local/emboss --without-x
        % make
        % sudo make install

4.  EMBOSSインストール時にパスを通してあればそのまま使えるはず:

        % fdnadist --help

{{%div class="note"%}}
`make` すると
`/usr/local/emboss/include/ajreg.h:19: fatal error: pcre_config.h: No such file or directory`
のように怒られる場合がある。
どっちかっつーとEMBOSSの問題だと思うんだけど、
そのファイルは `/usr/local/emboss/include/` ではなく
その下の `epcre/` にインストールされている。
`configure` かコンパイルの時に `-I` オプションか何かを指定してやるのが
正しい解決なのかもしれないけど、簡単にシムリンクを張ってやり過ごす

```sh
% cd /usr/local/emboss/include/
% sudo ln -s epcre/* ./
```
{{%/div%}}

## GEMBASSY

<http://www.g-language.org/gembassy/>
