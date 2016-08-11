+++
title = 'stringr'
subtitle = "Rの文字列をまともな方法で処理する"
tags = ["r", "hadley"]
[menu.main]
  parent = "rstats"
  weight = -60
+++

-   <https://cran.r-project.org/web/packages/stringr/>
-   <https://github.com/hadley/stringr>
-   <http://www.rdocumentation.org/packages/stringr>

R標準の `base` パッケージが提供する関数とほとんど同じ機能のように見えるものもあるが、
統一的なインターフェイスに合理的な挙動なのでプログラムの中で使いやすい。

-   `factor` と `character` を同じように扱う
-   引数オブジェクトの各要素の名前や位置を保持する
    -   長さゼロのオブジェクトを引数として与えた場合には長さゼロの結果を返す
    -   引数オブジェクトに `NA` が含まれる場合はその部分の結果を `NA` とする
-   対象文字列が一貫して第一引数
-   [stringi](http://www.gagolewski.com/software/stringi/) を使って動くため高速
-   [ICU正規表現](http://userguide.icu-project.org/strings/regexp)

今や `stringr` は [stringi](http://www.gagolewski.com/software/stringi/) のラッパーだし、
どちらもほぼ同じインターフェイスなので、
もし前者に不足があれば後者を直接使えばよいが、
普通に使う分にはそんな場面には出くわさない。
むしろ、機能がある程度絞られているほうが取っ付き易いし、
`str_*` のほうが `stri_*` よりも1文字短いので、
基本的には `stringr` を使っとけばよい。

Rの中から `install.package('stringr')` でインストールし、
使う前に `library(stringr)` でパッケージを読み込む。

## Functions

### Basic Operation

`str_length(string)`
:   文字列の長さを数える。
    `base::nchar(x)` と相同だが、`NA` に対して `2` ではなく `NA` を返す。

`str_sub(string, start=1, end=-1)`
:   文字列を部分的に参照・変更する。
    `base::substr()` と相同だが、負数で末尾からの位置を指定できる。

`str_c(..., sep='', collapse=NULL)`
:   文字列を結合する。
    デフォルトの `sep` がスペースじゃないので `base::paste0()` に近い。

`str_split(string, pattern, n=Inf)`
:   文字列を分割する。
    `base::strsplit(x, split)` と相同だが、
    最大 `n` 個に分割するということを指定できる。
    空文字で帳尻合わせしてちょうど `n` 個にする `str_split_fixed()` もある。
    `string` と `pattern` の要素数が噛み合わないときにちゃんと警告が出る。

`str_dup(string, times)`
:   指定した回数だけ文字列を繰り返して結合。
    `str_dup('#', 79)` とかで結果出力に区切りを入れたり。

### Pattern Matching

`str_detect(string, pattern)`
:   マッチする箇所があるかどうか `logical` を返す。
    `base::grepl(pattern, x)` と相同。

`str_count(string, pattern)`
:   マッチする箇所の数を返す。

`str_locate(string, pattern)`
:   マッチする最初の箇所の `start`, `end` 位置を行列で返す。

`str_extract(string, pattern)`, `str_extract_all(string, pattern)`
:   マッチした部分文字列を取り出す。しなかった要素には `NA`。
    `base::grep(pattern, x, value=TRUE)` はマッチする要素のみ、元の形で返す。

`str_match(string, pattern)`, `str_match_all(string, pattern)`
:   マッチした部分文字列を取り出し、後方参照を含む行列を返す。
    `str_extract(string, pattern)` と同じ結果全体 `\0` が1列目で、
    カッコでマッチさせた `\1` 以降の結果が2列目以降に入る。

`str_replace(string, pattern, replacement)`
:   マッチしなかった部分をそのままに、マッチした部分を置換する。
    `base::sub(pattern, replacement, x)` と相同。
    `base::gsub()` のように全てのマッチを置換するには `str_replace_all()` 。

------------------------------------------------------------------------

引数`pattern`はデフォルトで[ICU正規表現](http://userguide.icu-project.org/strings/regexp)。
以下の関数を通して渡すことでそれを変更できる。

`stringr::fixed(string)`
:   そのままの文字としてマッチさせる

`stringr::ignore.case(string)`
:   大文字と小文字の違いを無視してマッチさせる

`stringr::perl(string)`
:   Perl拡張の正規表現として扱う

### Formatting

`str_trim(string, side="both")`
:   空白文字を除去する。
    Python でいうところの `str.strip()`。

`str_pad(string, width, side="left", pad=" ")`
:   余白を作る。
    幅を `width` に伸ばして `side` に寄せて空白を `pad` で埋める。

`str_wrap(string, width=80, indent=0, exdent=0)`
:   指定した幅で折り返す。
    `indent` は先頭行の左余白。
    `exdent` はそれ以外の行の左余白。
