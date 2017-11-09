+++
title = 'stringr'
subtitle = "Rの文字列をまともな方法で処理する"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -60
+++

<a href="http://stringr.tidyverse.org/">
<img src="http://stringr.tidyverse.org/logo.png" align="right">
</a>

R標準の`base`パッケージが提供する関数でも文字列処理は可能だが、
`stringr`のほうが統一的なインターフェイスに合理的な挙動で使いやすい。

-   `factor` と `character` を同じように扱う
-   引数オブジェクトの各要素の名前や位置を保持する
    -   長さゼロのオブジェクトを引数として与えた場合には長さゼロの結果を返す
    -   引数オブジェクトに `NA` が含まれる場合はその部分の結果を `NA` とする
-   対象文字列が一貫して第一引数で、パターンが二番目
-   何をやる関数なのか名前から分かりやすい\
    (標準が覚えにくすぎ: `grep`, `grepl`, `regexpr`, `gregexpr`, `regexec`)
-   [ICU4C](http://site.icu-project.org/)
    (via [stringi](http://www.gagolewski.com/software/stringi/)) を使って動くため高速
-   [ICU正規表現](http://userguide.icu-project.org/strings/regexp) の仕様が明確

今や `stringr` は [stringi](http://www.gagolewski.com/software/stringi/) のラッパーだし、
どちらもほぼ同じインターフェイスなので、
もし前者に不足があれば後者を直接使えばよいが、
普通に使う分にはそんな場面には出くわさない。
むしろ、機能がある程度絞られているほうが取っ付き易いし、
`str_*` のほうが `stri_*` よりも1文字短いので、
基本的には `stringr` を使っとけばよい。

[tidyverse](https://github.com/tidyverse/tidyverse) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <http://r4ds.had.co.nz/strings.html>
-   <https://cran.r-project.org/web/packages/stringr/>
-   <https://github.com/tidyverse/stringr>
-   <http://www.rdocumentation.org/packages/stringr>

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

`str_split(string, pattern, n=Inf, simplify=FALSE)`
:   文字列を分割してlistを返す。
    `base::strsplit(x, split)` と相同だが、
    最大 `n` 個に分割するということを指定できる。
    空文字で帳尻合わせしてちょうど `n` 個にするショートカットが
    `str_split_fixed(string, pattern, n)` 。
    `string` と `pattern` の要素数が噛み合わないときにちゃんと警告が出る。
    `simplify=TRUE` とするとmatrixで返す。

`str_dup(string, times)`
:   指定した回数だけ文字列を繰り返して結合。
    `str_dup('#', 79)` とかで結果出力に区切りを入れたり。

### Pattern Matching

`str_count(string, pattern)`
:   マッチする箇所の数を返す。

`str_detect(string, pattern)`
:   マッチするかどうか `logical` を返す。
    `base::grepl(pattern, x)` と相同。

`str_extract(string, pattern)`, `str_extract_all(string, pattern)`
:   マッチした部分文字列を取り出す。しなかった要素には `NA`。
:   数値＋単位のような文字列から数値部分だけを抜き出すには
    [readr::parse_number()]({{< relref "readr.md#parse" >}})
    が便利。

`str_subset(string, pattern)`
:   `x[str_detect(x, pattern)]` のショートカット。
    マッチする要素だけ元の形で返すので
    `str_extract()` より `base::grep(pattern, x, value=TRUE)` に近い。

`str_which()`
:   マッチする要素のインデックスを整数で返す
    `which(str_detect(x, pattern))` のショートカット。
    `base::grep(pattern, x)` と相同。

`str_locate(string, pattern)`
:   マッチする最初の箇所の `start`, `end` 位置を行列で返す。

`str_match(string, pattern)`, `str_match_all(string, pattern)`
:   マッチした部分文字列を取り出し、後方参照を含む行列を返す。
    `str_extract(string, pattern)` と同じ結果全体 `\0` が1列目で、
    カッコでマッチさせた `\1` 以降の結果が2列目以降に入る。

`str_replace(string, pattern, replacement)`
:   マッチしなかった部分をそのままに、マッチした部分を置換する。
    `base::sub(pattern, replacement, x)` と相同。
    `base::gsub()` のように全てのマッチを置換するには `str_replace_all()` 。

------------------------------------------------------------------------

上記関数の`pattern`引数は普通に文字列を渡すと正規表現として解釈してくれるが、
下記の関数を通して渡すことでその挙動を変更することができる。

`stringr::regex(pattern, ignore_case=FALSE, multiline=FALSE, comments=FALSE, dotall=FALSE, ...)`
:   デフォルトの[ICU正規表現](http://userguide.icu-project.org/strings/regexp)。
    複数行ファイルに対するマッチではこの関数を通して挙動をいじることになる。

`stringr::fixed(pattern)`
:   正規表現ではなくそのままの文字としてマッチさせる

`stringr::boundary(type='character', skip_word_none=NA, ...)`
:   境界に対するマッチ。
    `type`の選択肢は `character`, `line_break`, `sentence`, `word`.

`stringr::coll(pattern, ignore_case=FALSE, locale=NULL, ...)`
:   よくわからないけど非ascii対策？


### Formatting

`str_to_upper()`, `str_to_lower()`, `str_to_title()`
:   大文字・小文字の変換

`str_interp(string, env=parent.frame())`
:   `sprintf()` と相同。
    文字列の中の `$[format]{expr}` がR表現として評価される。
    `[format]`部分は`sprintf()`と同じ形式で、省略可。
    `env` はlistやdata.frameでもよい。
:   e.g., `stringr::str_interp('Mean sepal width is $[.3f]{mean(Sepal.Width)}.', iris)`
:   新しい [`library(glue)`](http://glue.tidyverse.org/) も良さそう。

`str_pad(string, width, side="left", pad=" ")`
:   余白を作る。
    幅を `width` に伸ばして `side` に寄せて空白を `pad` で埋める。

`str_trim(string, side="both")`
:   空白文字を除去する。
    Python でいうところの `str.strip()`。

`str_trunc(string, width, side=c('right', 'left', 'center'), ellipsis='...')`
:   一定の長さを超えたら捨てて `...` にする。

`str_wrap(string, width=80, indent=0, exdent=0)`
:   指定した幅で折り返す。
    `indent` は先頭行の左余白。
    `exdent` はそれ以外の行の左余白。

文字列と数値の型変換はstringrの管轄外なので、標準の
`as.character()` や `as.double()` などを使うか、
[`readr::parse_*()`系の関数]({{< relref "readr.md#parse" >}})
を使う。


## 関連書籍

<a href="https://www.amazon.co.jp/Data-Science-Transform-Visualize-Model/dp/1491910399/ref=as_li_ss_il?s=books&ie=UTF8&qid=1508340700&sr=1-3&linkCode=li3&tag=heavywatal-22&linkId=6a53371cc80e2d6d7fc50fae5d8b862d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1491910399&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=1491910399" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/R%E3%81%A7%E3%81%AF%E3%81%98%E3%82%81%E3%82%8B%E3%83%87%E3%83%BC%E3%82%BF%E3%82%B5%E3%82%A4%E3%82%A8%E3%83%B3%E3%82%B9-Hadley-Wickham/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
