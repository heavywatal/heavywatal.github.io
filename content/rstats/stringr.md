+++
title = 'stringr'
subtitle = "Rの文字列をまともな方法で処理する"
tags = ["r", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -60
+++



<a href="https://stringr.tidyverse.org/">
<img src="/_img/hex-stickers/stringr.webp" align="right" width="120" height="139">
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
-   [ICU4C](https://icu.unicode.org/)
    (via [stringi](https://stringi.gagolewski.com/)) を使って動くため高速
-   [ICU正規表現](https://unicode-org.github.io/icu/userguide/strings/regexp.html) の仕様が明確

今や `stringr` は [stringi](https://stringi.gagolewski.com/) のラッパーだし、
どちらもほぼ同じインターフェイスなので、
もし前者に不足があれば後者を直接使えばよいが、
普通に使う分にはそんな場面には出くわさない。
むしろ、機能がある程度絞られているほうが取っ付き易いし、
`str_*` のほうが `stri_*` よりも1文字短いので、
基本的には `stringr` を使っとけばよい。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages("tidyverse")` で一括インストール、
`library(tidyverse)` で一括ロード。

-   <https://r4ds.hadley.nz/strings.html>
-   [base R関数との対応表](https://stringr.tidyverse.org/articles/from-base.html)

## Functions

### Basic Operation

`str_length(string)`
:   文字列の長さを数える。
    `base::nchar(x)` と相同。
    
    ```r
    str_length(c("NA", NA))
    ```
    
    ```
    [1]  2 NA
    ```

`str_sub(string, start = 1, end = -1)`
:   文字列を部分的に参照・変更する。
    `base::substr()` と相同だが、負数で末尾からの位置を指定できる。
    ただしRのインデックスは1始まりで終端も含むの。
    `str_sub<-` が定義されているので置換にも使える。
    
    ```r
    str_sub("supercalifragilisticexpialidocious", 10, -15)
    ```
    
    ```
    [1] "fragilistic"
    ```

`str_flatten(string, collapse = "")`
:   文字列vectorを1つの文字列に結合する。
    `base::paste0(string, collapse = "")` と同等だが `NA` を扱える。
    
    ```r
    str_flatten(c("Dragon", NA, "Force"))
    ```
    
    ```
    [1] NA
    ```
    
    ```r
    str_flatten(c("Dragon", NA, "Force"), na.rm = TRUE)
    ```
    
    ```
    [1] "DragonForce"
    ```
    
    ```r
    paste0(c("Dragon", NA, "Force"), collapse = "")
    ```
    
    ```
    [1] "DragonNAForce"
    ```

`str_c(..., sep = "", collapse = NULL)`
:   複数の引数で与えた文字列を結合する。
    デフォルトの `sep` がスペースじゃないので `base::paste0()` に近い。
    
    ```r
    str_c(c("Dragon", "Hammer"), c(NA, ""), c("Force", "Fall"))
    ```
    
    ```
    [1] NA           "HammerFall"
    ```
    
    ```r
    paste0(c("Dragon", "Hammer"), c(NA, ""), c("Force", "Fall"))
    ```
    
    ```
    [1] "DragonNAForce" "HammerFall"   
    ```

`str_split(string, pattern, n = Inf, simplify = FALSE)`
:   文字列を分割してlistを返す `base::strsplit(x, split)` の改良版。
    `string` と `pattern` の要素数が噛み合わないときにちゃんと警告が出る。
    最大 `n` 個に分割するということを指定できる。
    `simplify = TRUE` とするとmatrixで返す。
    
    ```r
    str_split(c("DragonForce", "HammerFall"), "(?<=[a-z])(?=[A-Z])")
    ```
    
    ```
    [[1]]
    [1] "Dragon" "Force" 
    
    [[2]]
    [1] "Hammer" "Fall"  
    ```
    
    ```r
    str_split(c("DragonForce", "HammerFall"), "(?<=[a-z])(?=[A-Z])", simplify = TRUE)
    ```
    
    ```
         [,1]     [,2]   
    [1,] "Dragon" "Force"
    [2,] "Hammer" "Fall" 
    ```
    
    ```r
    str_split_1("DragonForce", "(?<=[a-z])(?=[A-Z])")
    ```
    
    ```
    [1] "Dragon" "Force" 
    ```
    
    ```r
    str_split_i(c("DragonForce", "HammerFall"), "(?<=[a-z])(?=[A-Z])", 1)
    ```
    
    ```
    [1] "Dragon" "Hammer"
    ```
:   `str_split_1(string, pattern)` は1つの文字列を受け取ってvectorを返す簡略版。
:   `str_split_i(string, pattern, i)` は分割してできたi番目の要素だけをvectorで返す亜種。
:   `str_split_fixed(string, pattern, n)` は
    `n` 必須、 `simplify = TRUE` 固定でmatrixを返すショートカット。
:   data.frame内の文字列を分割したい場合は
    [`tidyr::separate*()`]({{< relref "tidyr.md" >}}) 系の関数を使う。

`str_dup(string, times)`
:   指定した回数だけ文字列を繰り返して結合。
    `base::strrep()` と同等。
    
    ```r
    str_dup("pizza", 10)
    ```
    
    ```
    [1] "pizzapizzapizzapizzapizzapizzapizzapizzapizzapizza"
    ```

### Pattern Matching

`str_count(string, pattern)`
:   マッチする箇所の数を返す。

`str_detect(string, pattern, negate = FALSE)`
:   マッチするかどうか `logical` を返す。
    `negate = TRUE` で結果を反転。
    `base::grepl(pattern, x)` と相同。
:   正規表現を覚えてなくても始まりと終わりだけ手軽にマッチできる
    `str_starts()`, `str_ends()` もある。

`str_extract(string, pattern)`, `str_extract_all(string, pattern)`
:   マッチした部分文字列を取り出す。しなかった要素には `NA`。
:   数値＋単位のような文字列から数値部分だけを抜き出すには
    [readr::parse_number()]({{< relref "readr.md#parse" >}})
    が便利。

`str_subset(string, pattern, negate = FALSE)`
:   `x[str_detect(x, pattern)]` のショートカット。
    マッチする要素だけ元の形で返すので
    `str_extract()` より `base::grep(pattern, x, value = TRUE)` に近い。

`str_which(string, pattern, negate = FALSE)`
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
    `str_remove()` はマッチした部分を消すためのショートカット。

------------------------------------------------------------------------

上記関数の`pattern`引数は普通に文字列を渡すと正規表現として解釈してくれるが、
下記の関数を通して渡すことでその挙動を変更することができる。

`stringr::regex(pattern, ignore_case = FALSE, multiline = FALSE, comments = FALSE, dotall = FALSE, ...)`
:   デフォルトの[ICU正規表現](https://unicode-org.github.io/icu/userguide/strings/regexp.html)。
    複数行ファイルに対するマッチではこの関数を通して挙動をいじることになる。

`stringr::fixed(pattern)`
:   正規表現ではなくそのままの文字としてマッチさせる

`stringr::boundary(type = "character", skip_word_none = NA, ...)`
:   境界に対するマッチ。
    `type`の選択肢は `character`, `line_break`, `sentence`, `word`.

`stringr::coll(pattern, ignore_case = FALSE, locale = NULL, ...)`
:   よくわからないけど非ascii対策？


### Formatting

`str_to_upper()`, `str_to_lower()`, `str_to_title()`, `str_to_sentence()`
:   大文字・小文字の変換

`str_glue(..., .sep = "", .envir = parent.frame())`
:   渡された文字列の中の `{R表現}` を評価して埋め込む。
    `sprintf()` よりも使い方が簡単。ライバルは `paste0()` とか。
    `str_interp()` はこれに取って代わられた。
    
    ``` r
    str_glue("fruit[1] is {fruit[1]}.")
    ```
    
    ```
    fruit[1] is apple.
    ```
:   data.frameの流れるpipe上では `str_glue_data()` が便利:
    
    ``` r
    mtcars |> str_glue_data("mean(disp) is {mean(disp)}.")
    ```
    
    ```
    mean(disp) is 230.721875.
    ```
:   本家の [`library(glue)`](https://glue.tidyverse.org/) にはほかのオプションもある。
    それでもPythonのf-stringのような簡易フォーマッタは無くてちょっと不便。

`str_pad(string, width, side = c("left", "right", "both"), pad = " ")`
:   文字列の幅を `width` に伸ばして `side` 側を `pad` で埋める。
    
    ``` r
    str_pad(c("9", "10"), 3L, pad = "0")
    ```
    
    ```
    [1] "009" "010"
    ```

`str_trim(string, side = "both")`
:   端の空白文字を除去する。
    Python でいうところの `string.strip()`。
:   `str_squish()` は両端trimしたうえに内部の連続する空白文字を1つに縮める亜種。
    
    ``` r
    str_trim("   trim   me   ")
    ```
    
    ```
    [1] "trim   me"
    ```
    
    ``` r
    str_squish("   trim   me   ")
    ```
    
    ```
    [1] "trim me"
    ```

`str_trunc(string, width, side = c("right", "left", "center"), ellipsis = "...")`
:   一定の長さを超えたら捨てて `...` にする。

`str_wrap(string, width = 80, indent = 0, exdent = 0)`
:   指定した幅で折り返す。
    `indent` は先頭行の左余白。
    `exdent` はそれ以外の行の左余白。

文字列と数値の型変換はstringrの管轄外なので、標準の
`as.character()` や `as.double()` などを使うか、
[`readr::parse_*()`系の関数]({{< relref "readr.md#parse" >}})
を使う。


## Rの文字列と正規表現

ダブルクォーテーションで挟んで作る。
文字列の中に `"` を含む場合はシングルクォーテーションで挟む。

```r
s = "This is a string."
s = 'This is a string with "double quotes".'
```

### エスケープシーケンス

バックスラッシュを使って改行 `\n` やタブ `\t` などの制御文字を表現できる。
バックスラッシュ自体を表すためには `\\` のように重ねる必要がある。


``` r
string = "x\ty\n0\t1\n"
print(string)
```

```
[1] "x\ty\n0\t1\n"
```

``` r
cat(string)
```

```
x	y
0	1
```

``` r
readr::read_tsv(I(string))
```

```
  x y
1 0 1
```

See [`?Quotes`](https://stat.ethz.ch/R-manual/R-patched/library/base/html/Quotes.html)



### 正規表現

[ICU正規表現](https://unicode-org.github.io/icu/userguide/strings/regexp.html)からよく使うやつを抜粋。

| メタ文字 | 意味 |
| ---- | ---- |
| `\d` | 数字   |
| `\s` | 空白   |
| `\w` | 英数字 |
| `.`  | 何でも |
| `^`  | 行頭   |
| `$`  | 行末   |

`\D`, `\S`, `\W` のように大文字にすると反転してそれ以外にマッチ。

| 演算子 | 意味 |
| ---- | ---- |
| `?`  | 0回か1回 |
| `*`  | 0回以上繰り返し |
| `+`  | 1回以上繰り返し |
| `{n,m}` | n回以上m回以下 |
| `XXX(?=YYY)`  | YYYに先立つXXX |
| `(?<=YYY)XXX`  | YYYに続くXXX |


### 生文字列

数字にマッチする正規表現を書こうとして `pattern = "\d"` とすると怒られる。
先述のようにバックスラッシュそのものを表すには二重にしておく必要があるため。
```r
"\d"
# Error: '\d' is an unrecognized escape in character string starting (<input>:1:3)

"\\d"
# Good.
```

エスケープシーケンスを無効にした生文字列(raw string)を用いることでバックスラッシュを重ねずに済む。
PythonやC++などでは前からあったけどRでもようやく4.0.0 から使えるようになった。

```r
pattern = "\\d"
pattern = r"(\d)"
pattern = R"(\d)"
pattern = r"---(\d)---"
pattern = r"---[\d]---"
pattern = r"---{\d}---"
stringr::str_count("1q2w3e4r", pattern)
```


## 関連書籍

<a href="https://www.amazon.co.jp/dp/1492097403?&linkCode=li3&tag=heavywatal-22&linkId=163b4c2d2d4f43d197e985a033d397c1&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=1492097403&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=1492097403" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311814X?&linkCode=li3&tag=heavywatal-22&linkId=a289b1f9dbb4f189b4209b374662d6f7&language=ja_JP&ref_=as_li_ss_il" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22&language=ja_JP" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&language=ja_JP&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
