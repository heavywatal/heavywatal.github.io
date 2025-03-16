+++
title = 'GenomicRanges'
subtitle = "ゲノム上の範囲やアノテーション"
tags = ["r", "bioconductor"]
[menu.main]
  parent = "rstats"
  weight = -40
+++





<https://bioconductor.org/packages/GenomicRanges>

ゲノム上の座標・範囲を表現するための型として S4 class `GenomicRanges` を提供する。
これは整数の範囲を扱う
[`IRanges`](https://bioconductor.org/packages/IRanges)
を拡張したもの。


## Installation

```r
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(GenomicRanges)
```


## `IRanges`

両端の値を含む閉区間 `[start, end]` なのが厄介。
end = start + width - 1 の関係。
width = 0 のときだけ例外的に end < start になり、両端を含まない。


``` r
ir1 = IRanges(start = 1:4, width = 3:0)
ir2 = IRanges(start = 1:4, end = 3)
ir3 = IRanges(end = 3, width = 3:0)
stopifnot(identical(ir1, ir2), identical(ir1, ir3))
ir1
```

```
IRanges object with 4 ranges and 0 metadata columns:
      start end width
  [1]     1   3     3
  [2]     2   3     2
  [3]     3   3     1
  [4]     4   3     0
```

``` r
start(ir1)
```

```
[1] 1 2 3 4
```

``` r
end(ir1)
```

```
[1] 3 3 3 3
```

``` r
width(ir1)
```

```
[1] 3 2 1 0
```


## `GenomicRanges`


``` r
gr = GRanges(
  seqnames = c("chr2", "chr1", "chr1"),    # 染色体の名前
  ranges = IRanges(101:103, width = 100),  # 座標
  strand = c("-", "+", "*"),
  score = 51:53,                           # 任意のelementMetadata列
  GC = 0.1 * 5:7,                          # 任意のelementMetadata列
  seqinfo = NULL,
  seqlengths = c(chr1 = 249250621, chr2 = 243199373)
)
gr
```

```
GRanges object with 3 ranges and 2 metadata columns:
      seqnames    ranges strand |     score        GC
         <Rle> <IRanges>  <Rle> | <integer> <numeric>
  [1]     chr2   101-200      - |        51       0.5
  [2]     chr1   102-201      + |        52       0.6
  [3]     chr1   103-202      * |        53       0.7
  -------
  seqinfo: 2 sequences from an unspecified genome
```

上記のように自分で作る機会はほぼ無くて、
`GenomicRanges::makeGRangesFromDataFrame(df)`,
`rtracklayer::import("annotation.gff3")`,
`GenomicFeatures::genes(txdb)`,
のような形で読み込むことが多い。

個々の区間の情報へのアクセスはIRangeと同じ
`start(gr)`, `end(gr)`, `width(gr)`
に加えて:

``` r
seqnames(gr)
```

```
factor-Rle of length 3 with 2 runs
  Lengths:    1    2
  Values : chr2 chr1
Levels(2): chr1 chr2
```

``` r
ranges(gr)
```

```
IRanges object with 3 ranges and 0 metadata columns:
      start end width
  [1]   101 200   100
  [2]   102 201   100
  [3]   103 202   100
```

``` r
strand(gr)
```

```
factor-Rle of length 3 with 3 runs
  Lengths: 1 1 1
  Values : - + *
Levels(3): + - *
```

``` r
mcols(gr)
```

```
DataFrame with 3 rows and 2 columns
      score        GC
  <integer> <numeric>
1        51       0.5
2        52       0.6
3        53       0.7
```

`S4Vectors::mcols()` で参照・代入するのは区間ごとのメタデータ。

データセット全体のメタデータとして染色体の長さなども扱う:

``` r
seqinfo(gr)
```

```
Seqinfo object with 2 sequences from an unspecified genome:
  seqnames seqlengths isCircular genome
  chr1      249250621         NA   <NA>
  chr2      243199373         NA   <NA>
```

``` r
seqlevels(gr)
```

```
[1] "chr1" "chr2"
```

``` r
seqlengths(gr)
```

```
     chr1      chr2 
249250621 243199373 
```

``` r
isCircular(gr)
```

```
chr1 chr2 
  NA   NA 
```

``` r
genome(gr)
```

```
chr1 chr2 
  NA   NA 
```



## Functions

多くの関数は `IRanges` と `GenomicRanges` で同様に動作する。
より単純な前者で例を示し、後者固有の話は適宜挟む。


``` r
ir = IRanges(c(1, 8, 14, 15, 19, 34, 40), width = c(12, 6, 6, 15, 6, 2, 7))
ir
```

```
IRanges object with 7 ranges and 0 metadata columns:
      start end width
  [1]     1  12    12
  [2]     8  13     6
  [3]    14  19     6
  [4]    15  29    15
  [5]    19  24     6
  [6]    34  35     2
  [7]    40  46     7
```

![plot of chunk iranges-orig](data:image/png;base64,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)

### Intra-range methods

個々の区間を操作

`shift(x, shift = 0L, use.names = TRUE)`
: ずらす。

`narrow(x, start = NA, end = NA, width = NA, use.names = TRUE)`
: 狭める。startには正、endには負の値を与える。
: start = 1, end = -1 のとき何もしないというのがわかりにくすぎて怖い。
  
  ``` r
  identical(ir, narrow(ir, 1, -1))
  ```
  
  ```
  [1] TRUE
  ```
: `threebands(x, start = NA, end = NA, width = NA)`
  は削られる両端部分も含めて
  `$left`, `$middle`, `$right` のリストで返してくれる亜種。
: 逆に広げるには `flank()` を `punion()` するか、
  `+` 演算子で両側に広げてから片側を `narrow()` するか？

`resize(x, width, fix = "start", use.names = TRUE, ...)`
: 幅を変える。
: `fix`: start, end, center

`flank(x, width, start = TRUE, both = FALSE, use.names = TRUE, ...)`
: start上流もしくはend下流の領域。両方いっぺんには取れない。
: `both = TRUE` は start (or end) を起点に両側という意味であり、範囲の両側ではない。
  
  ``` r
  flank(ir, 1, both = TRUE)
  ```
  
  ```
  IRanges object with 7 ranges and 0 metadata columns:
        start end width
    [1]     0   1     2
    [2]     7   8     2
    [3]    13  14     2
    [4]    14  15     2
    [5]    18  19     2
    [6]    33  34     2
    [7]    39  40     2
  ```
: `promoters(x, upstream=2000, downstream=200, use.names=TRUE, ...)`
  はstart起点に上下異なる幅で取ってこられる亜種。

`reflect(x, bounds, use.names = TRUE)`
: `bounds` の裏から見た相対位置にする。
  
  ``` r
  reflect(ir, IRanges(1, 1000))
  ```
  
  ```
  IRanges object with 7 ranges and 0 metadata columns:
        start  end width
    [1]   989 1000    12
    [2]   988  993     6
    [3]   982  987     6
    [4]   972  986    15
    [5]   977  982     6
    [6]   966  967     2
    [7]   955  961     7
  ```
: いかにも負のstrandの座標処理に使えそうだがなぜかGenomicRangesには未対応。
  自分で書くならこんな感じか:
  
  ``` r
  y = gr
  bounds = IRanges(start = 1, width = seqlengths(gr)[as.vector(seqnames(gr))])
  ranges(y)[strand(gr) == "-"] = reflect(ranges(gr), bounds)[strand(gr) == "-"]
  y
  ```
  
  ```
  GRanges object with 3 ranges and 2 metadata columns:
        seqnames              ranges strand |     score        GC
           <Rle>           <IRanges>  <Rle> | <integer> <numeric>
    [1]     chr2 243199174-243199273      - |        51       0.5
    [2]     chr1             102-201      + |        52       0.6
    [3]     chr1             103-202      * |        53       0.7
    -------
    seqinfo: 2 sequences from an unspecified genome
  ```
: `reverse(x)` は `reflect(x, range(x))` のショートカット。
  区間を逆順に並べる `rev()` とは違う。
  

`restrict(x, start = NA, end = NA, keep.all.ranges = FALSE, use.names = TRUE)`
: `start` から `end` までの範囲のみ残して外を捨てる。境界含む。
: `end = 14` で15から始まる区間が取れてきちゃうのはバグじゃない？
  
  ``` r
  restrict(ir, 10, 14)
  ```
  
  ```
  IRanges object with 4 ranges and 0 metadata columns:
        start end width
    [1]    10  12     3
    [2]    10  13     4
    [3]    14  14     1
    [4]    15  14     0
  ```

足し算・引き算は両側に伸縮:


``` r
IRanges(101:200) + 100
```

```
IRanges object with 1 range and 0 metadata columns:
      start end width
  [1]     1 300   300
```

掛け算はズームイン・ズームアウト:

``` r
IRanges(101:200) * 2
```

```
IRanges object with 1 range and 0 metadata columns:
      start end width
  [1]   126 175    50
```

``` r
IRanges(101:200) * 0.5
```

```
IRanges object with 1 range and 0 metadata columns:
      start end width
  [1]    51 250   200
```

### Inter-range methods

区間の集合を操作

`range(x, ..., with.revmap = FALSE, na.rm = FALSE)`
: 端から端まで1つの区間として返す。
  
  ``` r
  range(ir)
  ```
  
  ```
  IRanges object with 1 range and 0 metadata columns:
        start end width
    [1]     1  46    46
  ```

`reduce(x, drop.empty.ranges = FALSE, min.gapwidth = 1L, with.revmap = FALSE)`
: 重なっている区間をつなげて平らにする。
: `with.revmap = TRUE` とすると入力した区間がどこに含まれるかをmcolsに保持する。
  
  ``` r
  reduce(ir, with.revmap = TRUE)
  ```
  
  ```
  IRanges object with 3 ranges and 1 metadata column:
        start end width |    revmap
    [1]     1  29    29 | 1,2,3,...
    [2]    34  35     2 |         6
    [3]    40  46     7 |         7
  ```
  
  ![plot of chunk iranges-reduce](data:image/png;base64,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)
: `min.gapwidth` を大きくすると離れた区間もつなげられる。
  例えば `reduce(ir, min.gapwidth = 10L)` で全部つながって `range(ir)` と同じになる。

`gaps(x, start=NA, end=NA, ...)`
: `IRanges::setdiff(IRanges(start, end), x)` と同等。
  start/end 省略時は `range(x)` からの差分。
  
  ``` r
  gaps(ir)
  ```
  
  ```
  IRanges object with 2 ranges and 0 metadata columns:
        start end width
    [1]    30  33     4
    [2]    36  39     4
  ```
: GRangesに対しては染色体全体からの差分がstrandごとに計算される:
  
  ``` r
  gaps(gr)
  ```
  
  ```
  GRanges object with 9 ranges and 0 metadata columns:
        seqnames        ranges strand
           <Rle>     <IRanges>  <Rle>
    [1]     chr1         1-101      +
    [2]     chr1 202-249250621      +
    [3]     chr1   1-249250621      -
    [4]     chr1         1-102      *
    [5]     chr1 203-249250621      *
    [6]     chr2   1-243199373      +
    [7]     chr2         1-100      -
    [8]     chr2 201-243199373      -
    [9]     chr2   1-243199373      *
    -------
    seqinfo: 2 sequences from an unspecified genome
  ```


`disjoin(x, ...)`
: 全ての始点・終点を使い、重なり無しで最多の小区間にして返す。
  
  ``` r
  disjoin(ir)
  ```
  
  ```
  IRanges object with 10 ranges and 0 metadata columns:
         start end width
     [1]     1   7     7
     [2]     8  12     5
     [3]    13  13     1
     [4]    14  14     1
     [5]    15  18     4
     [6]    19  19     1
     [7]    20  24     5
     [8]    25  29     5
     [9]    34  35     2
    [10]    40  46     7
  ```
  
  ![plot of chunk iranges-disjoin](data:image/png;base64,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)

`disjointBins(x, ...)`
: 重ならないように描画するときのy座標に使える。
  
  ``` r
  disjointBins(ir)
  ```
  
  ```
  [1] 1 2 1 2 3 1 1
  ```

`union(x, y)`, `intersect(x, y)`, `setdiff(x, y)`
: 結果は `reduce()` 済みの区間。
: element-wise の `punion()`, `pintersect()`, `psetdiff()`, `pgap()` もある。

`GenomicRanges::subtract(x, y, minoverlap = 1L, ...)`
: `setdiff()` と似ているが元の区間ごとの `GRangesList` にして返す。


### Overlaps

`type`
: any: ちょっとでも重なっているか距離が `maxgap` 以下ならOK。
: start/end/equal: 開始/終了/両方の距離が `maxgap` 以下ならOK。
: within: queryがすっぽり含まれていればOK。 `maxgap` の意図・挙動は不明。

`findOverlaps(query, subject, maxgap = -1L, minoverlap = 0L, type, select, ...)`
: `select`: subject側に複数マッチした場合に何を返すか。
  デフォルトの "all" ならHits型オブジェクト、
  "first" や "last" なら整数vectorでインデックスを返す。
  "arbitrary" の挙動は謎だが乱数を振ったりするわけではなさそう。
  
  ``` r
  findOverlaps(ir, ir)
  ```
  
  ```
  Hits object with 15 hits and 0 metadata columns:
         queryHits subjectHits
         <integer>   <integer>
     [1]         1           1
     [2]         1           2
     [3]         2           1
     [4]         2           2
     [5]         3           3
     ...       ...         ...
    [11]         5           3
    [12]         5           4
    [13]         5           5
    [14]         6           6
    [15]         7           7
    -------
    queryLength: 7 / subjectLength: 7
  ```
  
  ``` r
  findOverlaps(ir, ir, select = "first")
  ```
  
  ```
  [1] 1 1 3 3 3 6 7
  ```
: `mergeByOverlaps(query, subject, ...)` と
  `findOverlapPairs(query, subject, ...)`
  はHits型とは違う謎の形式で返す亜種。

`countOverlaps(query, subject, maxgap = -1L, minoverlap = 0L, type, ...)`
: いくつsubjectと重なるか、queryと同じ長さの整数vectorを返す。
  
  ``` r
  countOverlaps(ir, ir)
  ```
  
  ```
  [1] 2 2 3 3 3 1 1
  ```

`overlapsAny(query, subject, maxgap = -1L, minoverlap = 0L, type, ...)`
: ひとつでもsubjectと重なるものがあるか、queryと同じ長さの論理vectorを返す。
  
  ``` r
  overlapsAny(ir, ir[4])
  ```
  
  ```
  [1] FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE
  ```
: `%over%`, `%within%` `%outside%`
  という演算子で短く書けるけど可読性低下とconflictが怖いので使わない。

`subsetByOverlaps(x, ranges, maxgap = -1L, minoverlap = 0L, type, invert = FALSE, ...)`
: `x[overlapsAny(x, ranges)]` のショートカット。
  

`overlapsRanges(query, subject, hits = NULL, ...)`
: subsetByOverlapsした上、重なっている部分のみ切り詰めて返す。

`poverlaps(query, subject, maxgap = 0L, minoverlap = 1L, type, ...)`
: `min()` に対する `pmin()` のように、element-wiseに比較する。
  
  ``` r
  poverlaps(ir, rev(ir))
  ```
  
  ```
  [1] FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE
  ```

`coverage(x, shift=0L, width=NULL, weight=1L, method=c("auto", "sort", "hash", "naive"))`
: 重なりの数をRle型で。
  
  ``` r
  coverage(ir, shift = 200, width = 800)
  ```
  
  ```
  integer-Rle of length 800 with 13 runs
    Lengths: 200   7   5   2   4   1   5   5   4   2   4   7 554
    Values :   0   1   2   1   2   3   2   1   0   1   0   1   0
  ```
: `cvg(x, from = NA, to = NA, weight = 1L, varname = "cvg", collapse = FALSE, ...)`
  も似たような関数だがドキュメント無し。
: `slice(x, ...)` を使うと「カバレッジがこれ以上の区間」みたいなのを簡単に取れる:
  
  ``` r
  coverage(ir) |> IRanges::slice(lower = 2, rangesOnly = TRUE)
  ```
  
  ```
  IRanges object with 2 ranges and 0 metadata columns:
        start end width
    [1]     8  12     5
    [2]    15  24    10
  ```


### Neighboring

`nearest(x, subject, select = c("arbitrary", "all"))`
: 最も近いもの、あるいはoverlapしているもの。
  
  ``` r
  nearest(ir, ir)
  ```
  
  ```
  [1] 1 1 2 3 3 6 7
  ```

`precede(x, subject, select = c("first", "all"))`
: subjectのうちどれの上流にあるか。重なるものは除外。
  
  ``` r
  precede(ir, ir)
  ```
  
  ```
  [1]  3  3  6  6  6  7 NA
  ```

`follow(x, subject, select = c("last", "all"))`
: subjectのうちどれの下流にあるか。重なるものは除外。
  
  ``` r
  follow(ir, ir)
  ```
  
  ```
  [1] NA NA  2  2  2  4  6
  ```

`distanceToNearest(x, subject, select = c("arbitrary", "all"))`
: 最も近いsubjectへの距離。overlapしている場合はゼロ。
  
  ``` r
  distanceToNearest(ir, ir[1:3])
  ```
  
  ```
  Hits object with 7 hits and 1 metadata column:
        queryHits subjectHits |  distance
        <integer>   <integer> | <integer>
    [1]         1           1 |         0
    [2]         2           1 |         0
    [3]         3           2 |         0
    [4]         4           3 |         0
    [5]         5           3 |         0
    [6]         6           3 |        14
    [7]         7           3 |        20
    -------
    queryLength: 7 / subjectLength: 3
  ```
: `distance(x, y)` はelement-wiseに計算。
  
  ``` r
  distance(ir, rev(ir))
  ```
  
  ```
  [1] 27 20  0  0  0 20 27
  ```


### `ignore.strand`

GenomicRangesの操作はデフォルトでseqnamesとstrandごとに分けて行われる。
strandを無視するには `ignore.strand = TRUE` を渡す。


``` r
reduce(gr)
```

```
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   102-201      +
  [2]     chr1   103-202      *
  [3]     chr2   101-200      -
  -------
  seqinfo: 2 sequences from an unspecified genome
```

``` r
reduce(gr, ignore.strand = TRUE)
```

```
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   102-202      *
  [2]     chr2   101-200      *
  -------
  seqinfo: 2 sequences from an unspecified genome
```
