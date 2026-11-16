+++
title = 'GenomicRanges'
subtitle = "ゲノム上の範囲やアノテーション"
tags = ["r", "bioconductor"]
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
