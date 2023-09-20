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


```r
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

```r
start(ir1)
```

```
[1] 1 2 3 4
```

```r
end(ir1)
```

```
[1] 3 3 3 3
```

```r
width(ir1)
```

```
[1] 3 2 1 0
```


## `GenomicRanges`


```r
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
`rtraclayer::import("annotation.gff3")`,
`GenomicFeatures::genes(txdb)`,
のような形で読み込むことが多い。

個々の区間の情報へのアクセスはIRangeと同じ
`start(gr)`, `end(gr)`, `width(gr)`
に加えて:

```r
seqnames(gr)
```

```
factor-Rle of length 3 with 2 runs
  Lengths:    1    2
  Values : chr2 chr1
Levels(2): chr1 chr2
```

```r
ranges(gr)
```

```
IRanges object with 3 ranges and 0 metadata columns:
      start end width
  [1]   101 200   100
  [2]   102 201   100
  [3]   103 202   100
```

```r
strand(gr)
```

```
factor-Rle of length 3 with 3 runs
  Lengths: 1 1 1
  Values : - + *
Levels(3): + - *
```

```r
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

```r
seqinfo(gr)
```

```
Seqinfo object with 2 sequences from an unspecified genome:
  seqnames seqlengths isCircular genome
  chr1      249250621         NA   <NA>
  chr2      243199373         NA   <NA>
```

```r
seqlevels(gr)
```

```
[1] "chr1" "chr2"
```

```r
seqlengths(gr)
```

```
     chr1      chr2 
249250621 243199373 
```

```r
isCircular(gr)
```

```
chr1 chr2 
  NA   NA 
```

```r
genome(gr)
```

```
chr1 chr2 
  NA   NA 
```



## Functions

多くの関数は `IRanges` と `GenomicRanges` で同様に動作する。
より単純な前者で例を示し、後者固有の話は適宜挟む。


```r
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

![plot of chunk iranges-orig](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAZAAAABkCAMAAACfFZZFAAAACXBIWXMAAA9hAAAPYQGoP6dpAAAC5VBMVEUAAAD///8AAP+AgP+AgICqqqpVVf+AgIBAgP+ZmZkzZsyAgICSkpJJbduAgIBAYN85ceNNZuaAgIA7YthJbduPj49AcN88aeGOjo6GhoZDa+SMjIxAZuaGhoZDZN6KiopAat89ZuCJiYmOjo6JiYk+auWIiIhCa96Hh4dAaN+Hh4dEaeGKiopCZuKOjo6Hh4eGhoZBad6GhoZAZt+IiIiFhYWLi4uIiIhBZ96Kioo+aOBCa+CHh4dAZ+KHh4eJiYmGhoZBaOGIiIhAauGKioqIiIhCa+KKioqHh4dAaN+JiYmJiYmJiYlDZ+GIiIiHh4eIiIhBZ+CIiIhAaOGIiIhCauGJiYlAauGHh4dCaeFAauKIiIhBaeCJiYmIiIiJiYlCaOGJiYmHh4dAauJCaeJAauCHh4eHh4eIiIhCaeKIiIiJiYmIiIhCaeGIiIiIiIiHh4eIiIiHh4eHh4eIiIhCaOGIiIhBaeKJiYlCaeCJiYlBaeGJiYmHh4eHh4dBauKIiIhAaeKHh4eIiIhBaeFAauGJiYmIiIhBaOFBaeFAaeKJiYmIiIiIiIhCaOFBauGIiIhBaOCIiIhBaeCHh4lAaeGGiIpBaeGIiIhBaOGIiIhBaeCIiIhCaOGHh4eChZBBauGJiYlBaOBBaeGIiIiIiIh/hZSIiIiIiIiIiIiHh4d9g5ZBaeGIiIiHh4dCaeGJiYmHh4dBaeBBaeGIiIhAauGIiIhBaOGIiIh5gZpAaeGIiIiIiIhBaeGJiYlBaeGJiYlBaeKIiIh2gKBBauFBaeF0gKGHh4dBaeCIiIiIiIhzgKKIiIhBaeBxfqRCaeGIiIiIiIhBaeGIiIhBaeFBaeGHh4dBaeGIiIiIiIhBaeGIiIiIiIiHh4eIiIiIiIhtfKpBaeGIiIhBaeGIiIiIiIhBaeFBaeFse6tBaeGIiIiIiIhse6xBaeFBaeGIiIhBaeFqe6xBaeGIiIhre61BaeEAAAAGfE2qAAAA83RSTlMAAQECAgMDBAQFBQYHBwgICQoMDQ4QEBESExMUFBUXGBgZGhscHR4fICAiIiMjJCQmJygoKywuLy8wMTIzNDU4OTs8PD0+Pj9AQEFDRUVHSElKS0xNTVBUVVVXWFpbXF1dX2BgYWNkZmdoaWprbW9yc3R1d3h4enp7fH1+f4CCgoODhIWGh4iJiYqLjI2PkJGSlZaWl5eYmZqanp6fn6KjpKamqKmrq62usLGxsbKzs7W1t7i5ury8vr6+v8HDxMTGx8jLy8/Q0dHS1NTV2dra293d3t7f4OLk5+fr7Ozv8fHx8vP09fb3+Pj5+vr6/P39/v6UpsZkAAAAAWJLR0T23NtKYQAABAhJREFUeNrt1/l33FMcxvFnYkpTUyqIohVbUxSpLVQl9jVhUFVCEFURrWWotcqU2ovp1Ja2lpbaxVprUGupBrVTW6tRNQj95v7OOTqZcfQkd+SZb2dun9cf8Jz7Oe/zPXMGIiIiIiIiq9AhTbsjNxVf/YSdmcPhkFqvCjmp91vG1vIauCNngxxq7L0Gd+RskOONvRa4I2eD9PvRWLsB7sjZINi3xdj5aVovuCN3gwCbbWulN1ySy0FWSwpCoiCOUhASBXGUgpAoiKMUhERBHKUgJAriKAUhURBHKQiJgjhKQUgUxFEKQqIgjlIQEgVxFCnIukMrbAztg3RbVPhnUAHyAifIRUuNnWUTCpC0ySvGTx/ugnxACTLc2KtB0mPGX5/lxTdCCTLZ2LsLK/Rcany2E/IAJcgFxt4VSPrY+Ov3DZAHKEEGfGFsfT0ISaOMv25GPqAEwcAJT8+28czEHZBy1LTZ/rlvZF78hOh/CIuCOEpBSBTEUQpCoiCOUhASBXGUgpAoiKMUhERBHKUgJAriKAUhURBHKQiJgjhKQUgUxFEKQqIgjvpPkMOvb8yCcXsgpV/DHY02bhtVDF+td8rkxn9M3B9dWn/klMZM3DJiTWQeZJzJjt+GIWnzFmPr7Y3goz4vmZSGLnu8bjJ17xoZByn+w2TJu/+r+Vnw0QiTZlFPdO40k7k9Mw5SYbJlWcf3OtPYmwIfXWLSDUDnrjOZOzHjIL0WmSx5HkkNxl4NfHSQSfNpATo3zGRs+UCbII/H0z36a3tW/HBPPOn2D9ptvX9r3E9vtnf4+YF4V95pz9CfL8b/5UysRP8ZnqwiTViZtUIrPPVtiGZHLxqiudHbKkTTOitEc7pXFeqOQnSqaSFoSr3xoLnW6w+a1vtBc6p3ILpFQfIqSNUJoCmqGwKaverWAU3twaApqyuBiIgARTOWPFkCiqpzaINF0xe/dxhpq/yN1mfLiZc+eETaFt/LdcHKL3ug+wKl8yK0web6HmWtW1O2ihYPDlYmQrRLD/DCqS2+0oUBoLkS3Xe350VYg32XBIBJZ1O29psO4KMy1qWFC+aEU1t8R84CEK0Hw/kR1mB5M4ALo7THlX9XyLp07MmxcGqLb/RVAOqjYIhEmIOFn+9D2tq0zTuX9bAt5wZj4dQW35jLAZwxnhSEOFgyJx5gbQW2/2ZvcMae2xmxcGqL75ipAK6sJQWhDQbP+6qa+bhLL+OMVceAWDi1xbfNJwBe3Y0UhDY49c61WVtjRwOIjuGMXdOWSHhtD3Vs8QXmHhuonh8kBWENbjcvANbWrgs2xuDvS3iX3nR02hbfhg8nHukLUhDW4Ele4m8XU7Zw3PxfXhhCvHRSOLklIiIiIiIiIiLZ8xefdJVMvduB9AAAAABJRU5ErkJggg==)

### Intra-range methods

個々の区間を操作

`shift(x, shift = 0L, use.names = TRUE)`
: ずらす。

`narrow(x, start = NA, end = NA, width = NA, use.names = TRUE)`
: 狭める。startには正、endには負の値を与える。
: start = 1, end = -1 のとき何もしないというのがわかりにくすぎて怖い。
  
  ```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
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


```r
IRanges(101:200) + 100
```

```
IRanges object with 1 range and 0 metadata columns:
      start end width
  [1]     1 300   300
```

掛け算はズームイン・ズームアウト:

```r
IRanges(101:200) * 2
```

```
IRanges object with 1 range and 0 metadata columns:
      start end width
  [1]   126 175    50
```

```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
  gaps(ir)
  ```
  
  ```
  IRanges object with 2 ranges and 0 metadata columns:
        start end width
    [1]    30  33     4
    [2]    36  39     4
  ```
: GRangesに対しては染色体全体からの差分がstrandごとに計算される:
  
  ```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
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
  
  ```r
  countOverlaps(ir, ir)
  ```
  
  ```
  [1] 2 2 3 3 3 1 1
  ```

`overlapsAny(query, subject, maxgap = -1L, minoverlap = 0L, type, ...)`
: ひとつでもsubjectと重なるものがあるか、queryと同じ長さの論理vectorを返す。
  
  ```r
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
  
  ```r
  poverlaps(ir, rev(ir))
  ```
  
  ```
  [1] FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE
  ```

`coverage(x, shift=0L, width=NULL, weight=1L, method=c("auto", "sort", "hash", "naive"))`
: 重なりの数をRle型で。
  
  ```r
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
  
  ```r
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
  
  ```r
  nearest(ir, ir)
  ```
  
  ```
  [1] 1 1 2 3 3 6 7
  ```

`precede(x, subject, select = c("first", "all"))`
: subjectのうちどれの上流にあるか。重なるものは除外。
  
  ```r
  precede(ir, ir)
  ```
  
  ```
  [1]  3  3  6  6  6  7 NA
  ```

`follow(x, subject, select = c("last", "all"))`
: subjectのうちどれの下流にあるか。重なるものは除外。
  
  ```r
  follow(ir, ir)
  ```
  
  ```
  [1] NA NA  2  2  2  4  6
  ```

`distanceToNearest(x, subject, select = c("arbitrary", "all"))`
: 最も近いsubjectへの距離。overlapしている場合はゼロ。
  
  ```r
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
: `distance(x, y)` はelement-wizeに計算。
  
  ```r
  distance(ir, rev(ir))
  ```
  
  ```
  [1] 27 20  0  0  0 20 27
  ```


### `ignore.strand`

GenomicRangesの操作はデフォルトでseqnamesとstrandごとに分けて行われる。
strandを無視するには `ignore.strand = TRUE` を渡す。


```r
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

```r
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
