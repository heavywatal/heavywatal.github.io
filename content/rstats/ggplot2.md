+++
title = 'ggplot2'
subtitle = "きれいなグラフを簡単に合理的に"
tags = ["r", "graph", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -79
+++

<a href="https://ggplot2.tidyverse.org/">
<img src="https://ggplot2.tidyverse.org/logo.png" align="right" width="120" height="139">
</a>

"The **G**rammer of **G**raphics" という体系に基づいて設計されたパッケージ。
単にいろんなグラフを「描ける」だけじゃなく「一貫性のある文法で合理的に描ける」。

Rのグラフ描画システムには`graphics`と`grid`の2つが存在しており、
R標準の`boxplot()`や`hist()`などは前者の上に、
本項で扱う`ggplot2`は後者の上に成り立っている。
使い方が全く異なるので、前者を知らずにいきなりggplot2から始めても大丈夫。

[tidyverse](https://tidyverse.tidyverse.org/) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。

- [初学者向け講義資料2018](/slides/nagoya2018/2-ggplot.html)
- <https://ggplot2.tidyverse.org>
- <http://r4ds.had.co.nz/data-visualisation.html>
- <http://r4ds.had.co.nz/graphics-for-communication.html>
- <http://www.cookbook-r.com/Graphs/>
- [version 2.0での変更点](https://blog.rstudio.com/2015/12/21/ggplot2-2-0-0/)
- [version 3.0での変更点](https://www.tidyverse.org/articles/2018/07/ggplot2-3-0-0/)

## 基本的な使い方: 指示を `+` していく

- `ggplot()` このデータでよろしく
- `geom_*()` 点や線をよろしく
- `theme_*()` 軸とか背景の見た目をよろしく

ggplot2についてくる`mpg`データを例に:
```r
library(tidyverse)
ggplot(data = mpg) +
  geom_point(mapping = aes(x = displ, y = cty)) +
  theme_classic(base_size = 20, base_family = "Helvetica")
```

途中経過をオブジェクトとして取っておける:
```r
p0 = ggplot(mpg, aes(x = displ, y = cty))
p1 = p0 + geom_point()
p2 = p1 + theme_classic(base_size = 20, base_family = "Helvetica")
p3 = p2 + stat_smooth(method = lm, formula = y ~ log(x))
print(p3)
```

画像ファイルに保存するところまできっちり書く:
```r
ggsave("mpg-displ-cty.png", p3, width = 4, height = 4, dpi=300)
```

`ggplot()` に渡すデータは、1行が1観測、1列が1変数という形の
[**整然データ**]({{< relref "programming.md#tidyverse" >}})。<br>
例: `mpg`, `mtcars`, `diamonds`


## [Aesthetic mapping](https://ggplot2.tidyverse.org/reference/aes_group_order.html)

データと見せ方を紐付ける。
`aes(colour = Species)` のように `aes()` 内で列名を指定すると、
その変数に応じて色やサイズなどを変えることができる。
言い換えると、データの値を色やサイズに変換する(スケールする)ことに相当する。
`aes()` の外で指定するとデータによらず全体に反映される:

```r
# データによって点のサイズ・色・形を変える
p0 + geom_point(mapping = aes(x = displ, y = cty, size = cyl,
                              colour = class, shape = drv))

# サイズは常に6、色はオレンジ、不透明度は0.4
p0 + geom_point(mapping = aes(x = displ, y = cty),
                size = 6, colour = "darkorange", alpha = 0.4)
```

- [色・透明度を変える](https://ggplot2.tidyverse.org/reference/aes_colour_fill_alpha.html)
  - `colour`: 点や線の色
  - `fill`: 塗りつぶしの色
  - [`alpha`](https://ggplot2.tidyverse.org/reference/scale_alpha.html): 不透明度 (0が透明、1が不透明)
- [大きさ・形を変える](https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html)
  - [`size`](https://ggplot2.tidyverse.org/reference/scale_size.html): 点や線のサイズ
  - [`shape`](https://ggplot2.tidyverse.org/reference/scale_shape.html): 点の形
  - [`linetype`](https://ggplot2.tidyverse.org/reference/scale_linetype.html): 線
- 単にグループ分けする
  - `group`: 色や形はそのままにグループ分け。反復試行の折れ線グラフなどに。


### `scale_*()` 関数で調整

変数をどうやって色やサイズに反映させるか、
各項目に対応する `scale_*()` 関数で調整する。

```r
ggplot(mpg) +
  geom_point(aes(x = displ, y = cty, colour = class)) +
  scale_colour_brewer(palette = "Spectral")
```

[`scale_*_viridis_c`](https://ggplot2.tidyverse.org/reference/scale_viridis.html)
:   for `color`, `fill`
:   色覚多様性が考慮されたパレットなので、連続値ならまずこれを使う。
    元は[別パッケージ](https://github.com/sjmgarnier/viridis)が必要だったけど
    ggplot2 v3.0 から標準装備になった。
:   `option=`: `viridis`, `magma`, `inferno`, `plasma`
:   連続値には `_c`、離散値には `_d`。

[`scale_*_brewer`](https://ggplot2.tidyverse.org/reference/scale_brewer.html)
:   for `colour`, `fill`
:   いい感じに考えられたパレット [Colorbrewer](http://colorbrewer2.org/)
    から選んで指定するだけなので楽ちん。
    特に離散値でお世話になる。
    利用可能なパレットは `RColorBrewer::display.brewer.all()` でも一覧できる。
:   `scale_*_brewer(..., palette='Blues', direction=1)`: 離散値
:   `scale_*_distiller(..., palette='Blues', direction=-1)`: 連続値

[`scale_*_gradient`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
:   for `colour`, `fill`
:   グラデーションの基準となる色を指定する。
:   `scale_*_gradient(..., low, high, ...)`: 普通の連続値に
:   `scale_*_gradient2(..., low, mid, high, midpoint=0, ...)`: ある中央値を挟んで上下に分けたいとき
:   `scale_*_gradientn(..., colours, values=NULL, ...)`: 多色のヒートマップなどに e.g., `colours=c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')`

[`scale_*_identity`](https://ggplot2.tidyverse.org/reference/scale_identity.html)
:   for `colour`, `fill`, `size`, `shape`, `linetype`, `alpha`
:   色の名前やサイズなどを示す列が予めデータに含まれている場合にそのまま使う。
    [点や線の種類](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html)

[`scale_*_manual`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
:   for `colour`, `fill`, `size`, `shape`, `linetype`, `alpha`
:   対応関係を `values` 引数で直に指定する。

[`scale_size`](https://ggplot2.tidyverse.org/reference/scale_size.html)
:   デフォルトではpointの面積を値にほぼ比例させるが、面積0にはならない。
    値0に面積0を対応させるには `scale_size_area()` を使う。
    半径を比例させるには `scale_radius()` があるけど要注意。

[`scale_alpha(..., range=c(0.1, 1))`](https://ggplot2.tidyverse.org/reference/scale_alpha.html)\
[`scale_linetype`](https://ggplot2.tidyverse.org/reference/scale_linetype.html)\
[`scale_shape`](https://ggplot2.tidyverse.org/reference/scale_shape.html)


### スケール共通オプション

値と見え方の対応関係が凡例(legend/colourbar)として表示される。
`scale_*()` 関数に以下のオプションを指定することで設定を変えられる。
[連続値の場合](https://ggplot2.tidyverse.org/reference/continuous_scale.html) と
[離散値の場合](https://ggplot2.tidyverse.org/reference/discrete_scale.html)
で微妙に意味が変わるけどだいたいこんな感じ:

- `name`: 凡例のタイトル。複数スケールで同じ名前にすると凡例が統合される。
- `breaks`: 目盛りや凡例に登場させる値。
- `labels`: breaksの値に対応する文字列。デフォルトはbreaksそのまま。
- `na.value`: 欠損値のときどうするか
- `limits`: 数値なら最大値と最小値のvector。
  文字列なら表示したいすべての値(順序も反映される)。
- `guide`: 文字列で `"legend"` か `"colourbar"`。
  さらに細かく制御したい場合は
  [`guide_legend()`](https://ggplot2.tidyverse.org/reference/guide_legend.html) や
  [`guide_colourbar()`](https://ggplot2.tidyverse.org/reference/guide_colourbar.html)
  で。
  消したい場合は `"none"` か `FALSE` を渡せる。


### 変数によってパネルを分割する

多変量データを俯瞰するには、データに応じたパネル分割も便利。
色・サイズなどと合わせれば、x軸y軸プラス3次元程度はパッと可視化できることになる。

[`facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
:   1変数で分割して並べる
    ```r
    facet_wrap(facets, nrow=NULL, ncol=NULL, scales='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE,
               dir='h', strip.position='top')

    p1 + facet_wrap(~ class, ncol = 4L)
    ```

[`facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)
:   2変数で分割して縦横に並べる
    ```r
    facet_grid(facets, margins=FALSE, scales='fixed', space='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE)

    p1 + facet_grid(cyl ~ class)
    ```
    3変数以上にしたい場合は `+` で追加できる。
    1変数でいい場合は `. ~ class` のように片方をドットにする。

[ファセットラベルの調整](https://ggplot2.tidyverse.org/reference/labellers.html)
:   デフォルトでは値だけがfacetラベルに表示されるが、
    変数名を同時に表示したり、数式を表示したりもできる。
    見た目の調整はテーマの `strip.*` で。


### 内部変数を使う

https://github.com/hadley/ggplot2-book/blob/master/layers.rmd#generated-variables

ヒストグラムや箱ヒゲなどの表示に必要な計算は
ggplot内部で `stat_*()` を通して行われる。
そうした値 (computed variables) の一部は
`..count..` や `..density..` みたいに
ピリオドで囲まれた特殊な名前で参照することができる。


## 座標軸やタイトルを変更

[軸の区切りを変更したり対数にしたり](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
:   `gp + scale_x_continuous(breaks=seq(10, 100, by=10))`
:   `gp + scale_y_log10("Beer consumption")`
:   `gp + scale_y_reverse()`
:   上記の<a href="#スケール共通オプション">スケール共通オプション</a>に加えて:
    - `expand`: デフォルトでは値域よりも少し余裕を持たせてあるが
      `c(0, 0)` を与えるとピッタリになる。`geom_tile()`を使うときなどに。
    - `trans`: 数値の変換。exp, identity, log, log10, reciprocal, reverse など。
      文字列変数の順序を変えたい場合は `limits` のほうを使う。
    - `position`: top, bottom, left, right
    - `sec.axis`: 第二軸

[描画する範囲を指定](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
:   `gp + ylim(0, 42) + xlim("b", "c", "d")`
:   `gp + coord_cartesian(xlim = NULL, ylim = NULL)`
:   前者はデータそのものを切るが、後者はデータを変えずに描画領域だけ切る

[X軸とY軸の比率を固定](https://ggplot2.tidyverse.org/reference/coord_fixed.html)
:   `gp + coord_fixed(ratio=1)`

[XY軸の反転](https://ggplot2.tidyverse.org/reference/coord_flip.html)
:   `gp + coord_flip()`

[極座標](https://ggplot2.tidyverse.org/reference/coord_polar.html)
:   パイチャートも作れるらしい

[座標変換](https://ggplot2.tidyverse.org/reference/coord_trans.html)
:   `gp + coord_trans(x='log10', y='sqrt')`
:    表示する座標を変換する。
     stat前に適用される `scale_x_*` とは微妙に違う。

[軸ラベルとタイトル](https://ggplot2.tidyverse.org/reference/labs.html)
:   `gp + labs(x="time", y="weight", title="growth", tag="A")`
:   `gp + xlab("time") + ylab("weight") + ggtitle("growth")`


## `theme`: 背景やラベルの調整

<https://ggplot2.tidyverse.org/reference/ggtheme.html>

### 既成テーマ

`theme_grey(base_size=11, base_family='')`, `theme_gray(...)`
:   灰色背景に白い格子。
    ggplotらしいデフォルトだが、論文には使いにくい。

`theme_bw(base_size=11, base_family='')`
:   黒枠白背景にうっすら灰色格子

`theme_linedraw(base_size=11, base_family='')`
:   細いけど濃い色の `panel.grid`

`theme_light(base_size=11, base_family='')`
:   それを薄くした感じ

`theme_minimal(base_size=11, base_family='')`
:   外枠なしの `theme_bw`

`theme_classic(base_size=11, base_family='')`
:   xy軸がL字に描かれているだけで枠もグリッドも無し

`theme_void(base_size=11, base_family='')`
:   完全に枠なし

これらをカッコ無しでコンソールに打ち込むと、
下記の各エレメントの設定方法やデフォルト値を知ることができる。

引数として `base_family` に"Helvetica Neue"などのフォントを指定できる。
Macなら"HiraKakuProN-W3"を指定すれば日本語でも文字化けしなくなるはず。
テーマを構成する `axis.text` などはこの設定を継承するが、
`geom_text()` などプロット内部の要素には引き継がれないことに注意。

ほかにもいろんなテーマが
[ggthemes](https://cran.r-project.org/web/packages/ggthemes/vignettes/ggthemes.html)
というパッケージで提供されている。


### 設定項目

<https://ggplot2.tidyverse.org/reference/theme.html>

`theme()` 関数に項目と値を指定したものを、
ほかのレイヤーと同じようにどんどん足しながら変更していく。
テーマは直線、長方形、文字の3種類のエレメントからなり、
それらの性質を変更する場合は `element_***()` を介して行う

```r
## ベースとなるテーマを先に適用してから
gp = gp + theme_bw(base_family='HiraKakuProN-W3', base_size=14)
gp = gp + theme(legend.position='bottom')
gp = gp + theme(plot.background=element_rect(fill="transparent"))
```

全体
: `line`: (`element_line`)
: `rect`: (`element_rect`)
: `text`: (`element_text`)
: `title`: (`element_text`; inherits from `text`)
: `aspect.ratio`:

軸タイトル、軸ラベル、目盛
: `axis.title`: (`element_text`; inherits from `text`)\
  &emsp;`__.x`, `__.x.top`, `__.y`, `__.y.right`
: `axis.text`: (`element_text`; inherits from `text`)\
  &emsp;`__.x`, `__.x.top`, `__.y`, `__.y.right`
: `axis.ticks`: (`element_line`; inherits from `line`)\
  &emsp;`__.x`, `__.y`
: `axis.ticks.length`: (`unit`)
: `axis.line`: (`element_line`; inherits from `line`)\
  &emsp;`__.x`, `__.y`

凡例
: `legend.background`: (`element_rect`; inherits from `rect`)
: `legend.margin`: (`margin`)
: `legend.spacing`:(`unit`)\
   &emsp;`__.x`, `__.y`
: `legend.key`: (`element_rect`; inherits from `rect`)
: `legend.key.size`: (`unit`)
: `legend.key.height`: (`unit`; inherits from `legend.key.size`)
: `legend.key.width`: (`unit`; inherits from `legend.key.size`)
: `legend.text`: (`element_text`; inherits from `text`)
: `legend.text.align`: (number from `0` (left) to `1` (right))
: `legend.title`: (`element_text`; inherits from `title`)
: `legend.title.align`: (number from `0` (left) to `1` (right))
: `legend.position`: (`"left"`, `"right"`, `"bottom"`, `"top"`, `"none"` `c(0, 1)`)
: `legend.direction`: (`"horizontal"` or `"vertical"`)
: `legend.justification`: (`"center"` or `c(0, 1)` のような数値でアンカー位置を指定)
: `legend.box`: (`"horizontal"` or `"vertical"`)
: `legend.box.just`: (`"top"`, `"bottom"`, `"left"`, or `"right"`)
: `legend.box.margin`: (`margin`)
: `legend.box.background: (`element_rect`; inherits from `rect`)
: `legend.box.spacing`:(`unit`)

プロット領域の背景、余白、格子
: `panel.background`: (`element_rect`; inherits from `rect`)
: `panel.border`: (`element_rect`; inherits from `rect`; should be used with `fill=NA`)
: `panel.spacing`: (`unit`; `facet_*` の間隔)\
  &emsp;`__.x`, `__.y`
: `panel.grid`: (`element_line`; inherits from `line`)
: `panel.grid.major`: (`element_line`; inherits from `panel.grid`)\
  &emsp;`__.x`, `__.y`
: `panel.grid.minor`: (`element_line`; inherits from `panel.grid`)\
  &emsp;`__.x`, `__.y`

全体の背景、タイトル、余白
: `plot.background`: (`element_rect`; inherits from `rect`)
: `plot.title`: (`element_text`; inherits from `title`)
: `plot.subtitle`: (`element_text`; inherits from `title`)
: `plot.caption`: (`element_text`; inherits from `title`)
: `plot.margin`: (`unit` with the sizes of the top, right, bottom, and left margins)

`facet` したときのラベル
: `strip.background`: (`element_rect`; inherits from `rect`)
: `strip.placement`: ('inside', 'outside')
: `strip.text`: (`element_text`; inherits from `text`)\
  &emsp;`__.x`, `__.y`

その他
: `complete`: 部分的な変更か、完全なテーマか (`FALSE`)
: `validate`: 毎回チェックするか (`TRUE`)


### エレメント

`element_rect(fill, colour, size, linetype, inherit.blank)` --- 長方形
:   `fill`: 塗りつぶしの色
:   `colour`: 枠の色

`element_line(colour, size, linetype, lineend, arrow, inherit.blank)` --- 直線

`element_text(family, face, colour, size, hjust, vjust, angle, lineheight, margin)` --- 文字
:   `family`: フォントファミリー。 空なら `theme_bw(base_family=...)` などの指定を継承。
:   `face`: (`'plain'`, `'italic'`, `'bold'`, `'bold.italic'`)
:   `hjust`, `vjust`: 水平位置と垂直位置の寄せ方をそれぞれ `[0, 1]` の実数で。
:   `angle`: 角度 `[0, 360]`
:   `margin`: スペース調整を関数 `margin(top, right, bottom, left)` 越しに。

`element_blank()` --- 空
:   消したい要素にはこれを指定する

    ```r
    gp = gp + theme(axis.ticks=element_blank())
    gp = gp + theme(panel.grid=element_blank())
    ```

`rel(x)`
:   デフォルトからの相対値で `size` 引数を指定したいときに。

`grid::unit(x, units, data=NULL)`
:   こちらは絶対指定。`grid` パッケージに入ってる。
    units で使いそうなのは
    `cm`, `mm`, `inches`, `points`, `lines`, `native`

`grid::arrow(angle, length, ends, type)`
:   `axis.line` の `element_line()` にこれを与えて軸を矢印にするとか。

`margin(t=0, r=0, b=0, l=0, unit='pt')`
:   marginクラス

`calc_element(element, theme, verbose=FALSE)`
:   継承などを考慮した上でelementがどんな値にセットされるか確かめる。


## ファイルに書き出す

RStudioやQuartzの保存ダイアログを利用して書き出すこともできるけど、
それだとサイズ調整やファイル名入力を手作業でやることになってしまう。
`ggsave()` をスクリプトに書いておけば何回でも同じ設定で出力できる。

```r
ggsave(filename, plot = last_plot(), device = NULL, path = NULL,
       scale = 1,  width = NA, height = NA,
       units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, ...)
```

- 画像形式はファイル名の拡張子から自動的に判別される
  (e.g., `iris.png`, `iris.pdf`)。
- `width`や`height`を小さくするほど、文字・点・線などの要素が相対的に大きくなる。
- `dpi`を変えることで、見た目のバランスを保ったまま解像度を変えられる。
  (これはPNGなどラスタ形式だけの話。PDFなどのベクタ形式なら気にしなくていい)
- タイトルや軸ラベルの文字サイズを変えたいときはテーマの
  `theme_bw(base_size = 42)` や各要素の `element_text(size = 42)` を使う。
- scaleやunitsを使うのは慣れてからで十分。

```r
# 7inch x 300dpi = 2100px四方 (デフォルト)
ggsave("mpg1.png", p1) # width = 7, height = 7, dpi = 300

# 4     x 300    = 1200  全体7/4倍ズーム
ggsave("mpg2.png", p1, width = 4, height = 4) # dpi = 300

# 2     x 600    = 1200  全体をさらに2倍ズーム
ggsave("mpg3.png", p1, width = 2, height = 2, dpi = 600)

# 4     x 300    = 1200  テーマを使って文字だけ拡大
ggsave("mpg4.png", p1 + theme_bw(base_size = 22), width = 4, height = 4)
# 7inch x 300dpi = 2100px四方 (デフォルト)
ggsave("mpg1.png", p1) # width = 7, height = 7, dpi = 300
# 4     x 300    = 1200  全体7/4倍ズーム
ggsave("mpg2.png", p1, width = 4, height = 4) # dpi = 300
# 2     x 600    = 1200  全体をさらに2倍ズーム
ggsave("mpg3.png", p1, width = 2, height = 2, dpi = 600)
# 4     x 300    = 1200  テーマを使って文字だけ拡大
ggsave("mpg4.png", p1 + theme_bw(base_size = 22), width = 4, height = 4)
```


## プロットの種類

[散布図](https://ggplot2.tidyverse.org/reference/geom_point.html)
:   `gp + geom_point(size=2, alpha=0.3)`
:   重なった点をランダムにばらかしたいときは
    [`geom_jitter()`](https://ggplot2.tidyverse.org/reference/geom_jitter.html)
:   [点の形(shape)一覧](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#sec:shape-spec)

[折れ線グラフ](https://ggplot2.tidyverse.org/reference/geom_path.html)
:   `gp + geom_path(size=2, linetype="dashed")` データ順に結ぶ
:   `gp + geom_line()` x軸上の順で結ぶ
:   `gp + geom_step()` 階段状に結ぶ
:   [線の種類(linetype)一覧](https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#sec:line-type-spec)

[面グラフ](https://ggplot2.tidyverse.org/reference/geom_ribbon.html)
:   `gp + geom_ribbon()` --- yminからymaxの面
:   `gp + geom_area()` --- 0からyの面

[ヒストグラム、密度曲線](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
:   `gp + geom_histogram()` --- 棒グラフ(連続値を`stat_bin()` で区切って)
:   `gp + geom_bar()` --- 棒グラフ(離散値を`stat_count()`で数えて)
:   `gp + geom_freqpoly()` --- 折れ線
:   `gp + geom_density()` --- 密度推定されたスムーズな線
:   `gp + geom_bin2d()` --- 二次元ヒストグラム
:   `gp + geom_hex()` --- 六角形版二次元ヒストグラム

[棒グラフ](https://ggplot2.tidyverse.org/reference/geom_bar.html)
:   `gp + geom_col()`
:   グループ分けする場合のオプション:
    - [`position='stack'`](https://ggplot2.tidyverse.org/reference/position_stack.html):
      縦に積み重ねる (デフォルト)
    - [`position='dodge'`](https://ggplot2.tidyverse.org/reference/position_dodge.html):
      横に並べる
    - `position='fill'`: 縦に積み重ね、高さを1に揃えて割合を示す

[箱ひげ図](https://ggplot2.tidyverse.org/reference/geom_boxplot.html)
:   `gp + geom_boxplot()`
:   `gp + geom_violin()`

[ヒートマップ](https://ggplot2.tidyverse.org/reference/geom_tile.html)
:   `gp + geom_tile(aes(fill=z))`
:   `gp + geom_raster(aes(fill=z))` --- 各タイルの大きさを揃える制約のため高速

[エラーバー](https://ggplot2.tidyverse.org/reference/geom_linerange.html)
:   `gp + geom_errorbar(aes(ymax = y + se, ymin = y - se), width = 0.1)`
:   `gp + geom_linerange(...)`
:   `gp + geom_pointrange(...)`

[関数](https://ggplot2.tidyverse.org/reference/stat_function.html)
:   `ggplot(data.frame(x=c(-4, 4)), aes(x)) + stat_function(fun=dnorm, args=c(0, 1), n=200)`

------------------------------------------------------------------------

[回帰曲線](https://ggplot2.tidyverse.org/reference/geom_smooth.html)
:   `gp + geom_smooth(method=glm, method.args=list(family=poisson), se=FALSE)`

[切片と傾きで直線を描く](https://ggplot2.tidyverse.org/reference/geom_abline.html)
:   `gp + geom_abline(intercept=3, slope=5)`
:   `gp + geom_hline(yintercept=7) + geom_vline(xintercept=11)`

[始点と終点で曲線や矢印を描く](https://ggplot2.tidyverse.org/reference/geom_segment.html)
:   `gp + geom_curve(aes(x, y, xend, yend), curvature = -0.2)`
:   `gp + geom_segment(aes(x, y, xend, yend), arrow=arrow())`
:   矢印の調整は [`grid::arrow()`](https://www.rdocumentation.org/packages/grid/topics/arrow)

[文字列や図形を書き加える](https://ggplot2.tidyverse.org/reference/annotate.html)
:   `gp + annotate("text", x=1:4, y=4:1, label=sprintf("x = %d", 1:4))`
:   テーマの `base_family` は引き継がれないので `family=` で指定すべし。
:   数式を表示するには `label="italic(N[t])"` のような文字列で渡して `parse=TRUE`。
:   データ点に対応する文字列を添えるには
    `gp + geom_text(aes(label=foo))` のほうが適している。
    オプションで `nudge_x=2, nudge_y=2` などとすれば点と重ならないようにずらせる。
    [`position_nudge()`](https://ggplot2.tidyverse.org/reference/position_nudge.html)


## Extensions

- <https://ggplot2.tidyverse.org/articles/extending-ggplot2.html>
- <https://www.ggplot2-exts.org/>

ggplotを拡張するための仕組みがversion 2.0から正式に導入され、
ユーザーが独自の stats や geom を作って登録することが容易になった。

### `gridExtra`

<https://github.com/baptiste/gridextra/wiki>

データによって自動的にパネルを分割するには
`facet_grid()` や `facet_wrap()` を使えばよいが、
関係ない複数の図を1枚に描きたい場合は `grid` や `gtable` の機能を使う必要がある。
`gridExtra` はそのへんの操作を手軽にできるようにしてくれるパッケージ。

"grob" は "grid graphical object" の略。
ggplotオブジェクトと同じように `ggsave()` に渡して保存可能。

```r
grob = gridExtra::arrangeGrob(p1, p2, nrow=2, ncol=1, bottom='Time')
grid.newpage()
grid.draw(grob)
```

複数ページのPDFに書き出したい場合は
[purrr]({{< relref "purrr.md" >}}) などを使って list of ggplots を作っておき、
`marrangeGrob()` に渡す。
```r
.grobs = purrr::map(.dataframes, my_ggplot_func)
.gtable = gridExtra::marrangeGrob(.grobs, nrow=4, ncol=3)
ggsave('multi_page.pdf', .gtable, width=7, height=9.9)
```

### `cowplot`

- https://github.com/wilkelab/cowplot
- https://cran.r-project.org/web/packages/cowplot

ggplotを学術論文向けにカスタマイズしやすくする。
主な利用目的はgridExtraと同じでggplotを並べる機能。
論文figureのようなA, B, Cラベルをオプションで簡単に付けられるのが良い。

`cowplot::plot_grid()`
:   `facet_wrap()`のように、ざっと並べるのに便利。
    もちろん入れ子も可能。

    ```r
    cowplot::plot_grid(..., plotlist=NULL,
        align=c('none', 'h', 'v', 'hv'),
        nrow=NULL, ncol=NULL,
        scale=1, rel_widths=1, rel_heights=1,
        labels=NULL, label_size=14,
        hjust=-0.5, vjust=1.5)
    ```
    さらに細やかな制御をしたいときは以下の関数を個別に重ねていく。

`cowplot::ggdraw(plot=NULL, xlim=c(0, 1), ylim=c(0, 1))`
:  これの後ろに `+` 演算子で `draw_***()` を足していく。

`cowplot::get_legend(plot)`
:  凡例が共通する図を並べるとき、代表のやつをこれで取っておいてあとで並べる。

`draw_figure_label()`
`draw_grob()`
`draw_label()`
`draw_line()`
`draw_plot()`
`draw_plot_label()`
`draw_text()`

`theme_cowplot()`
`ggsave2()`

パッケージ読み込みと同時に勝手にテーマを変更したり、
`ggsave()` 関数を上書きしたりという問題が過去にはあったが、今は大丈夫。

できあがった図を並べるための新しいパッケージとして
[patchwork](https://github.com/thomasp85/patchwork)
がすごくエレガントで期待大。
ただし、演算子を多用するスタイルは忘れやすく検索しにくい諸刃の剣。


## 関連書籍

<a href="https://www.amazon.co.jp/dp/331924275X/ref=as_li_ss_il?ie=UTF8&qid=1508463669&sr=8-2&keywords=ggplot2&linkCode=li3&tag=heavywatal-22&linkId=4880b062b674960adf66dab4707eed01" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=331924275X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=331924275X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/4873116538/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=27fae9b2fd078aeba62baef727299bdd" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116538&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116538" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/dp/487311814X/ref=as_li_ss_il?ie=UTF8&qid=1508340144&sr=8-1&keywords=r&linkCode=li3&tag=heavywatal-22&linkId=4137d3d3351f8ccab5a93cefdc28fdec" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=487311814X&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=487311814X" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
