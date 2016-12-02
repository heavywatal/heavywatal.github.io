+++
title = 'ggplot2'
subtitle = "きれいなグラフを簡単に合理的に"
tags = ["r", "graph", "tidyverse"]
[menu.main]
  parent = "rstats"
  weight = -79
+++

Rのグラフ描画システムには`graphics`と`grid`の2つが存在しており、
R標準の`plot()`などは前者の上に、
本項で扱う`ggplot2`は後者の上に成り立っている。
使い方が全く異なるので、前者のことを知る必要はない。

[tidyverse](https://github.com/tidyverse/tidyverse) に含まれているので、
`install.packages('tidyverse')` で一括インストール、
`library(tidyverse)` で一括ロード。

- <http://r4ds.had.co.nz/data-visualisation.html>
- <http://docs.ggplot2.org/>
- <http://www.cookbook-r.com/Graphs/>
- <http://www.rdocumentation.org/packages/ggplot2>
- <https://github.com/hadley/ggplot2-book>

## 基本的な使い方

R に入ってるお馴染みサンプルデータ `iris` を使って

```r
## ライブラリの読み込み
library(tidyverse)

## データと全体設定を持ったggplotオブジェクトを作る
gp = ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length, colour=Species))

## グラフのレイヤーを重ねる
gp = gp + geom_point(size=3, alpha=0.7)

## 描画してみる
print(gp)

## さらにグラフを重ねたり、タイトルやテーマの設定をしたり
gp = gp + geom_smooth(method=glm, family=gaussian)
gp = gp + labs(title="Iris Sepal")
gp = gp + theme_bw()
gp = gp + theme(panel.grid.minor=element_blank())
print(gp)

## ファイルに保存
ggsave("iris_sepal.png", gp)
```

{{%div class="note"%}}
`ggplot()` に渡すデータを整形するには
[dplyr]({{< relref "dplyr.md" >}}) と [tidyr]({{< relref "tidyr.md" >}}) を使う。
{{%/div%}}

## プロットの種類

[version 2.0での変更点](https://blog.rstudio.org/2015/12/21/ggplot2-2-0-0/)

[散布図](http://docs.ggplot2.org/current/geom_point.html)
:   `gp + geom_point(size=2, alpha=0.3)`

[折れ線グラフ](http://docs.ggplot2.org/current/geom_path.html)
:   `gp + geom_path(size=2, linetype="dashed")` データ順に結ぶ\
    `gp + geom_line()` x軸上の順で結ぶ\
    `gp + geom_step()` 階段状に結ぶ

[ヒストグラム、密度曲線](http://docs.ggplot2.org/current/geom_histogram.html)
:   `gp + geom_histogram()` --- 棒グラフ(連続値を`stat_bin()` で切って)
:   `gp + geom_bar()` --- 棒グラフ(離散値を`stat_count()`で数えて)
:   `gp + geom_freqpoly()` --- 折れ線
:   `gp + geom_density()` --- 密度推定されたスムーズな線

[棒グラフ](http://docs.ggplot2.org/current/geom_bar.html)
:   `gp + geom_col()`\
    同じxに対して複数グループのyが存在するとき、
    `position='dodge'` にすると横並び
    (デフォルトは縦積みの`'stack'`)。

[箱ひげ図](http://docs.ggplot2.org/current/geom_boxplot.html)
:   `gp + geom_boxplot()`\
    `gp + geom_violin()`

[ヒートマップ](http://docs.ggplot2.org/current/geom_tile.html)
:   `gp + geom_tile(aes(fill=z))`\
    `gp + geom_raster(aes(fill=z))`\
    後者は各タイルの大きさがすべて同じ場合の特殊ケースで、高速。

[エラーバー](http://docs.ggplot2.org/current/geom_linerange.html)
:   `limits = aes(ymax=height+se, ymin=height-se)`\
    `gp + geom_errorbar(limits, width=0.1)`\
    `gp + geom_pointrange(limits)`\
    `gp + geom_crossbar(limits, width=0.2)`

[関数](http://docs.ggplot2.org/current/stat_function.html)
:   `ggplot(data.frame(x=c(-4, 4)), aes(x)) + stat_function(fun=dnorm, args=c(0, 1), n=200)`

------------------------------------------------------------------------

[回帰曲線](http://docs.ggplot2.org/current/geom_smooth.html)
:   `gp + geom_smooth(method=glm, family=poisson, se=FALSE)`

[切片と傾きで直線を描く](http://docs.ggplot2.org/current/geom_abline.html)
:   `gp + geom_abline(intercept=3, slope=5)`\
    `gp + geom_hline(yintercept=7) + geom_vline(xintercept=11)`

[始点と終点で曲線や矢印を描く](http://docs.ggplot2.org/current/geom_segment.html)
:   `gp + geom_curve(aes(x, y, xend, yend), curvature = -0.2)`\
    `gp + geom_segment(aes(x, y, xend, yend), arrow=arrow())`\
    矢印の調整は [grid::arrow()](https://www.rdocumentation.org/packages/grid/topics/arrow)

[文字列や図形を書き加える](http://docs.ggplot2.org/current/annotate.html)
:   `gp + geom_text(aes(y=y_val+10), label=y_val)`\
    `gp + annotate("text", x=1:4, y=4:1, label=sprintf("x = %d", 1:4))`\
    テーマの `base_family` は引き継がれないので `family=` で指定すべし。\
    数式を使う場合は文字列 `label='italic(N[t])'` のような文字列で渡して `parse=TRUE`。

## プロットの調整

### データによって色やサイズを変える

`colour`: 点や線の色\
`fill`: 塗りつぶしの色\
`size`: 点や線のサイズ <http://docs.ggplot2.org/current/scale_size.html>\
`shape`: 点の形 <http://docs.ggplot2.org/current/scale_shape.html>\
`linetype`: 線 <http://docs.ggplot2.org/current/scale_linetype.html>\
`alpha`: 不透明度 (0が透明、1が不透明) <http://docs.ggplot2.org/current/scale_alpha.html>

`aes()` の中にはデータによって描き分けたい列名を指定し、
外にはデータによらない値を指定する:

    # サイズは常に3、色は species という列のデータによって変える
    geom_point(aes(colour=species), size=3)+
    scale_colour_brewer(palette='Set1')

    # 色は常に赤、サイズは frequency という列の値に比例
    geom_point(aes(size=frequency), colour='red')

<http://docs.ggplot2.org/current/scale_brewer.html>
:   [Colorbrewer](http://colorbrewer2.org/) で定義されているパレットを使う。 色盲対策もできるし、原色も出てこないのでオススメ。\
    cf. [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer/topics/ColorBrewer)\
    `scale_{colour/fill}_brewer(..., type='seq', palette=1, direction=1)`: 離散値\
    `scale_{colour/fill}_distiller(...)`: 連続値

<http://docs.ggplot2.org/current/scale_gradient.html>
:   グラデーションの基準色を指定する。\
    `scale_{colour/fill}_gradient(..., low, high, ...)`: 普通の連続値に\
    `scale_{colour/fill}_gradient2(..., low, mid, high, midpoint=0, ...)`: ある中央値を挟んで上下に分けたいとき\
    `scale_{colour/fill}_gradientn(..., colours, ...)`: 多色のヒートマップなどに e.g., `colours=c('#000000', '#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000')`

<http://docs.ggplot2.org/current/scale_identity.html>
:   データフレームに入ってる値をそのまま使う\
    `scale_{colour/fill/size/shape/linetype/alpha}_identity(..., guide='none')`

<http://docs.ggplot2.org/current/scale_manual.html>
:   値を直に指定する\
    `scale_{colour/fill/size/shape/linetype/alpha}_manual(..., values)`

legend/colourbarのタイトルを変更したい場合は上記関数に `name='New Title'` を指定する。
さらに細かく制御したい場合は
[`guide_legend()`](http://docs.ggplot2.org/current/guide_legend.html) や
[`guide_colourbar()`](http://docs.ggplot2.org/current/guide_colourbar.html)
を引数 `guide` に渡す。


### 内部変数を使う

https://github.com/hadley/ggplot2-book/blob/master/layers.rmd#generated-variables

ヒストグラムや箱ヒゲなどの表示に必要な計算は
ggplot内部で `stat_*()` を通して行われる。
そうした値 (computed variables) の一部は
`..count..` や `..density..` みたいに
ピリオドで囲まれた特殊な名前で参照することができる。


### パネルの分割

年ごとや種ごとに傾向を見たいときなど、データに応じてパネルを分割して並べる。

[`facet_wrap()`](http://docs.ggplot2.org/current/facet_wrap.html)
:   1変数で分割して並べる
    ```r
    facet_wrap(facets, nrow=NULL, ncol=NULL, scales='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE,
               dir='h', strip.position='top')

    ggplot(iris, aes(Petal.Length, Sepal.Length))+geom_point()+
        facet_wrap(~Species, nrow=2)
    ```

[`facet_grid()`](http://docs.ggplot2.org/current/facet_grid.html)
:   2変数で分割して縦横に並べる
    ```r
    facet_grid(facets, margins=FALSE, scales='fixed', space='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE)

    ggplot(iris, aes(Petal.Length, Sepal.Length))+geom_point()+
        facet_grid(. ~ Species)
    ```

    1変数でいい場合は片方をドット `.` で固定できる。

[ファセットラベルの調整](http://docs.ggplot2.org/current/labellers.html)
:   デフォルトでは値だけがfacetラベルに表示されるが、
    変数名を同時に表示するなど細かい調整も可能。


### 軸やタイトルを変更

[軸の区切りを変更したり対数にしたり](http://docs.ggplot2.org/current/scale_continuous.html)
:   `gp + scale_x_continuous(breaks=seq(10, 100, by=10))`\
    `gp + scale_y_log10("Beer consumption")`\
    オプション: name, breaks, labels, na.value, limits, trans, expand\
    デフォルトでは値域よりも少し余裕を持たせてあるが、 `geom_tile()` などでピッタリにしたいときは軸ごとに `expand=c(0, 0)` とか。
    `position='right'`とかで軸の位置を変更できる。
    `sec.axis`オプションで反対側に別の軸を追加できる。

[描画する範囲を指定](http://docs.ggplot2.org/current/coord_cartesian.html)
:   `gp + ylim(0, 42) + xlim("b", "c", "d")`\
    `gp + coord_cartesian(xlim = NULL, ylim = NULL)`\
    前者はデータそのものを切るが、後者はデータを変えずに描画領域だけズームする

[X軸とY軸の比率を固定](http://docs.ggplot2.org/current/coord_fixed.html)
:   `gp + coord_fixed(ratio=1)`

[XY軸の反転](http://docs.ggplot2.org/current/coord_flip.html)
:   `gp + coord_flip()`

[極座標](http://docs.ggplot2.org/current/coord_polar.html)
:   パイチャートも作れるらしい

[座標変換](http://docs.ggplot2.org/current/coord_trans.html)
:   `gp + coord_trans(x='log10', y='sqrt')`\
    表示の座標だけ変更する。
    stat前に適用される `scale_x_*` とかとはちょいと違う。

[軸ラベルとタイトル](http://docs.ggplot2.org/current/labs.html)
:   `gp + labs(x="time", y="weight", title="growth")`\
    `gp + xlab("time") + ylab("weight") + ggtitle("growth")`

## `theme`: 背景やラベルの調整

<http://docs.ggplot2.org/current/ggtheme.html>

### 既成テーマ

`theme_grey(base_size=12, base_family='')`, `theme_gray(...)`
:   灰色背景に白い格子。`ggplot` らしいデフォルト。

`theme_bw(...)`
:   黒枠白背景にうっすら灰色格子

`theme_linedraw(...)`
:   細いけど濃い色の `panel.grid`

`theme_light(...)`
:   それを薄くした感じ

`theme_minimal(...)`
:   外枠なしの `theme_bw`

`theme_classic(...)`
:   xy軸がL字に描かれているだけで枠もグリッドも無し

`theme_void(...)`
:   完全に枠なし

これらをカッコ無しでコンソールに打ち込むと、
下記の各エレメントの設定方法やデフォルト値を知ることができる。

引数として `base_family` に"Helvetica Neue"などのフォントを指定できる。
Macなら"HiraKakuProN-W3"を指定すれば日本語でも文字化けしなくなるはず。
(`grDevices::quartzFonts()` などを `.Rprofile` に書いてフォントを登録するには
[config]({{< relref "config.md#rprofile" >}}) を参照)

ほかにもいろんなテーマが
[ggthemes](https://cran.r-project.org/web/packages/ggthemes/vignettes/ggthemes.html)
というパッケージで提供されている。

{{%div class="note"%}}
この設定はテーマを構成する `axis.text` などには引き継がれるが
`geom_text()` などプロット内部の要素には引き継がれない。
{{%/div%}}


### 設定項目

<http://docs.ggplot2.org/current/theme.html>

<http://docs.ggplot2.org/dev/vignettes/themes.html>

`theme()` 関数に項目と値を指定したものを、
ほかのレイヤーと同じようにどんどん足しながら変更していく。
テーマは直線、長方形、文字の3種類のエレメントからなり、
それらの性質を変更する場合は `element_***()` を介して行う

```r
## ベースとなるテーマを先に適用してから
gp = gp + theme_bw(base_family='HiraKakuProN-W3')
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
: `axis.title`: (`element_text`; inherits from `text`)\\
  &emsp;`__.x`, `__.x.top`, `__.y`, `__.y.right`
: `axis.text`: (`element_text`; inherits from `text`)\\
  &emsp;`__.x`, `__.x.top`, `__.y`, `__.y.right`
: `axis.ticks`: (`element_line`; inherits from `line`)\\
  &emsp;`__.x`, `__.y`
: `axis.ticks.length`: (`unit`)
: `axis.line`: (`element_line`; inherits from `line`)\\
  &emsp;`__.x`, `__.y`

凡例
: `legend.background`: (`element_rect`; inherits from `rect`)
: `legend.margin`: (`margin`)
: `legend.spacing`:(`unit`)\\
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
: `panel.spacing`: (`unit`; `facet_*` の間隔)\\
  &emsp;`__.x`, `__.y`
: `panel.grid`: (`element_line`; inherits from `line`)
: `panel.grid.major`: (`element_line`; inherits from `panel.grid`)\\
  &emsp;`__.x`, `__.y`
: `panel.grid.minor`: (`element_line`; inherits from `panel.grid`)\\
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
: `strip.text`: (`element_text`; inherits from `text`)\\
  &emsp;`__.x`, `__.y`

その他
: `complete`: 部分的な変更か、完全なテーマか (`FALSE`)
: `validate`: 毎回チェックするか (`TRUE`)


### エレメント

`element_rect(fill, colour, size, linetype, inherit.blank)` --- 長方形
:   `fill`: 塗りつぶしの色\
    `colour`: 枠の色

`element_line(colour, size, linetype, lineend, arrow, inherit.blank)` --- 直線

`element_text(family, face, colour, size, hjust, vjust, angle, lineheight, margin)` --- 文字
:   `family`: フォントファミリー。 空なら `theme_bw(base_family=...)` などの指定を継承。\
    `face`: (`'plain'`, `'italic'`, `'bold'`, `'bold.italic'`)\
    `hjust`, `vjust`: 水平位置と垂直位置の寄せ方をそれぞれ `[0, 1]` の実数で。\
    `angle`: 角度 `[0, 360]`\
    `margin`: スペース調整を関数 `margin(top, right, bottom, left)` 越しに。

`element_blank()` --- 空
:   消したい要素にはこれを指定する
    ```r
    gp = gp + theme(axis.title=element_blank())
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

[線の種類](http://www.cookbook-r.com/Graphs/Shapes_and_line_types/) --- `linetype=`
:   `0`: `'blank'`\
    `1`: `'solid'`\
    `2`: `'dashed'`\
    `3`: `'dotted'`\
    `4`: `'dotdash'`\
    `5`: `'longdash'`\
    `6`: `'twodash'`

## ファイルに書き出す

```r
ggsave(filename = default_name(plot), plot = last_plot(),
    device = default_device(filename), path = NULL, scale = 1,
    width = par("din")[1], height = par("din")[2],
    units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, ...)
```

サイズを変えるには scale, width, height, units, dpi をいじる。
デフォルトはそれぞれ 1, 7, 7, 'in', 300 で、1辺177.8mmの正方形に2100x2100ピクセル並べるイメージ。
A4 PDFにするなら `ggsave('fig1.pdf', gp, width=8.27, height=11.7)` とか
`ggsave('fig1.pdf', gp, width=210, height=297, units='mm')` とか。
画像形式はファイル名の拡張子から自動的に判別される。

{{%div class="note"%}}
`width` と `height` が省略されると `par("din")` が呼び出される。
Rscript から実行してるのに白いQuartzウィンドウが出現したり、
空の `Rplots.pdf` ができたりするのはこのせい。
{{%/div%}}

MacのQuartzで調整して、そこで見たままを保存したい場合は

```r
dev.off()
quartz(width=9.9, height=7)
print(gp)
quartz.save('plot.png')
```

## Extensions

ggplotを拡張するための仕組みがversion 2.0から正式に導入され、
ユーザーが独自の stats や geom を作って登録することが容易になった。

<http://www.ggplot2-exts.org/>

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
ラベル(e.g., A, B, C)をオプションで簡単に付けられるのが良い。

{{%div class="warning"%}}
`library(cowplot)` せずに `loadNamespace('cowplot')`
して常に名前空間つきで呼んだほうが安全。

-   `cowplot::ggsave()`が定義されてて`ggplot2`と衝突するのが気持ち悪い。
-   `.onAttach()` で基本テーマが勝手に変更されるのが気持ち悪い。
    `theme_cowplot()` そのものはシンプルで良くできているが、
    カスタマイズのベースとするには扱いにくい
    (`rect`に変なデフォルト値が設定されてしまう、など)。
{{%/div%}}

`cowplot::plot_grid()`
:   `facet_wrap()`のように、ざっと並べるのに便利。

    ```r
    corplot::plot_grid(..., plotlist=NULL,
        align=c('none', 'h', 'v', 'hv'),
        nrow=NULL, ncol=NULL,
        scale=1, rel_widths=1, rel_heights=1,
        labels=NULL, label_size=14,
        hjust=-0.5, vjust=1.5)
    ```
    さらに細やかな制御をしたいときは以下の関数を個別に重ねていく。

`cowplot::ggdraw(plot=NULL, xlim=c(0, 1), ylim=c(0, 1))`
:  これの後ろに `+` 演算子で `draw_***()` を足していく。

`draw_figure_label()`
`draw_grob()`
`draw_label()`
`draw_line()`
`draw_plot()`
`draw_plot_label()`
`draw_text()`


### `GGally`

<http://cran.r-project.org/web/packages/GGally/>

`graphics::pairs()` のような plot matrix を `ggplot2` で作るためのパッケージ。
そもそもデータをザッと俯瞰するためのものであり、
レイヤーを重ねて変更を加えていくようなものでもないので、
`ggplot2` である恩恵はそこまでない。
でもまあ `pairs(iris)` よりは `ggpairs(iris, colour=Species)`
のほうが見やすいのは確か。

```r
GGally::ggpairs(data, columns=1:ncol(data), title='',
    upper=list(), lower=list(), diag=list(), params=NULL, ...,
    axisLabels='internal', legends=FALSE, verbose=FALSE)
```

## 関連書籍

<a href="https://www.amazon.co.jp/R%E3%82%B0%E3%83%A9%E3%83%95%E3%82%A3%E3%83%83%E3%82%AF%E3%82%B9%E3%82%AF%E3%83%83%E3%82%AF%E3%83%96%E3%83%83%E3%82%AF-_ggplot2%E3%81%AB%E3%82%88%E3%82%8B%E3%82%B0%E3%83%A9%E3%83%95%E4%BD%9C%E6%88%90%E3%81%AE%E3%83%AC%E3%82%B7%E3%83%94%E9%9B%86-Winston-Chang/dp/4873116538/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=911c7022439b7b328c58e0168e4c66e7" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873116538&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4873116538" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/%E3%82%B0%E3%83%A9%E3%83%95%E3%82%A3%E3%83%83%E3%82%AF%E3%82%B9%E3%81%AE%E3%81%9F%E3%82%81%E3%81%AER%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%9F%E3%83%B3%E3%82%B0-H-%E3%82%A6%E3%82%A3%E3%83%83%E3%82%AB%E3%83%A0/dp/4621061356/ref=as_li_ss_il?ie=UTF8&linkCode=li3&tag=heavywatal-22&linkId=5ef26b015fb6a8298b410cf9feafb728" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621061356&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li3&o=9&a=4621061356" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
