+++
title = 'ggplot2'
subtitle = "きれいなグラフを簡単に合理的に"
[menu.main]
  parent = "rstats"
+++

-   <http://ggplot2.org/>
-   <http://docs.ggplot2.org/>
-   <http://www.cookbook-r.com/Graphs/>
-   <http://www.rdocumentation.org/packages/ggplot2>

R で以下のコマンドを実行してインストール

```r
install.packages("ggplot2")
```

## 基本的な使い方

R に入ってるお馴染みサンプルデータ `iris` を使って

```r
## ライブラリの読み込み
library(ggplot2)

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
`ggplot()` に渡すデータを整形するには [tidyr]({{< relref "tidyr.md" >}}) が便利。
{{%/div%}}

## プロット

[散布図](http://docs.ggplot2.org/current/geom_point.html)
:   `gp + geom_point(size=2, alpha=0.3)`

[折れ線グラフ](http://docs.ggplot2.org/current/geom_line.html)
:   `gp + geom_line(size=2, linetype="dashed")`\
    `gp + geom_path()`\
    前者はx軸の小さい順に結び、後者はデータ順に結ぶ

[ヒストグラム、密度曲線](http://docs.ggplot2.org/current/geom_histogram.html)
:   `gp + geom_histogram(fill=..count..)`\
    `gp + geom_density(alpha = 0.2)`

[棒グラフ](http://docs.ggplot2.org/current/geom_bar.html)
:   `gp + geom_bar(stat='identity')`\
    stat を指定しないとヒストグラムになってしまう。 `position='dodge'` にすると横並び (デフォルト: `'stack'`)。

[箱ひげ図](http://docs.ggplot2.org/current/geom_boxplot.html)
:   `gp + geom_boxplot()`\
    `gp + geom_violin()`

[ヒートマップ](http://docs.ggplot2.org/current/geom_tile.html)
:   `gp + geom_tile(aes(fill=z))`\
    `gp + geom_raster(aes(fill=z))`\
    後者は各タイルの大きさがすべて同じ場合の特殊ケースで、高速。

[エラーバー](http://docs.ggplot2.org/current/geom_errorbar.html)
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
    矢印の調整は [grid::arrow()](http://www.inside-r.org/r-doc/grid/arrow)

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
    cf. [RColorBrewer](http://www.inside-r.org/packages/cran/RColorBrewer/docs/ColorBrewer)\
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

### パネルの分割

年ごとや種ごとに傾向を見たいときなど、データに応じてパネルを分割して並べる。

<http://docs.ggplot2.org/current/facet_wrap.html>
:   1変数で分割して並べる

    ```r
    facet_wrap(facets, nrow=NULL, ncol=NULL, scales='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE)

    ggplot(iris, aes(Petal.Length, Sepal.Length))+geom_point()+
        facet_wrap(~Species, nrow=2)
    ```

<http://docs.ggplot2.org/current/facet_grid.html>
:   2変数で分割して縦横に並べる

    ```r
    facet_grid(facets, margins=FALSE, scales='fixed', space='fixed',
               shrink=TRUE, labeller='label_value',
               as.table=TRUE, switch=NULL, drop=TRUE)

    ggplot(iris, aes(Petal.Length, Sepal.Length))+geom_point()+
        facet_grid(. ~ Species)
    ```

    1変数でいい場合は片方をドット `.` で固定できる。

<http://docs.ggplot2.org/current/labellers.html>
:   デフォルトでは値だけがfacetラベルに表示されるが、
    変数名を同時に表示するなど細かい調整も可能。

### 軸やタイトルを変更

[軸の区切りを変更したり対数にしたり](http://docs.ggplot2.org/current/scale_continuous.html)
:   `gp + scale_x_continuous(breaks=seq(10, 100, by=10))`\
    `gp + scale_y_log10("Beer consumption")`\
    オプション: name, breaks, labels, na.value, limits, trans, expand\
    デフォルトでは値域よりも少し余裕を持たせてあるが、 `geom_tile()` などでピッタリにしたいときは軸ごとに `expand=c(0, 0)` とか。

[描画する範囲を指定](http://docs.ggplot2.org/current/coord_cartesian.html)
:   `gp + ylim(0, 42) + xlim("b", "c", "d")`\
    `gp + coord_cartesian(xlim = NULL, ylim = NULL)`\
    前者はデータそのものを切るが、後者はデータを変えずに描画領域だけズームする

[XY軸の反転](http://docs.ggplot2.org/current/coord_flip.html)
:   `gp + coord_flip()`

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
この設定はテーマを構成する `axis.text` などには引き継がれるが
`geom_text()` などプロット内部の要素には引き継がれない。

{{%div class="note"%}}
[config]({{< relref "config.md" >}})

`.Rprofile` で `grDevices::quartzFonts()` を登録しておけば
その名前でも指定できる。
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

`line`: (`element_line`)\
`rect`: (`element_rect`)\
`text`: (`element_text`)\
`title`: (`element_text`; inherits from `text`)

軸タイトル、軸ラベル、目盛

`axis.title`: (`element_text`; inherits from `text`)\
`axis.title.x`: (`element_text`; inherits from `axis.title`)\
`axis.title.y`: (`element_text`; inherits from `axis.title`)\
`axis.text`: (`element_text`; inherits from `text`)\
`axis.text.x`: (`element_text`; inherits from `axis.text`)\
`axis.text.y`: (`element_text`; inherits from `axis.text`)\
`axis.ticks`: (`element_line`; inherits from `line`)\
`axis.ticks.x`: (`element_line`; inherits from `axis.ticks`)\
`axis.ticks.y`: (`element_line`; inherits from `axis.ticks`)\
`axis.ticks.length`: (`unit`)\
`axis.ticks.margin`: (`unit`)\
`axis.line`: (`element_line`; inherits from `line`)\
`axis.line.x`: (`element_line`; inherits from `axis.line`)\
`axis.line.y`: (`element_line`; inherits from `axis.line`)

凡例

`legend.background`: (`element_rect`; inherits from `rect`)\
`legend.margin`: (`unit`)\
`legend.key`: (`element_rect`; inherits from `rect`)\
`legend.key.size`: (`unit`)\
`legend.key.height`: (`unit`; inherits from `legend.key.size`)\
`legend.key.width`: (`unit`; inherits from `legend.key.size`)\
`legend.text`: (`element_text`; inherits from `text`)\
`legend.text.align`: (number from `0` (left) to `1` (right))\
`legend.title`: (`element_text`; inherits from `title`)\
`legend.title.align`: (number from `0` (left) to `1` (right))\
`legend.position`: (`"left"`, `"right"`, `"bottom"`, `"top"`, `"none"`; プロット領域内での位置を `c(0, 1)` のような数値で)\
`legend.direction`: (`"horizontal"` or `"vertical"`)\
`legend.justification`: (`"center"` or `c(0, 1)` のような数値でアンカー位置を指定)\
`legend.box`: (`"horizontal"` or `"vertical"`)\
`legend.box.just`: (`"top"`, `"bottom"`, `"left"`, or `"right"`)

プロット領域の背景、余白、格子

`panel.background`: (`element_rect`; inherits from `rect`)\
`panel.border`: (`element_rect`; inherits from `rect`; should be used with `fill=NA`)\
`panel.margin`: (`unit`; `facet_*` の間隔)\
`panel.grid`: (`element_line`; inherits from `line`)\
`panel.grid.major`: (`element_line`; inherits from `panel.grid`)\
`panel.grid.minor`: (`element_line`; inherits from `panel.grid`)\
`panel.grid.major.x`: (`element_line`; inherits from `panel.grid.major`)\
`panel.grid.major.y`: (`element_line`; inherits from `panel.grid.major`)\
`panel.grid.minor.x`: (`element_line`; inherits from `panel.grid.minor`)\
`panel.grid.minor.y`: (`element_line`; inherits from `panel.grid.minor`)

全体の背景、タイトル、余白

`plot.background`: (`element_rect`; inherits from `rect`)\
`plot.title`: (`element_text`; inherits from `title`)\
`plot.margin`: (`unit` with the sizes of the top, right, bottom, and left margins)

`facet` したときのラベル

`strip.background`: (`element_rect`; inherits from `rect`)\
`strip.text`: (`element_text`; inherits from `text`)\
`strip.text.x`: (`element_text`; inherits from `strip.text`)\
`strip.text.y`: (`element_text`; inherits from `strip.text`)

### エレメント

`element_line(colour, size, linetype, lineend)` --- 直線

`element_rect(fill, colour, size, linetype)` --- 長方形
:   `fill`: 塗りつぶしの色\
    `colour`: 枠の色

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

<http://ggplot2-exts.github.io/>

### `gridExtra`

<https://github.com/baptiste/gridextra/wiki>

データによって自動的にパネルを分割するには
`facet_grid()` や `facet_wrap()` を使えばよいが、
関係ない複数の図を1枚に描きたい場合は `grid` や `gtable` の機能を使う必要がある。
`gridExtra` はそのへんの操作を手軽にできるようにしてくれるパッケージ。

```r
grob = gridExtra::arrangeGrob(p1, p2, nrow=2, ncol=1, bottom='Time')
grid.newpage()
grid.draw(grob)
```

"grob" は "grid graphical object" の略。
ggplotオブジェクトと同じように `ggsave()` に渡して保存可能。

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

## 書籍

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4873116538/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51S2-F8zkRL._SX180_.jpg" alt="Rグラフィックスクックブック ―ggplot2によるグラフ作成のレシピ集" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4621061356/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41r4XY9%2B-PL._SX160_.jpg" alt="グラフィックスのためのRプログラミング" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4320019059/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/41QhVxe5ajL._SX160_.jpg" alt="Rグラフィックス ―Rで思いどおりのグラフを作図するために―" /></a>
