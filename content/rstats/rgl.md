+++
title = 'rgl'
subtitle = "3Dグラフ描画"
tags = ["r", "graph"]
[menu.main]
  parent = "rstats"
  weight = 1
+++

-   [Project Homepage](https://dmurdoch.github.io/rgl/)
-   <https://cran.r-project.org/web/packages/rgl/>

長らく2種類の書き方が混在していたが、
バージョン 1.0 から[`rgl.*()` 関数の使用が非推奨](https://dmurdoch.github.io/rgl/articles/deprecation.html)となり一本化された。
ドキュメントも一新されて使いやすくなってきた。

## プロット

### Primitive shapes

<https://dmurdoch.github.io/rgl/reference/primitives.html>

`rgl::points3d(x, y = NULL, z = NULL, ...)`
:   散布図。

`rgl::lines3d(x, y = NULL, z = NULL, ...)`
:   折れ線

`rgl::segments3d(x, y = NULL, z = NULL, ...)`
:   線分

`rgl::triangles3d(x, y = NULL, z = NULL, ...)`
:   3点を結ぶ面

`rgl::quads3d(x, y = NULL, z = NULL, ...)`
:   4点を結ぶ三角形2つ

### Other shapes

`rgl::spheres3d(x, y = NULL, z = NULL, radius = 1, fastTransparency = TRUE, ...)`
:   球体。

`rgl::surface3d(x, y = NULL, z = NULL, ..., normal_x = NULL, normal_y = NULL, normal_z = NULL, texture_s = NULL, texture_t = NULL, flip = FALSE)`
:   地形図のような局面

`rgl::plot3d(x, ...)`
:   `type =`引数で上記の様々な形を描ける高次関数。
    ホントはあんまり使いたくないけど、
    `xlim`, `ylim`, `zlim` オプションを受け付ける関数がこれしかないようなので、
    境界を指定しつつ球体を描きたい場合は `spheres3d()` ではなく
    `plot3d(type = "s")` を使うしかないっぽい。


## 背景や軸などの調整

`rgl::title3d(main, sub, xlab, ylab, zlab, line, level, floating, ...)`
:   これを使うとmainとsubも視点によって動いてしまう。
    `bgplot3d({plot.new(); title("main")})` なら固定背景に書ける。

`rgl::mtext3d(text, edge, at = NULL, line = 0, level = 0, floating = FALSE, pos = NA, ...)`

`rgl::bg3d(color, sphere = FALSE, back = "lines", fogtype = "none", fogScale = 1, col, ...)`

`rgl::light3d(theta = 0, phi = 15, x = NULL, y = NULL, z = NULL, viewpoint.rel = TRUE, ambient = "#FFFFFF", diffuse = "#FFFFFF", specular = "#FFFFFF")`

`rgl::par3d(..., no.readonly = FALSE, dev = cur3d(), subscene = currentSubscene3d(dev))`

[`rgl::material3d(..., id = NULL)`](https://dmurdoch.github.io/rgl/reference/material.html)
: プロットに渡せるオプション(`color`など)はここで確認

### 軸

`rgl::axis3d(edge, at = NULL, labels = TRUE, tick = TRUE, line = 0, pos = NULL, nticks = 5, ...)`
:   `xyz` と `+-` の組み合わせで軸1本を指定して描く。
    `x` は　`x--` と等価。

`rgl::box3d(...)`
:   12辺の箱を描く。

`rgl::bbox3d(xat = NULL, yat = NULL, zat = NULL, xunit = "pretty", yunit = "pretty", zunit = "pretty", expand = 1.03, draw_front = FALSE, xlab = NULL, ylab = NULL, zlab = NULL, xlen = 5, ylen = 5, zlen = 5, marklen = 15, marklen.rel = TRUE, ...)`
:   手前の辺が自動で消えるような箱を描く。

`rgl::axes3d(edges = "bbox", labels = TRUE, tick = TRUE, nticks = 5, box = FALSE, expand = 1.03, ...)`
:   上記の3つをまとめる関数。分かりにくいので使わないほうがいい。
    `edges = "bbox"` の場合 `tick = FALSE` は効かないので `xlen = 0, ylen = 0, zlen = 0` とする必要がある。

`rgl::view3d(theta = 0, phi = 15, fov = 60, zoom = 1, scale = par3d("scale"), interactive = TRUE, userMatrix)`
:   `theta`: 0のとき正面がxy平面。観察者が地球の公転と同じ方向に動くのが正。\
    `phi` [-90, 90]: 0のとき視点が水平面(xz平面)上。観察者が上に動くのが正。\
    `fov` [0, 179]: 0のとき無限遠から見たような平行投影。


### 複数の図をまとめる

```r
# レイアウトを指定
mfrow3d(nr, nc, byrow = TRUE, parent = NA, sharedMouse = FALSE, ...)
layout3d(mat, widths, heights, parent = NA, sharedMouse = FALSE, ...)

# 次のsubsceneに移動
next3d(current = NA, clear = TRUE, reuse = TRUE)
```

これらはなぜかグローバルスコープでしか動作しない。
つまり、関数やループ内に入れるとサイズなどがうまく反映されない。


## 出力

### デバイスの起動と終了

`rgl::open3d(..., params = get3dDefaults(), useNULL = rgl.useNULL(), silent = FALSE)`
: 明示的に新しいデバイスを開く。
  何も無い状態で`plot3d()`などが呼ばれたら勝手に開かれる。
  サイズ指定は`windowRect = c(0, 0, 600, 600)`のような引数で。

`rgl::close3d(dev = cur3d(), silent = TRUE)`
: デバイスを閉じる。

`rgl::clear3d(type = c("shapes", "bboxdeco", "material"), defaults = getr3dDefaults(), subscene = 0)`

### Display

<https://dmurdoch.github.io/rgl/dev/articles/rgl.html#default-display>

デフォルトでは独立のウィンドウ(XQuartz, X11など)が立ち上がる。

`options(rgl.useNULL = TRUE)`
: 独立ウィンドウが開くのを抑制。

`rglwidget()`
: WebGLに変換してRStudio, [VSCode]({{< relref "vscode.md" >}}), ウェブブラウザなど出力。

`options(rgl.printRglwidget = TRUE)`
: 自動的に `rglwidget()` を呼ぶ。使わないほうが無難。


### ファイルに書き出す

`rgl::scene3d(minimal = TRUE)`
: rglネイティブな形での全構成要素リスト。

`rgl::snapshot3d(filename, fmt = "png", top = TRUE, ..., scene, width, height, webshot)`
: PNGのみ。
  `top = FALSE`にしてはダメ。謎。

`rbl.postscript(filename, fmt = "eps", drawText = TRUE)`
:   ps, eps, tex, pdf, svg をサポート。
    透過や `bgplot3d` は反映されないらしいので注意。

`rgl::writeWebGL()`
:   deprecatedだから代わりに `rglwidget()` を使えとのことだがそちらにファイル書き出し機能は無い。


### HTMLに埋め込む

[Quarto/RMarkdown/knitr]({{< relref "knitr.md" >}})でコードを書く。
パッケージを読み込み、hookを設定しておく(`rgl::setupKnitr()` を使う手もある):

````markdown
```{r library}
options(rgl.useNULL = TRUE)
library(rgl)
knitr::knit_hooks$set(webgl = rgl::hook_webgl)
knitr::knit_hooks$set(rgl = rgl::hook_rgl)
```
````

WebGLを出力したい場合は `webgl = TRUE`、
PNG静止画を出力したい場合は `rgl = TRUE`:

````markdown
```{r plot, webgl = TRUE}
rgl::box3d()
rgl::title3d("main", "sub", "x", "y", "z")
```
````

`rgl::rglwidget()` を明示的に呼ぶならhookは不要。
複数描画したいときは `htmltools::tagList()` に詰める手もある:

````markdown
```{r widget}
purrr::map(seq_len(3), ~{
  rgl::box3d()
  rgl::rglwidget(width = 200, height = 200)
}) |> htmltools::tagList()
```
````

WebGLへの変換は出力先がHTMLであることを条件にしているらしく、
`rmarkdown::render()` なら上記コードで問題ないが `knitr::knit()` はダメ。
強制的に変換する手段はあるのかな...?


### アニメーション

`rgl::spin3d(axis = c(0, 0, 1), rpm = 5, dev = cur3d(), subscene)`

`rgl::par3dinterp(times = NULL, userMatrix, scale, zoom, FOV, method, extrapolate)`

`rgl::play3d(f, duration = Inf, dev = cur3d(), ..., startTime = 0)`

`rgl::movie3d(f, duration, dev = cur3d(), ..., fps = 10, movie = "movie", frames = movie, dir = tempdir(), covert, clean, verbose, top, type, startTime, webshot)`

```r
rgl::view3d(-25, 15, 40)
rgl.bringtotop()
.anime = rgl::spin3d(axis = c(0, 1, 0), rpm = 15)
rgl::play3d(.anime)
rgl::movie3d(.anime, duration = 4, fps = 16, movie = "basename", dir = "~/tmp")
```


## ほかに3Dグラフを描けそうな手段:

- [plotly](https://plotly.com/r/3d-charts/)
- [threejs](https://bwlewis.github.io/rthreejs/)
- [VisPy](https://vispy.org/)
