+++
title = 'rgl'
subtitle = "3Dグラフ描画"
tags = ["r", "graph"]
[menu.main]
  parent = "rstats"
+++

3D visualization device system (OpenGL).

-   [Project Homepage](http://rgl.neoscientists.org/)
-   [R-tips](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/57.html)
-   <https://cran.r-project.org/web/packages/rgl/vignettes/rgl.html>
-   <http://cran.r-project.org/web/packages/rgl/>
-   <http://www.rdocumentation.org/packages/rgl>

## デバイスの起動と終了

```r
rgl::open3d()  # open new device
rgl.close()    # close current device
rgl.quit()     # shutdown rgl device system
```

`rgl::clear3d(type=c('shapes', 'bboxdeco', 'material'), defaults, subscene=0)`

## プロット

`rgl::points3d()`
:   散布図。

`rgl::spheres3d()`
:   球。

`rgl::lines3d()`
:   折れ線

`rgl::segments3d()`
:   線分

`rgl::triangles3d()`
:   3点を結ぶ面

`rgl::quads3d()`
:   4点を結ぶ三角形2つ

`rgl::surface3d()`, `rgl::terrain3d()`
:   地形図のような局面

## 背景や軸などの調整

`rgl::title3d(main, sub, xlab, ylab, zlab, line=NA, ...)`
:   これを使うとmainとsubも視点によって動いてしまう。
    `bgplot3d({plot.new(); title('main')})` なら固定背景に書ける。

`rgl::mtext3d(text, edge, line=0, at=NULL, pos=NA, ...)`

`rgl::bg3d()`

`rgl::light3d()`

`rgl::par3d()`

`rgl::material3d()`

### 軸

`rgl::axis3d(edge, at=NULL, labels=TRUE, tick=TRUe, line=0, pos=NULL, nticks=5, ...)`
:   `xyz` と `+-` の組み合わせで軸1本を指定して描く。
    `x` は　`x--` と等価。

`rgl::box3d(...)`
:   12辺の箱を描く。

`rgl::bbox3d(xat=NULL, yat, zat, xunit='pretty', yunit, zunit, expand=1.03, draw_front=FALSE)`
:   手前の辺が自動で消えるような箱を描く。

`rgl::axes3d(edges='bbox', labels=TRUE, tick=TRUE, nticks=5, box=FALSE, expand=1.03, ...)`
:   上記の3つをまとめる関数。分かりにくいので使わないほうがいい。
    `edges='bbox'` の場合 `tick=FALSE` は効かないので `xlen=0, ylen=0, zlen=0` とする必要がある。

`rgl::view3d(theta=0, phi=15, fov=60, zoom=1, scale=par3d("scale"), interactive=TRUE, userMatrix)`
:   `theta`: 0のとき正面がxy平面。観察者が地球の公転と同じ方向に動くのが正。\
    `phi` [-90, 90]: 0のとき視点が水平面(xz平面)上。観察者が上に動くのが正。\
    `fov` [0, 179]: 0のとき無限遠から見たような平行投影。

## ファイルに書き出す

`rgl::snapshot3d(...)`
:   PNGのみ

`rbl.postscript(filename, fmt='eps', drawText=TRUE)`
:   ps, eps, tex, pdf, svg をサポート。
    透過や `bgplot3d` は反映されないらしいので注意。

`rgl::writeWebGL(dir='webGL', filename, template, prefix, snapshot, commonParts, reuse, font, width, height)`
:   ディレクトリ構造無しの単発HTMLでいい場合は
    `writeWebGL('.', 'rgl.html')` のように指定する。

## アニメーション

`rgl::spin3d(axis=c(0, 0, 1), rpm=5)`

`rgl::par3dinterp(times=NULL, userMatrix, scale, zoom, FOV, method, extrapolate)`

`rgl::play3d(f, duration=Inf, ...)`

`rgl::movie3d(f, duration, ..., fps=10, movie="movie", frames=movie, dir=tempdir(), ...)`

```r
## 角度をセット
rgl::view3d(-25, 15, 40)

## 最前面に持ってくる
rgl.bringtotop()

## アニメーション関数を作る
.anime = rgl::spin3d(axis=c(0, 1, 0), rpm=15)

## X11で再生
rgl::play3d(.anime)

## GIFアニメとして保存
rgl::movie3d(.anime, duration=4, fps=16, movie="basename", dir="~/tmp")
```