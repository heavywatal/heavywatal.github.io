+++
date = 2016-05-23T17:51:48+09:00
tags = ["math", "graph"]
title = "数理形態学"
subtitle = "Mathematical morphology"

[menu.main]
  parent = "bio"
+++

2D/3D Cellular Automaton上の個体・細胞の分布を評価したい。
そのためには白黒の二値画像処理の手法が結構使える。

図形 *X*
: それぞれのノード(画素)の在・不在情報の集合。
  $x \in X$

構造要素 (Structuring Element: SE)
: さまざまな処理を施すために用いられる単位図形のようなもの。
  例えば、原点とそのムーア近傍。
  $b \in B$

voxel
: 3D空間における単位。2Dでいうpixel。

## 基本処理

### Translation 平行移動

<div>\[
X_b = \{x + b \mid x \in X\}
\]</div>

### Dilation 膨張

<div>\[
X \oplus B = \bigcup_{b \in B} X_b
\]</div>

*X* と *B* のMinkowski和。
*X* を *B* の範囲でずらしながらunionを取ったもの。
国土を*X* 、半径12海里の円をSEとした、領空みたいなイメージ。

### Erosion 浸食

<div>\[
X \ominus B = \bigcap_{b \in B} X_b
\]</div>

*X* と *B* のMinkowski差。
*X* を *B* の範囲でずらしながらintersectを取ったもの。
SEを消しゴムとして *X* の外周上を走らせ、削るイメージ。

### Opening

<div>\[
X \circ B = (X \ominus B) \oplus B
\]</div>

浸食してから膨張する。
*X* からハミ出ないようにSEを滑らせた軌跡に相当する。
トゲの先端や *X* 外部のチリなど、SEより小さい構造が削られて小さくなる。
特定の形を持ったSEを使えば、それを含む領域だけを抽出するのにも使える。

元画像との差分 $X - (X \circ B)$ は **Top Hat** と呼ばれ、
トゲの先っちょや背景のノイズ成分が得られる。

### Closing

<div>\[
X \bullet B = (X \oplus B) \ominus B
\]</div>

膨張してから浸食する。
*X* の外部をOpeningすることと同義。
*X* 内部のヒビやチリなど、SEより小さい構造が塗りつぶされ、大きくなる。

元画像との差分 $(X \bullet B) - X$ は **Black Hat** と呼ばれ、
*X* 内のヒビやトゲの根元らへんが得られる。

## 応用

### Pattern Spectrum, サイズ分布

小さいSEから順に大きくしながら
Openingで削れた部分の面積を記録していく。
元画像の面積で割ったものはサイズ密度関数(size density function)と呼ばれる。
細かいギザギザを含む図形ほど小さいSEで削れる成分が多い。
要約統計量としてはモーメントやエントロピーが使える。

### Morphological gradient

<div>\[
(X \oplus B) - (X \ominus B)
\]</div>

dilationとerosionの差。
エッジ検出法のひとつ。
X上の境界が欲しい場合は$X - (X \ominus B)$。
背景側の境界が欲しい場合は$(X \oplus B) - X$。


## ノイズ除去

平滑化フィルタ
: SEを端から端まで動かしつつ、その中に含まれる画素の平均値を中央画素に適用していく。
  Gaussian filterのように、遠いものほど軽くなるように重み付けをする場合もある。
  いずれにせよ、エッジがボヤけてしまうのが問題。

Median filter
: 平均値ではなく中央値で置き換える。
  エッジは保存されるが、ソートを伴うので計算量は多め。


## ライブラリ

画像処理を施す

[scikit-image](http://scikit-image.org/)
: Pythonモジュール。
  [scipy.ndimage](http://docs.scipy.org/doc/scipy/reference/tutorial/ndimage.html)
  を更に拡張したもの。
  [numpy.array](http://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html)
  を使って表現されるので汎用関数の適用も容易。

[OpenCV (Open Source Computer Vision)](http://opencv.org/)
: C++、Pythonなど。信頼と実績があるらしく、書籍やネット上の情報も多い。

[CImg](http://cimg.eu/)
: C++。ヘッダひとつincludeするだけ。
  ドキュメントも良さげ。

[imager](http://dahtah.github.io/imager/)
: R。新しめでドキュメントも充実。内部でCImgを利用。
  CairoやX11に依存しているので、
  Rもそれらしくビルドされてる必要がある。

[mmand](https://github.com/jonclayden/mmand)
: R。READMEは良さげ。内部はRcpp。
  ほとんど他のライブラリに依存していないのでインストールしやすい。

[Morpho](https://github.com/zarquon42b/Morpho)
: R。ドキュメント不足。

----

- [非線形画像・信号処理 (モルフォロジの基礎と応用)](https://www.amazon.co.jp/dp/4621082949?&linkCode=ll1&tag=heavywatal-22&linkId=8d2ef0480c9752ec4c45ea04893e0fa2)
- [3次元画像処理](https://www.amazon.co.jp/dp/4862460844?&linkCode=ll1&tag=heavywatal-22&linkId=c4f0058e3e42e8d1ba850b1c9276e02d)
