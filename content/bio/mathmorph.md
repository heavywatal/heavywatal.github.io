+++
date = "2016-05-23T17:51:48+09:00"
tags = ["math", "graph"]
title = "数理形態学"
subtitle = "Mathematical morphology"

[menu.main]
  parent = "bio"
+++

図形 *X*
: binaryでいう1の画素の集合。
  $x \in X$

構造要素 (Structuring Element: SE)
: 原点と近傍の集合。
  $b \in B$

voxel
: 3D空間における単位。2Dでいうpixel。

## 基本処理

### Translation 平行移動

<div>$$\begin{split}
X_b = \{x + b \mid x \in X\}
\end{split}$$</div>

### Dilation 膨張

<div>$$\begin{split}
X \oplus B = \bigcup_{b \in B} X_b
\end{split}$$</div>

*X* と *B* のMinkowski和。
*X* を *B* の範囲でずらしながらunionを取ったもの。
国土を*X* 、半径12海里の円をSEとした、領空みたいなイメージ。

### Erosion 浸食

<div>$$\begin{split}
X \ominus B = \bigcap_{b \in B} X_b
\end{split}$$</div>

*X* と *B* のMinkowski差。
*X* を *B* の範囲でずらしながらintersectを取ったもの。
SEを消しゴムとして *X* の外周上を走らせ、削るイメージ。

### Opening

<div>$$\begin{split}
X \circ B = (X \ominus B) \oplus B
\end{split}$$</div>

浸食してから膨張する。
*X* からハミ出ないようにSEを滑らせた軌跡に相当する。
トゲの先端や *X* 外部のチリなど、SEより小さい構造が削られて小さくなる。
特定の形を持ったSEを使えば、それを含む領域だけを抽出するのにも使える。

### Closing

<div>$$\begin{split}
X \bullet B = (X \oplus B) \ominus B
\end{split}$$</div>

膨張してから浸食する。
*X* の外部をOpeningすることと同義。
*X* 内部のヒビやチリなど、SEより小さい構造が塗りつぶされ、大きくなる。

## 応用

### Pattern Spectrum, サイズ分布

小さいSEから順に大きくしながら
Openingで削れた部分の面積を記録していく。
元画像の面積で割ったものはサイズ密度関数(size density function)と呼ばれる。
細かいギザギザを含む図形ほど小さいSEで削れる成分が多い。
要約統計量としてはエントロピーが使える(ギザギザ持ちほど高い)。

### Top Hat

<div>$$\begin{split}
X - (X \circ B)
\end{split}$$</div>

元画像とopeningの差。
トンガリの先っちょや背景のノイズ成分が得られる。

### Black Hat

<div>$$\begin{split}
X - (X \bullet B)
\end{split}$$</div>

元画像とclosingの差。
*X* 内のヒビやトゲの根元らへんが得られる。

### Morphological gradient

<div>$$\begin{split}
(X \oplus B) - (X \ominus B)
\end{split}$$</div>

dilationとerosionの差。
エッジ検出法のひとつ。
X上の境界が欲しい場合は$X - (X \ominus B)$。
背景側の境界が欲しい場合は$(X \oplus B) - X$。


## ライブラリ

[scipy.ndimage](http://docs.scipy.org/doc/scipy/reference/tutorial/ndimage.html)
: SciPyのサブモジュール。`np.array`を利用する。

[scikit-image](http://scikit-image.org/)
: Pythonモジュール。
  SciPyコミュニティの人たちが作ってて`scipy.ndimage`の拡張だと言っているので、
  とりあえずこっちを使っておけばいいんじゃなかろうか。

[OpenCV (Open Source Computer Vision)](http://opencv.org/)
: C++、Pythonなど。信頼と実績があるらしく、書籍やネット上の情報も多い。

[CImg](http://cimg.eu/)
: C++。ヘッダひとつincludeするだけ。
  ドキュメントも良さげ。

[imager](http://dahtah.github.io/imager/)
: R。新しめでドキュメントも充実。内部でCImgを利用。

[mmand](https://github.com/jonclayden/mmand)
: R。READMEは良さげ。

[Morpho](https://github.com/zarquon42b/Morpho)
: R。ドキュメント不足。

----

<a  href="http://www.amazon.co.jp/gp/product/4621082949/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=4621082949&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621082949&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=4621082949" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a  href="http://www.amazon.co.jp/gp/product/4862460844/ref=as_li_ss_il?ie=UTF8&camp=247&creative=7399&creativeASIN=4862460844&linkCode=as2&tag=heavywatal-22"><img border="0" src="http://ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4862460844&Format=_SL250_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="http://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=as2&o=9&a=4862460844" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
