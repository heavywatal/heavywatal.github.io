+++
title = '複雑ネットワーク'
tags = ["math"]
[menu.main]
  parent = "bio"
+++

<a href="http://www.amazon.co.jp/exec/obidos/ASIN/4764903636/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/412zEm9DGuL._SX160_.jpg" alt="複雑ネットワーク―基礎から応用まで" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/1107626250/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51f41OuhBOL._SX160_.jpg" alt="Dynamical Processes on Complex Networks" /></a>
<a href="http://www.amazon.co.jp/exec/obidos/ASIN/B005PS507U/heavywatal-22/" rel="nofollow" target="_blank"><img src="http://ecx.images-amazon.com/images/I/51F6oA7jkwL._SX150_.jpg" alt="Analysis of Biological Networks (Wiley Series in Bioinformatics)" /></a>
## 用語

### Node / Vertex

ノード、頂点。
GRNでは遺伝子。

### Edge / Link

エッジ、枝。ノードとノードを結ぶ線。
方向性を持つ場合と持たない場合がある。
GRNでは転写促進・抑制の制御関係なので有向。

## ネットワークの特徴量

### Order

ネットワーク中のノードの数

### Size

ネットワーク中のエッジの数

### Degree

次数。ひとつのノードがもつエッジの数

### Average Degree

平均次数。
次数をネットワーク中のノードで平均したもの。

### Density / Connectance

密度あるいは結合度。
ネットワーク中に存在しているエッジの数を、最大可能エッジ数で割ったもの。
最大可能エッジ数は、ノード数と自己制御の有無によって決まる。

自己制御あり有向グラフの最大可能エッジ数: $V ^ 2$\
自己制御なし有向グラフの最大可能エッジ数: $V (V - 1)$\
自己制御なし無向グラフの最大可能エッジ数: $V (V - 1) / 2$

### Clustering Coefficient

クラスター係数。
あるノードから見て、隣接する2つのノード同士もエッジで繋がっていると三角形ができる。
この三角形が多いほどクラスター係数が大きくなる。
ネットワーク中のノードについて平均したのが平均クラスタ係数。

### The number of selfloops

自己制御数。
有向グラフにおいて、両端が同じノードに接続しているエッジの数。

### Degree assortativity

次数相関。
隣接する2つのノードの次数が似ているほど高くなる。
次数の高いハブ的なノードが同じようにハブ的なノードと接続しがちな場合、assortative。
逆に、ハブに対して次数の低いノードが接続しがちな場合、disassortative。
