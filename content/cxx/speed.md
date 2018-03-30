+++
title = 'C++高速化'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

## はじめに

-   速いプログラムで得られるメリットを超えるようなコストを払わないように。
    まずは動くプログラムを書いて目的を達成することが大事。
    自分律速じゃなくてプログラム律速だなと感じた段階でリファクタリングを考える。
-   プログラム本来の意図が読み取れなくなりそうなマニアックな高速化は避ける。
    清く正しくメンテナンスしやすいプログラムを書くほうが結局は生産的。
-   学習目的でない限り、車輪の再発明を避ける。
    やろうとしていることはきっと既に誰かが実現し、
    再利用可能な形で公開してくれているはず。
    まずは標準ライブラリとかBoostを探してみる。
    あとGitHubでスターが多いやつとか。

## 頑張れコンパイラ

Intelの icc でビルドされたプログラムは速いらしい。
gcc や clang の最適化技術も着々と進歩しており、
新しいコンパイラを使うほうがその恩恵を受けられる。

### 最適化オプション

<https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>

```sh
% g++ -O3 main.cpp
```

### コンパイル時定数 `constexpr`

コンパイル時に計算できるものを予め計算して定数にしておくことで、
実行時の計算を減らすことができる。
C++11では明示的に `constexpr` 修飾することができるので積極的に使うべし。
テンプレートメタプログラミングという手もある。
C++11を使えない場合でも、
定数同士の四則演算などはまとめておけばコンパイラが計算してくれるはず。

```c++
// 掛け算も割り算も毎回計算
return 4.0 * M_PI * std::pow(radius, 3) / 3.0;

// コンパイル時定数を掛けるだけ
constexpr double c = 4.0 * M_PI / 3.0;
return c * std::pow(radius, 3);
```

また、四則演算の中でも割り算は特に遅いらしいので、
逆数のかけ算になるような書き方のほうが高速になるかもしれない。


### インライン展開

関数の呼び出しにはコストがかかり、
コードが短くて何度も呼ばれる関数では特にそのコストの割合が高くて馬鹿にできない。
呼び出される予定のところにコンパイルの段階で関数を展開してやること
(=インライン展開) でそのコストを無くせる。

ファイルをまたいで使う関数の場合は定義宣言を
`ソース.cpp` ではなく `ヘッダ.h` に置き、頭に `inline` を付ける。
メンバ関数の場合はクラス定義の中で定義宣言にするだけでよい。
ただしそのように書いてもコンパイラへのヒントになるだけで、
インライン展開に必要でも十分でもない。

`-O2` では `-finline-small-functions` がオンになり、
`-O3` では `-finline-functions` がオンになる。
実行ファイルのサイズがデカくなるし、コンパイルは遅くなるので、バランスを考える。

```c++
inline int pow2(int x) {
    return x *= x;
}
```

### 関数オブジェクト

関数オブジェクトとは、`operator()` を定義した `class` または `struct` のこと。
`std::for_each()` や `std::transform()` の最後の引数に関数を渡す場合など、
何度も呼び出される小さい関数については、
普通の関数や関数ポインタで渡すよりも関数オブジェクトを使ったほうがよい。
コンパイラによる最適化がしやすいらしい。
メンバ変数として値を保持できるので、
引数のやり取りやメモリの確保・開放が少なくてすむという場面もありそう。

```c++
class isAboveThreshold {
  public:
    isAboveThreshold(const double t): threshold_(t);
    bool operator()(const double x) const {
        return x > threshold_;
    }
  private:
    const double threshold_;
};

std::transform(v.begin(), v.end(), result.begin(), isAboveThreshold(3.14));
```

C++11からはラムダ式が便利

```c++
const double threshold = 3.14;
std::transform(v.begin(), v.end(), result.begin(),
               [threshold](double x) -> bool {return x > threshold;});
```

## 余計な一時オブジェクトを作らない

普通にプログラムを書くと、思わぬところで余計な一時オブジェクトが作られることになる。
メモリもcpu時間ももったいないので、なるべく避けよう。

### `const`参照渡し or ポインタ渡し

オブジェクトの値渡しにはコピーのコストが生じるので、
STLコンテナや自作クラスなどは参照渡しのほうが高速。
ただし参照を解決するコストも無いわけではないので、
`int`や`double`のような単純な型はむしろ普通の値渡しがよい。

参照渡しされた仮引数の中身を関数内で変更すると、
呼び出し元の実引数の値も変更されることになる。
変更するつもりがない場合は明示的に `const` をつけて受け取るようにすると安全

```c++
void do_something(const std::vector<int>& huge_array) {
    // コピーせず本体への参照のみ受け取る。
    // constがついているので、変更しようとするとコンパイルエラー
}
```

引数の中身を変更するつもりがある場合はポインタを受け取るようにすると、
呼び出すときに `&` が必要となるので意図が伝わりやすい。

```c++
void do_something(std::vector<int>* huge_array) {
    // コピーせずポインタを受け取る。
    // huge_arrayの変更は呼び出し元にも影響する。
}

do_something(&some_array);  // ああ、この配列は変更されるんだな
```


### ムーブセマンティクス (C++11)

受け取ったものを関数の中で変更して返すような関数を作るとき、
引数の渡し方も受け取り方も3パターンずつ考えられる。

- 渡し方
    - lvalue (`vec0`)
    - xvalue (`std::move(vec0)`)
    - prvalue (`Vector{}`)
- 受け取り方
    - const lvalue参照 (`const Vector&`)
    - lvalue値 (`Vector`)
    - const lvalue参照 (`const Vector&`) とrvalue参照 (`Vector&&`) のオーバーロード

[それぞれの組み合わせでコピーとムーブが何回起こるかテストするコード]
(https://github.com/heavywatal/sketchbook/blob/master/c%2B%2B/move.cpp)

const lvalue参照で受け取ったものを変更するためにはコピーが必要になるので、
普通のlvalue渡しでは問題ないけどrvalue渡しには適さない。
rvalue参照受け取りでオーバーロードすると、
コピーせず1回のムーブで済ませられる。
ただし引数が増えると手に負えなくなる。
lvalue値受け取りの関数はrvalueを受け取るときにはコピーを生じない。
オーバーロードの場合と比べてムーブが余計に1回生じるが、
そのコストと定義の手軽さを天秤に掛けて考慮する価値はある。

`return`最適化にはムーブすらしない"copy elision"とムーブの2種類があって、
ローカル変数やprvalueを返す場合は両方を試み、
引数を返す場合は後者のみを試みるらしい。

関数から値を返すときにムーブ返ししたくなるところだが、
多くの場合コンパイラがうまいことやってくれるのでわざわざ
`return std::move(output);` などと書く必要はない。
やってしまうと、むしろコンパイラによる最適化を妨げるかもよと警告される。
ただし、返る変数と関数定義で型が一致せず暗黙の型変換が挟まる場合は、
明示的にムーブ返しする必要がある。

関連記事はたくさん見つかるが、特に読みやすく参考になったのはこちら:

> [本当は怖くないムーブセマンティクス - yohhoyの日記（別館）](http://yohhoy.hatenablog.jp/entry/2012/12/15/120839) \
> [参照渡し or 値渡し？ - yohhoyの日記](http://d.hatena.ne.jp/yohhoy/20120524/p1)


### `noexcept`

クラスに自前のデストラクタやコピーコンストラクタを書くと、
暗黙のムーブコンストラクタが作られなくなり、
ただ移動したいときにもコピーが行われてしまう。
自前のムーブコンストラクタを書いても、
`noexcept` が添えられていないと `std::move_if_noexcept()` による移動
(例えば `std::vector` のリアロケーション) ではコピーされてしまう。

[ムーブコンストラクタの定義の仕方でコピー・ムーブのされ方が変わることを確かめるコード]
(https://github.com/heavywatal/sketchbook/blob/master/c%2B%2B/move_if_noexcept.cpp)

コンストラクタ類を必要に応じて書き足していくと忘れそうなので、
`Cls (Cls&&) noexcept = default;`
のような明示的default/deleteを最初に全部書いてしまうほうがいいかも。
さらに `static_assert()` でコンパイル時に
`std::is_nothrow_move_constructible` などをチェックしておくと安心。

それ以外の部分でも「この関数は例外を投げない」
と宣言したほうがコンパイラにとって最適化しやすくなる。


## コンテナ

### 用途に合わせる

http://en.cppreference.com/w/cpp/container

`std::array`
:   Cの配列に便利なメンバ関数を持たせたようなクラス。
    長さはコンパイル時に固定で、
    メモリ上に連続した領域を確保する。
    要素へのアクセスは高速で、
    インデックス`[]` によるランダムアクセスも可能。

`std::vector`
:   可変長にした `array` 、つまり実行時に伸長・縮小可能。
    C++で配列作るならとりあえずこれ。
    途中に `insert()` するのは遅い
    （メモリ領域を別のところに再確保して全部コピーしなければならないので）。

`std::valarray`
:   要素ごとの四則演算など関数・演算子が定義済みの `vector` 亜種。
    便利なだけでなくきっと最適化もされやすい。
    ただし長さの変更など苦手な点もあるので使い所は限られる。
    本格的なベクタ演算・行列演算がしたければ
    [Eigen](http://eigen.tuxfamily.org/) や
    [Armadillo](http://arma.sourceforge.net/) などを使ったほうよさそう。

`std::deque`
:   `vector` とほぼ同じだが、`reserve()` ができない。
    代わりに、`push_front()` と `pop_front()` が使える。
    つまり、先頭に対する追加・削除が必要な場合だけ使うコンテナ。

`std::list`
:   飛び飛びのメモリ領域に散らばって存在できる配列クラス。
    これは、各要素の値とともに前後の要素を示すイテレータを格納することで実現されている。
    そのため、途中への `insert()` はイテレータを書き換えるだけなので高速。
    ただし、その分メモリは余分に食うし、ランダムアクセスできないし、イテレータで総なめするのも遅い。

`std::unordered_set`, `std::unordered_map`
:   順序を気にしないぶん `std::set` や `std::map` より高速。
    ただしハッシュ関数の準備など多少めんどい。


### メモリは一気に確保

`std::vector` の `push_back()` は勝手にメモリ領域を確保してくれるので、
大きさの心配をする必要がなくて便利。
これは多くの場合、領域が足りなくなる都度「別のところに倍の領域を確保してコピー」
という処理をすることによって行われる。
始め1、次2、4、8、16、、、という具合に。
なので、10000個の要素を格納すると分かっている場合には、始めからそれだけ確保しておくと速い。

```c++
std::vector<int> v;
v.reserve(10000);
// then push_back() many times
```


## ループ

### 継続条件

`for` や `while` を回すときの継続条件式は回る度に評価されるので意外とコストになるかも。
イテレータで回すときの `v.end()` も毎回呼び出されるので、
ループの中身が軽い処理の場合には無視できない差になるかもしれない。

```c++
for (size_t i=0; i < v.size(); ++i) {
    // v.size() many times!
}

for (size_t i=0, n=v.size(); i<n; ++i) {
    // v.size() only once
}

for (vector<int>::iterator it=v.begin(), v_end=v.end(); it!=v_end; ++it) {
    // v.end() only once
}

for (const auto& x: v) {
    // C++11 range-based for
}
```

### 前置インクリメント

後置インクリメント `i++` でも
前置インクリメント `++i` でも `for` ループの結果は変わらない。
が、`i++` だと前の値を記憶してからプラスする（一時オブジェクトが作られる）ので、
`++i` のほうがいいらしい。特にイテレータのとき。

```c++
for (std::vector<int>::iterator it=v.begin(); it!=v.end(); ++it) {
    // do something
}
```

### 出来る限り外で処理

`if-else` や `try-catch` のブロックをループの外で大きく取れないか、
一時変数の定義や演算などをループの外で予めやっておけないか、確認すべし。

## 入出力

### 標準入出力

Cストリーム(`std::printf` とか)とC++ストリーム(`std::cout` とか)
が混在するプログラムでもちゃんと関数が呼ばれた順に入出力を行うため、
デフォルトではこれら２つのストリームが同期するようになっている。
必要がなければ切っておく。

`std::cin` はデフォルトで `std::cout` に結びつけられてて、
`std::cin` される度に `std::flush` されてしまうらしいので、
そうならないように切り離す。

```c++
#include <iostream>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
}
```

## 関連書籍

<a href="https://www.amazon.co.jp/Effective-Modern-_C-11-14%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%A0%E3%82%92%E9%80%B2%E5%8C%96%E3%81%95%E3%81%9B%E3%82%8B42%E9%A0%85%E7%9B%AE/dp/4873117364/ref=as_li_ss_il?_encoding=UTF8&psc=1&refRID=2W1SJD1QE3Z43HGNHQPV&linkCode=li2&tag=heavywatal-22&linkId=b5ae72d9c28f98283a01af1b110e1c8c" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4873117364&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li2&o=9&a=4873117364" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Effective-3-ADDISON-WESLEY-PROFESSIONAL-COMPUTI/dp/4621066099/ref=as_li_ss_il?ie=UTF8&linkCode=li2&tag=heavywatal-22&linkId=920cd9757c2ac53cfd19d947c4d1d72b" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4621066099&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li2&o=9&a=4621066099" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/gp/product/4894714108/ref=as_li_ss_il?ie=UTF8&linkCode=li2&tag=heavywatal-22&linkId=dc3cda30d8eee4b5b5193d3835261ad6" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4894714108&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li2&o=9&a=4894714108" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
<a href="https://www.amazon.co.jp/Efficient-C-%E3%83%91%E3%83%95%E3%82%A9%E3%83%BC%E3%83%9E%E3%83%B3%E3%82%B9%E3%83%97%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%9F%E3%83%B3%E3%82%B0%E3%83%86%E3%82%AF%E3%83%8B%E3%83%83%E3%82%AF-%E3%83%96%E3%83%AB%E3%82%AB-%E3%83%80%E3%83%96/dp/4894712458/ref=as_li_ss_il?s=books&ie=UTF8&qid=1477818040&sr=1-1&keywords=efficient+c++&linkCode=li2&tag=heavywatal-22&linkId=728f05422ad565bce17f6ea3d1159c7d" target="_blank"><img border="0" src="//ws-fe.amazon-adsystem.com/widgets/q?_encoding=UTF8&ASIN=4894712458&Format=_SL160_&ID=AsinImage&MarketPlace=JP&ServiceVersion=20070822&WS=1&tag=heavywatal-22" ></a><img src="https://ir-jp.amazon-adsystem.com/e/ir?t=heavywatal-22&l=li2&o=9&a=4894712458" width="1" height="1" border="0" alt="" style="border:none !important; margin:0px !important;" />
