+++
title = '擬似乱数生成器'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

{{%div class="warning"%}}
-   `<cstdlib>` の `std::rand()` は乱数の質も悪く、速度も遅いので非推奨。
-   `<algorithm>` の `std::random_shuffle()` は引数省略で
    `std::rand()` が使われてしまうので非推奨。
{{%/div%}}

-   外部の生成器としてはSFMTやdSFMTが高速で高品質。
-   C++11 からは新しい `<random>` が標準ライブラリに追加されてまともに使えるっぽい。
-   C++11 の `<algorithm>` の `std::shuffle()`
    は生成器を明示的に渡す必要があるのでよい。

## Mersenne Twister

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html>

松本眞と西村拓士によって開発された高速・高品質な擬似乱数生成器。
標準の `<random>` でも利用可能。

### SFMT

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/>

Mersenne Twisterを松本眞と斎藤睦夫がさらに改良したもの。
2倍速く、より均等な分布になってるらしい。整数乱数。
標準には含まれず、ソースからのビルドが必要。

### dSFMT

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/#dSFMT>

SFMTの `double` 版。
整数乱数から除算で変換するよりも良質で高速。
v2.1からは整数も出力可能。
標準には含まれず、ソースからのビルドが必要。

{{%div class="note"%}}
SFMTやdSFMTを簡単に導入するためのインストーラを作った。

<https://github.com/heavywatal/sfmt-installer/releases>
{{%/div%}}

## Xorshift

[<http://www.jstatsoft.org/v08/i14>](http://www.jstatsoft.org/v08/i14)

周期が2\^128でSFMTほど良質ではないらしいが、生成は超高速で、
何より実装が簡単

```c++
// シードは４要素の配列。どっかで一度適当に定義すること。
extern unsigned int seed128[4];

// シードを与える関数
inline void init_xorshift(unsigned int s){
    for (unsigned int i=0; i<4; ++i) seed128[i]=s=1812433253U*(s^(s>>30))+i;
}

// 32bitの整数乱数を生成
inline unsigned int xorshift128(){
    unsigned int *a(seed128);
    unsigned int  t(a[0]^(a[0]<<11));
    a[0] = a[1]; a[1] = a[2]; a[2] = a[3];
    return a[3]=(a[3]^(a[3]>>19))^(t^(t>>8));
}
```

## `<random>`

-   <http://www.cplusplus.com/reference/random/>
-   <http://en.cppreference.com/w/cpp/numeric/random>
-   <https://sites.google.com/site/cpprefjp/reference/random>

C++11 ではまともに使える乱数ライブラリが追加された。
乱数生成エンジンと分布生成器を組み合わせて使う。

## Seed

### `/dev/urandom`

`/dev/random` は擬似乱数を生成するデバイスで、
環境ノイズから乱数を生成するため自然乱数に近い。
`/dev/urandom` は十分なノイズが蓄積していなくても
内部プールの再利用によってすぐに生成してくれるが、
そのために `/dev/random` ほど安全ではない。
長期に渡って使われる暗号鍵の生成以外の目的では
`/dev/urandom` の利用が推奨されている。
LinuxやOSXで利用可能。

```c++
##include <fstream>

unsigned int dev_urandom() {
    unsigned int x;
    try {
        std::ifstream fin("/dev/urandom", std::ios::binary | std::ios::in);
        fin.exceptions(std::ios::failbit);
        fin.read(reinterpret_cast<char*>(&x), sizeof(unsigned int));
    }
    catch (std::ios::failure& e) {throw std::ios::failure("/dev/urandom");}
    return x;
}
```

### `std::random_device`

C++11ではそのための関数が `<random>` に用意されている。
これも基本的には `/dev/urandom` から生成するらしい

```c++
std::random_device rd;
const std::random_device::result_type seed = rd();
```
