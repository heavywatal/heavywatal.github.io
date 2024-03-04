+++
title = '擬似乱数生成器'
tags = ["c++"]
[menu.main]
  parent = "cxx"
+++

-   `<cstdlib>` の `std::rand()` は乱数の質も悪く、速度も遅いので非推奨。
    C++11 から標準ライブラリに追加された `<random>` を使う。
-   `<algorithm>` の `std::random_shuffle()` は引数省略で
    `std::rand()` が使われてしまうので非推奨。
    C++11 で追加された `std::shuffle()` に生成器を明示的に渡して使う。
-   非標準の生成器としてはSFMTやdSFMTが高速で高品質。Xorshift系とPCGの動向にも注目。


## `<random>`

-   <https://cplusplus.com/reference/random/>
-   <https://en.cppreference.com/w/cpp/numeric/random>
-   <https://cpprefjp.github.io/reference/random.html>

C++11 ではまともに使える乱数ライブラリが追加された。
乱数生成エンジンと分布関数オブジェクトを組み合わせて使う。

```c++
#include <iostream>
#include <random>

int main() {
    // seed
    std::random_device rd;
    const auto seed = rd();

    // engine
    std::mt19937 rng(seed);

    // probability density distribution
    const double mean = 0.0;
    const double sd = 1.0;
    std::normal_distribution<double> dist(mean, sd);

    // generate!
    for (size_t i=0; i<8; ++i) {
        std::cout << dist(rng) << std::endl;
    }

    return 0;
}
```

## Mersenne Twister

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/mt.html>

松本眞と西村拓士によって開発された高速・高品質な擬似乱数生成器。
標準の `<random>` でも利用可能になっており、
パラメータ定義済みの `std::mt19937` がよく使われる。

### SFMT

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/>

Mersenne Twisterを松本眞と斎藤睦夫がさらに改良したもの。
SIMD命令を利用して、速度も品質も向上したらしい。
標準には含まれず、ソースからのビルドが必要。

### dSFMT

<http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/#dSFMT>

SFMTの `double` 版。
整数乱数から除算で変換するよりも良質で高速。
v2.1からは整数も出力可能。
標準には含まれず、ソースからのビルドが必要。

### インストール方法、使い方

<https://github.com/heavywatal/sfmt-class/>

SFMTやdSFMTを簡単に導入するためのインストーラを作って公開した。
C++標準 `<random>` の `std::mt19937`
と同じように使えるようにしたラッパークラス
(`wtl::sfmt19937`) も書いた。

## Xorshift

[Marsaglia (2003)](https://www.jstatsoft.org/article/view/v008i14)
から始まって、いくつかの改良版が派生している。
省メモリで高速。
周期は $2^{128}$ ほどでMTに比べれば短いけど、大概はこれで十分。
最新情報は
[xoroshiro.di.unimi.it](http://xoroshiro.di.unimi.it/)
で追えば良さそう。

## PCG

線形合同法(LCG)の出力をpermutationした[PCG](http://www.pcg-random.org/)も良さそう。
作者O'Neillのブログは読み応えがあっておもしろい。


## Seed

### `/dev/urandom`

`/dev/random` は擬似乱数を生成するデバイスで、
環境ノイズから乱数を生成するため自然乱数に近い。
`/dev/urandom` は十分なノイズが蓄積していなくても
内部プールの再利用によってすぐに生成してくれるが、
そのために `/dev/random` ほど安全ではない。
長期に渡って使われる暗号鍵の生成以外の目的では
`/dev/urandom` の利用が推奨されている。

```c++
#include <fstream>

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
これも基本的には `/dev/urandom` から生成するらしい。

```c++
std::random_device rd;
const std::random_device::result_type seed = rd();
```

大量のシミュレーションを回すときなど、
エントロピーがもっと必要な場合は `std::seed_seq` を利用する。
