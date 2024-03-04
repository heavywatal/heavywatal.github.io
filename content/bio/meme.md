+++
title = 'MEME'
subtitle = "モチーフ発見ツール"
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

- https://meme-suite.org/
- https://meme-suite.org/doc/overview.html
- [Bailey et al. 2009](https://academic.oup.com/nar/article/37/suppl_2/W202/1135092)
- [Bailey and Elkan 1994](https://www.ncbi.nlm.nih.gov/pubmed/7584402)

## インストール

https://meme-suite.org/doc/install.html

1.  ソースコードを[ダウンロード](https://meme-suite.org/doc/download.html)して展開:

        wget -O- https://meme-suite.org/meme-software/4.11.2/meme_4.11.2_1.tar.gz | tar xz

1.  `configure` してビルド:

        cd meme_4.11.2/
        ./configure --prefix=${HOME}/bio/local/meme --with-url=https://meme-suite.org --enable-build-libxml2 --enable-build-libxslt --enable-opt CC=clang
        make

1.  パスを通す:

        export PATH=${PATH}:/usr/local/meme/bin

## MEME

https://meme-suite.org/doc/meme.html

### 使い方

複数の配列が含まれるFASTAファイルを渡すだけ:

    meme sequences.fasta [options]

`-h`, `-version`
:   ヘルプ、バージョン表示

`-dna`, `-protein`
:   配列がDNAかタンパク質か (`-protein`)

`-maxsize`
:   入力ファイルの許容サイズ (`100000`)

`-nmotifs`, `-evt`
:   探索するモチーフ数を制御するため、
    個数そのものか E-value の上限を指定する。
    `-evt` を使うときは `-nmotifs` 大きめにしておく。
    (`-nmotifs 1`)

`-mod`
:   モチーフが配列上にどう分布しているか\
    `oops`: One Occurrence Per Sequence\
    `zoops`: Zero or OOPS\
    `anr`: Any Number of Repetitions

`-nsites`, `-minsites`, `-maxsites`
:   それぞれのモチーフがいくつ登場すると仮定するか
    (デフォルト値は `-mod` により異なる)

`-w`, `-minw`, `-max`
:   探索するモチーフの長さを指定
    (`-minw 8 -maxw 50`)

`-revcomp`
:   逆向きも考慮する

`-pal`
:   パリンドロームを探す

`-bfile <bfile>`
:   バックグラウンド配列を生成するマルコフ過程のパラメータを記述したファイルを指定。
    これを指定しない場合はトレーニング配列の塩基頻度のみを利用した0階マルコフ。
    FASTA配列からファイルを作ってくれるプログラム `fasta-get-markov` も用意されてる。
    <http://meme.nbcr.net/meme/doc/fasta-get-markov.html>:

        # order 0
        A 3.081e-01
        C 1.919e-01
        G 1.919e-01
        T 3.081e-01
        # order 1
        AA 1.078e-01
        AC 5.256e-02
        AG 5.908e-02
        AT 8.848e-02
        CA 6.519e-02
        CC 3.858e-02
        CG 2.908e-02
        CT 5.908e-02
        GA 6.239e-02
        GC 3.841e-02
        GG 3.858e-02
        GT 5.256e-02
        TA 7.284e-02
        TC 6.239e-02
        TG 6.519e-02
        TT 1.078e-01

一度適当に走らせてみて、出力結果
`meme.txt` の **COMMAND LINE SUMMARY** や
`meme.html` の **model parameters**
を見るとよい。デフォルト値もそこで分かる。


### スコア

[Bailey and Gribskov 1998](http://bioinformatics.oxfordjournals.org/content/14/1/48)

E-value
:   そのモチーフが同じサイズのランダムな配列の中にたまたま見つかる個数の期待値

Position p-value

Combined p-value

### モチーフの出力形式

LOGO
:   アルファベットの大きさで視覚的に表示

PSPM: position-specific probability matrix
:   ポジションごとの塩基・アミノ酸の相対的な頻度を実数[0, 1]の行列で表示。
    position weight matrix (PWM) と呼ぶことが多いような。

PSSM: position-specific scoring matrix
:   このあと MAST で使える形式の行列

BLOCKS, FASTA
:   そのモチーフを含む配列のID、開始位置、ヒットした領域の配列

Raw
:   モチーフにヒットした領域を切り出して並べただけ

regular expression
:   `[AT]` のように正規表現の文字集合を使った配列

## `DREME`

Discriminative Regular Expression Motif Elicitation

短いモチーフが得意で効率的。
background (negative) 配列を指定できる。
ChIP-seqデータではピーク周辺100bpくらいを使うべし。

## `MEME-ChIP`

長いモチーフが得意な `MEME` と
短いモチーフが得意な `DREME` を組み合わせて ensemble。

## `MAST`

<http://meme.nbcr.net/meme/doc/mast.html>

既知のモチーフ (`MEME` で発見されたとか) を配列データベースから検索する。

## References

- [バイオインフォマティクスのためのアルゴリズム入門](https://www.amazon.co.jp/dp/4320056507?&linkCode=ll1&tag=heavywatal-22&linkId=76478843e121eb67b247097dc03afec8)
- [バイオインフォマティクス](https://www.amazon.co.jp/dp/4621062514?&linkCode=ll1&tag=heavywatal-22&linkId=36ad1d7b5e767b5d4e8f068df2de77b4)
- [バイオインフォマティクス事典](https://www.amazon.co.jp/dp/4320056280?&linkCode=ll1&tag=heavywatal-22&linkId=1c78a849ed5ddcccea4826a023d95c10)
- [Pattern Recognition and Machine Learning](https://www.amazon.co.jp/dp/0387310738?&linkCode=ll1&tag=heavywatal-22&linkId=0dacd1cec1bcc3d73dc0a9f27d158183)
