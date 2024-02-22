+++
title = 'RepeatMasker'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<https://www.repeatmasker.org/>

## インストール

### プログラム + データベース

- <https://hub.docker.com/r/dfam/tetools>
- <https://github.com/Dfam-consortium/TETools/>

全部入りの [Dockerコンテナ]({{< relref "docker.md" >}}) を使うのが楽。
データベースも含めて20GB以上ダウンロードし、80GB以上ストレージを使う。

```sh
docker image pull dfam/tetools
docker container run -dit --mount type=bind,source="$PWD",target=/work --workdir /work --user "$(id -u):$(id -g)" --name dfamtet dfam/tetools
docker container exec dfamtet RepeatMasker | head
docker container exec dfamtet rmblastn -version
docker container exec dfamtet trf -v
```

```sh
apptainer pull dfam-tetools_1.sif docker://dfam/tetools:1
apptainer exec dfam-tetools_1.sif rmblastn -version
```

### プログラム

プログラムだけを入れるなら
[Homebrew]({{< relref "homebrew.md" >}}) でも可能だが、
numpyや古いPythonに依存していて不便:
```sh
brew info brewsci/bio/repeatmasker
```

依存プログラムも自動的に入る:
[RMBlast](https://www.repeatmasker.org/RMBlast.html),
[HMMER](http://hmmer.org/),
[Tandem Repeat Finder](https://tandem.bu.edu/trf/trf.html),
[h5py](https://docs.h5py.org/).


### データベース

単純な反復や一般的なアーティファクトはライブラリに組み込まれているが、
もっとしっかり使いたい場合はライブラリを更新して利用する。

1.  ディレクトリ移動:
    `cd $(brew --prefix)/opt/repeatmasker/libexec/`
1.  ファイルを追加・更新する。
    - [RepBase](https://www.girinst.org/server/RepBase/):
      アカデミックな用途なら無料で使えていたが、いつの間にか有料になっていた。
    - [Dfam](https://www.dfam.org/):
      使いたいものをダウンロードして `./Library/Dfam.*` を差し替える。<br>
      RepeatMasker 4.1.0 以前は `Dfam.hmm` や `Dfam.embl` を使っていたが、
      4.1.1 以降では `Dfam.h5` を使う。
      4.1.2 には始めから Dfam 3.3 curatedonly が付いてくる。
1.  変更を反映させる: `/usr/bin/perl ./configure <configure.input`

実行時に毎回 `-lib` オプションでファイルを指定する手もある。


## 使用方法

<https://www.repeatmasker.org/webrepeatmaskerhelp.html>

### コマンド

<https://github.com/rmhubley/RepeatMasker/blob/master/repeatmasker.help> :

```sh
RepeatMasker -help
RepeatMasker -pa 4 -qq -species oryza -dir . -xsmall -gff seq.fa
```

入力ファイルは圧縮 `.fa.gz` でもいいけど勝手に展開してしまうので注意。

`-engine [crossmatch|wublast|abblast|ncbi|hmmer|decypher]`

`-parallel 1`
:   並列化の恩恵は大きい

`-s` (slow), `-q` (quick), `-qq` (rush)
:   sensitivityとのトレードオフ

`-nolow` / `-low`
:   low complexity DNAをマスクしない

`-noint` / `-int`
:   interspersed repeatsをマスクしない

`-norna`
:   small RNAをマスクしない

`-div [number]`
:   コンセンサスからの分化度が指定したパーセント未満のやつだけマスク

`-lib [filename]`
:   自分で用意したライブラリを使う。
:   そのときのスコア閾値: `-cutoff 225`

`-species CLADE_NAME`
:   <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/>
    に出てくるクレード名で指定可能。
    大文字小文字は無視。
    chimpanzee, mammals みたいな英語名も一部可能。
    デフォルトは primate らしい。

`-frag 60000`
:   最大マスク長

`-nopost`
:   最後に自動で `PostProcess` を走らせない。

`-dir OUTDIR`
:   出力先。デフォルトはカレントではなくクエリと同じとこ。

`-xsmall`
: 反復配列を小文字にするsoft mask。
  デフォルトでは `N` に置き換えるhard mask。

`-gff`
:   GFFファイルも出力する。

`-small`
: "returns complete .masked sequence in lower case"
: 意味不明。ソースコードを眺めた感じでは、何もしてない。


### 生成物

計算途中の一時ファイルが `./RM_{pid}.{datetime}` に書き出される。

`${INFILE}.cat.gz`
:   RepeatMasker が出力する大元の結果。
    以下のファイルはこれを元に `ProcessRepeats` が作る。

`${INFILE}.masked`
:   見つかった箇所を `N` や小文字に置き換えた配列ファイル。
:   入力ファイルの折り返し幅を保持してくれない。
:   soft mask済み入力ファイルの小文字を保持してくれない。
    判定から外れた部分は大文字に戻される。
    追加マスクしたいならGFFとかを使って自分でやる必要がある。

`${INFILE}.tbl`
:   見つかった反復配列の要約

`${INFILE}.out`
:   アノテーション情報。
:   固定幅っぽい変なレイアウトの表で扱いにくい。ちょっと眺めるだけ。

`${INFILE}.out.gff`
: `-gff` オプションを付ければ作ってくれる。
  実質的に使える出力ファイルはこれだけかも。


## RepeatScout

<http://bix.ucsd.edu/repeatscout/>

反復配列を *de novo* で拾い、RepeatMaskerで利用可能なライブラリを生成する。

1.  RepeatScout本体をインストール:
    ```sh
    wget -O- http://bix.ucsd.edu/repeatscout/RepeatScout-1.0.5.tar.gz | tar xz
    cd RepeatScout-1/
    make
    ```
1.  L-mer の頻度テーブルをつくる:
    ```sh
    build_lmer_table -l 14 -sequence myseq.fa -freq lmer_table
    ```
1.  そのテーブルと配列から反復配列のFASTAを作る:
    ```sh
    RepeatScout -sequence myseq.fa -output rs_output.fa -freq lmer_table -l 14
    ```
1.  TRFとNSEGを呼び出して &gt;50% low-complexity なものを除外:
    ```sh
    cat rs_output.fa | filter-stage-1.prl >rs_filtered1.fa
    ```
    [NSEG](ftp://ftp.ncbi.nih.gov/pub/seg/nseg/) はビルド不可能なので
    `filter-stage-1.prl` を適当に書き換える必要がある。
1.  RepeatMaskerで位置と登場回数を調べる:
    ```sh
    RepeatMasker -parallel 4 -dir . -lib rs_filtered1.fa myseq.fa
    ```
1.  一定回数に満たないものを除外:
    ```sh
    cat rs_filtered1.fa | filter-stage-2.prl --thresh=10 --cat=myseq.fa.out >rs_filtered2.fa
    ```
1.  遺伝子領域のGFFなどを与え、mobile elementっぽくないものを除去:
    ```sh
    compare-out-to-gff.prl --gff=known_genes.gff --cat=myseq.fa.out --f=rs_filtered2.fa >lib.ref
    ```


## RepeatModeler

<https://www.repeatmasker.org/RepeatModeler/>

上記RepeatScout手順を簡単に実行するラッパー？
私はうまく使えた試しがない。

### 前準備

-   RepeatMasker, TRF, RepeatScout: 上記
-   CD-Hit: `brew install brewsci/bio/cd-hit`
-   RECON:
    ```sh
    wget -O- http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz | tar xz
    cd RECON-1.08/src/
    make
    make install
    ```
    `#include <ctypes.h>` を明示的に書けと怒られるので書く。
    `bin/` にPATHを通す。
-   RepeatModeler 本体:
    ```sh
    wget -O- https://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.3.tar.gz | tar xz
    cd RepeatModeler-2.0.3/
    ./configure
    ```
    Perlモジュールをインストールせよと言われたら `cpan JSON` とか適当に。


### 使い方

1.  カレントディレクトリにBLASTデータベースを構築:
    ```sh
    BuildDatabase -name Colletotrichum_orbiculare -engine ncbi path/to/Colletotrichum_orbiculare.fa
    ```
1.  本体を実行(かなり時間がかかる):
    ```sh
    RepeatModeler -engine ncbi -pa 4 -database Colletotrichum_orbiculare >run.out
    ```
    `RM_[PID].[DATE]/` に結果が書き出される。
1.  できあがった `consensi.fa.classified` をライブラリとして RepeatMasker を実行:
    ```sh
    RepeatMasker -lib consensi.fa.classified some_sequence.fa
    ```

## TEclass

<https://www.compgen.uni-muenster.de/tools/teclass>

トランスポゾンを分類する。

与えられたFASTAに含まれているのはTEだ、という仮定で分類するだけので、
単純反復配列や重複遺伝子などを予めしっかり除去しておく必要がある。

1.  本体をダウンロード:

        wget -O- http://www.compgen.uni-muenster.de/tools/teclass/download/TEclass-2.1.3.tar.gz | tar xz
        cd TEclass-2.1.3/
        less README

1.  周辺ライブラリを整備。
    でもとりあえず分類したいだけならblastclustなどは不要らしい:

        ./Download_dependencies.sh
        ./Compile_dependencies.sh
        ./Configure.pl

1.  pre-built classifiers (&gt;400MB) をダウンロード:

        cd path/to/TEclass/classifiers/
        wget -O- http://www.compgen.uni-muenster.de/tools/teclass/download/classifiers.tar.gz | tar xz

1.  実行: `./TEclassTest.pl file.fa`

1.  結果はひとつのディレクトリにまとめて書き出される
    -   `file.fa`: 元ファイルからTEだけ抜き出したもの？
    -   `file.fa.html`: 一覧
    -   `file.fa.lib`: RepeatMasker用？
    -   `file.fa.stat`: LTRなどがそれぞれいくつあったか集計
