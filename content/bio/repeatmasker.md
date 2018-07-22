+++
title = 'RepeatMasker'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

<http://www.repeatmasker.org/RMDownload.html>

## 前準備

[Homebrew]({{< relref "homebrew.md" >}}) を使う。
`brew tap brewsci/science brewsci/bio` しておく。

### プログラム

-   Perl 5.8 以上:

        % perl --version

-   検索エンジンをどれか1つ以上
    1.  [CrossMatch](http://www.phrap.org/)
    2.  [RMBlast](http://www.repeatmasker.org/RMBlast.html)
        (NCBI BLASTを RepeatMasker で使えるようにするためのラッパー)
        -   まず普通の [NCBI BLAST+]({{< relref "blast.md" >}}) をインストール。
            Homebrewでもできるけどビルドに時間かかるのでバイナリで。
        -   RMBlastをインストールし、元のBLASTの `bin` 内にシムリンクを張る:

                % brew install rmblast --without-blast
                % brew list rmblast
                % sudo ln -s $(brew --prefix)/bin/rmblastn /usr/local/ncbi/blast/bin/

    3.  [AB-BLAST](http://www.advbiocomp.com/blast.html) (有料)
    4.  [HMMER](http://hmmer.janelia.org/) & DFAM (human only):

            % brew install hmmer

-   [Tandem Repeat Finder](http://tandem.bu.edu/trf/trf.html) :

        % brew install trf

-   RepeatMasker 本体
    (どこでもいいけど今回は `/usr/local/` に入れる):

        % wget -O- http://www.repeatmasker.org/RepeatMasker-open-4-0-5.tar.gz | tar xz
        % sudo mv RepeatMasker /usr/local/

    Homebrewで入れる場合は `--without-configure` をつけて、
    後から自分で `./configure` する。

-   データベースを入れて `configure` (下記)

### データベース

RepeatMasker 用に加工された
[RepBase](http://www.girinst.org/repbase) ライブラリを用いる。
ただし取得には無料登録が必要で、承認まで時間がかかる。
また `@gmail.com` のような商用ドメインでは弾かれるので、
`@soken.ac.jp` のような所属機関アドレスが必要。

1.  `/usr/local/RepeatMasker/Libraries/` 直下に格納
    (念のためフォルダごと取っといてシムリンク。古いやつは適当に保持):

        % cd /usr/local/RepeatMasker/Libraries/
        % wget -O- --http-user="USERNAME" --http-password="PASSWORD" http://www.girinst.org/server/archive/RepBase20.04/protected/repeatmaskerlibraries/repeatmaskerlibraries-20140131.tar.gz | tar xz
        % mv Libraries repeatmaskerlibraries-20140131
        % mv RepeatMaskerLib.embl RepeatMaskerLib.embl.orig
        % ln -s repeatmaskerlibraries-20140131/RepeatMaskerLib.embl

2.  [Dfam](http://www.dfam.org/) のアップデートがあれば
    `Libraries/Dfam.hmm` ファイルを入れ替え:

        % mv Dfam.hmm Dfam.hmm.orig
        % wget ftp://selab.janelia.org/pub/dfam/Release/Dfam_1.4/Dfam.hmm.gz
        % gunzip Dfam.hmm.gz

3.  `configure` スクリプト実行:

        % cd ..
        % ./configure <~/bio/RepeatMasker.conf

    標準入力は毎回同じなのでファイルにしてしまったほうが楽チン

FASTA形式の反復配列を分類群単位で入手することもできる:

    % wget -o- --http-user="USERNAME" --http-password="PASSWORD" http://www.girinst.org/server/RepBase/protected/RepBase20.04.fasta/fngrep.ref

## 使用方法

<http://www.repeatmasker.org/webrepeatmaskerhelp.html>

### コマンド

<http://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html> :

    % /usr/local/RepeatMasker/RepeatMasker -h | less
    % less /usr/local/RepeatMasker/repeatmasker.help

`-engine [crossmatch|wublast|abblast|ncbi|hmmer|decypher]`

`-parallel 1`
:   並列化の恩恵は大きい

`-s` (slow), `-q` (quick), `-qq` (rush)
:   sensitivityとのトレードオフ

`-nolow`
:   low complexity DNAをマスクしない

`-noint`
:   interspersed repeatsをマスクしない

`norna`
:   small RNAをマスクしない

`-div [number]`
:   コンセンサスからの分化度が指定したパーセント未満のやつだけマスク

`-cutoff 225`

`-species CLADE_NAME`
:   <http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/>
    に出てくるクレード名で指定可能。
    大文字小文字は無視。
    chimpanzee, mammals みたいな英語名も一部可能。
    デフォルトは primate らしい。

`-frag 60000`
:   最大マスク長

`-nopost`
:   最後に自動で PostProcess を走らせない。

`-dir OUTDIR`
:   出力先。デフォルトはカレントではなくクエリと同じとこ。

`-gff`
:   GFFファイルも出力する。

### 生成物

`${INFILE}.masked`
:   見つかった箇所 `NNNNN` に置き換えたファイル

`${INFILE}.cat.gz`
:   RepeatMasker が出力する大元の結果。
    以下のファイルはこれを元に ProcessRepeats が作る。

`${INFILE}.tbl`
:   見つかった反復配列の要約

`${INFILE}.out`
:   アノテーション情報。
    TSVではなく固定幅に近い表で、例外も多くて使いにくい。

    {{%div class="warning"%}}
バグってるのでこのファイルを使ってはいけない！
    {{%/div%}}

## RepeatScout

<http://bix.ucsd.edu/repeatscout/>

反復配列を *de novo* で拾い、RepeatMaskerで利用可能なライブラリを生成する。

1.  RepeatScout本体をインストール:

        % brew install repeatscout --with-trf

2.  L-mer の頻度テーブルをつくる:

        % build_lmer_table -l 14 -sequence myseq.fa -freq lmer_table

3.  そのテーブルと配列から反復配列のFASTAを作る:

        % RepeatScout -sequence myseq.fa -output rs_output.fa -freq lmer_table -l 14

4.  TRFとNSEGを呼び出して &gt;50% low-complexity なものを除外:

        % cat rs_output.fa | filter-stage-1.prl >rs_filtered1.fa

    {{%div class="note"%}}
[NSEG](ftp://ftp.ncbi.nih.gov/pub/seg/nseg/) はビルド不可能なので
`filter-stage-1.prl` を適当に書き換える必要がある。
    {{%/div%}}

5.  RepeatMaskerで位置と登場回数を調べる:

        % RepeatMasker -parallel 4 -dir . -lib rs_filtered1.fa myseq.fa

6.  一定回数に満たないものを除外:

        % cat rs_filtered1.fa | filter-stage-2.prl --thresh=10 --cat=myseq.fa.out >rs_filtered2.fa

7.  遺伝子領域のGFFなどを与え、mobile elementっぽくないものを除去:

        % compare-out-to-gff.prl --gff=known_genes.gff --cat=myseq.fa.out --f=rs_filtered2.fa >lib.ref

## RepeatModeler

<http://www.repeatmasker.org/RepeatModeler.html>

上記RepeatScout手順を簡単に実行するラッパー？

{{%div class="warning"%}}
うまくできた試しがないので、RepeatScoutを直に動かしたほうが良い。
{{%/div%}}

### 前準備

-   RepeatMasker: 上記
-   RepeatScout: 上記
-   RECON:

        % brew install recon

-   RepeatModeler 本体:

        % wget -O- http://www.repeatmasker.org/RepeatModeler-open-1-0-8.tar.gz | tar xz
        % sudo mv RepeatModeler /usr/local/

    例によって `configure` も実行:

        % cd /usr/local/RepeatModeler
        % perl configure <~/bio/RepeatModeler.conf

### 使い方

1.  カレントディレクトリにBLASTデータベースを構築:

        % BuilDatabase -name Colletotrichum_orbiculare -engine ncbi path/to/Colletotrichum_orbiculare.fa

2.  本体を実行(かなり時間がかかる):

        % RepeatModeler -engine ncbi -pa 4 -database Colletotrichum_orbiculare >run.out

    `RM_[PID].[DATE]/` に結果が書き出される。

3.  できあがった `consensi.fa.classified` をライブラリとして RepeatMasker を実行:

        % RepeatMasker -lib consensi.fa.classified some_sequence.fa

## TEclass

<http://www.compgen.uni-muenster.de/tools/teclass/>

トランスポゾンを分類する。

与えられたFASTAに含まれているのはTEだ、という仮定で分類するだけので、
単純反復配列や重複遺伝子などを予めしっかり除去しておく必要がある。

1.  本体をダウンロード:

        % wget -O- http://www.compgen.uni-muenster.de/tools/teclass/download/TEclass-2.1.3.tar.gz | tar xz
        % cd TEclass-2.1.3/
        % less README

2.  周辺ライブラリを整備。
    でもとりあえず分類したいだけならblastclustなどは不要らしい:

        % ./Download_dependencies.sh
        % ./Compile_dependencies.sh
        % ./Configure.pl

3.  pre-built classifiers (&gt;400MB) をダウンロード:

        % cd path/to/TEclass/classifiers/
        % wget -O- http://www.compgen.uni-muenster.de/tools/teclass/download/classifiers.tar.gz | tar xz

4.  実行:

        % ./TEclassTest.pl file.fa

5.  結果はひとつのディレクトリにまとめて書き出される
    -   `file.fa`: 元ファイルからTEだけ抜き出したもの？
    -   `file.fa.html`: 一覧
    -   `file.fa.lib`: RepeatMasker用？
    -   `file.fa.stat`: LTRなどがそれぞれいくつあったか集計
