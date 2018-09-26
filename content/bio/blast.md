+++
title = 'BLAST'
tags = ["genetics"]
[menu.main]
  parent = "bio"
+++

Basic Local Alignment Search Tool

<http://blast.ncbi.nlm.nih.gov/>

## blast+ のインストール

[Homebrew/Linuxbrew]({{< relref "homebrew.md" >}}) を使って
`brew install blast` するのが楽チン。
以下に紹介するのはそれ以外の茨の道。

### Linux, Mac (binary)

1.  <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>
    から最新版の `ncbi-blast-*-x64-linux.tar.gz`
    あるいは `ncbi-blast-*-universal-macosx.tar.gz`
    をダウンロードして展開:

        % wget -O- ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.30+-x64-linux.tar.gz | tar xz

2.  しかるべきところに移動してシムリンク:

        % sudo mv ncbi-blast-2.2.30+ /usr/local/
        % sudo ln -s /usr/local/ncbi-blast-2.2.30+ /usr/local/ncbi-blast

    あとは `/usr/local/ncbi-blast/bin` にパスを通して使う。


### Mac (source)

`ncbi-blast-*.dmg` や `ncbi-blast-*-universal-macosx.tar.gz`
は古いGCCでUniversalビルドされてるっぽいので、
自分のマシンに最適化されるようにビルドしてみる。
ホントに速いかどうかは試してない。

<http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.Installation>

1.  Boost ライブラリをインストールする。このとき
    `regex`, `spirit`, `system`, `filesystem`, `test`, `thread`
    をビルド対象に含める。 cf. [Boost]({{< relref "boost.md" >}})
2.  <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>
    から最新版の `ncbi-blast-*-src.tar.gz` をダウンロード
3.  展開して移動:

        % tar xzf ncbi-blast-2.2.30+-src.tar.gz
        % cd ncbi-blast-2.2.30+-src/c++/

4.  `ncbi-blast-*-universal-macosx.tar.gz` に入ってる
    `ncbi_package_info` を参考に `configure` して `make`:

        % ./configure --without-debug --with-strip --without-pcre --with-mt --with-flat-makefile --with-64 --with-ncbi-public --without-ccache --without-caution --without-makefile-auto-update --with-projects=scripts/projects/blast/project.lst --with-internal --prefix=/usr/local/ncbi-blast --with-boost=/usr/local/boost-gcc CC=gcc-4.9 CXX=g++-4.9

    {{%div class="note"%}}
`clang` (および Xcode のニセ `gcc`)
では `configure` が通らないので本物の `gcc` を
[Homebrew]({{< relref "homebrew.md" >}}) などでインストールしておく。
    {{%/div%}}

5.  そのまま `make` してもダメらしいので
    `Makefile.mk` の `-std=gnu++11` を消す:

        % sed -i.orig -e 's/-std=[a-z0-9+]\+//' ReleaseMT/build/Makefile.mk

6.  ビルドしてインストール:

        % make
        % sudo make install

### 共通設定

1.  `.zshenv` 等でプログラムとデータベースのパスを設定:

        export PATH=/usr/local/ncbi-blast/bin:${PATH}
        export BLASTDB=${HOME}/db/blast

    {{%div class="note"%}}
`~/.ncbirc` でも `BLASTDB` を指定できるっぽいけど
普通にシェルの環境変数を使うほうがわかりやすいし、
`makeblastdb` の出力先として参照できるので便利。
    {{%/div%}}

2.  確認:

        % rehash
        % blastn -version

## ローカルデータベース構築

`-subject` オプションを使えばFASTAファイル同士でいきなり検索できるが、
何度もやるならデータベース化しておいたほうが効率いいはず。

<http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.makeblastdb_application_opt>

1.  データベースの置き場所をひとつ決め (e.g. `~/db/blast`)、
    上記のように環境変数 `BLASTDB` を設定しておく。
2.  FASTA形式の配列データを用意する (e.g. <ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/> から `yeast.aa.gz`)
3.  `makeblastdb` コマンドでBLASTデータベース形式に変換:

        % makeblastdb -in yeast.nt -out ${BLASTDB}/yeast.nt -dbtype nucl -parse_seqids -hash_index
        % makeblastdb -in yeast.aa -out ${BLASTDB}/yeast.aa -dbtype prot -parse_seqids -hash_index

    圧縮ファイルは直接読めないので、展開して標準入力に流し込む。
    このとき `-title` と `-out` は必須の引数になる:

        % gunzip -c mydata.fa.gz | makeblastdb -in - -title mydata -out ${BLASTDB}/mydata -dbtype nucl -hash_index

`-in` (stdin)

`-dbtype` (`prot`)
:   `prot` or `nucl`

`-title` (入力ファイル名)
:   どういうときに使われるか不明。
    標準入力を使う場合は必須オプション。

`-parse_seqids`
:   配列名が `gi|129295` みたいな特定の形式になってる場合にうまいことやる

`-hash_index`
:   よくわからんけど検索のスピードアップに繋がるっぽいオプション

`-out` (入力ファイル名)
:   出力先のベースネーム。
    標準入力を使う場合は必須オプション。

`-taxid`

## プログラム実行

### 共通オプション

<http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.options_common_to_all_blast>

`-db`

`-query` (stdin)

`-out` (stdout)

`-evalue` (10.0)
:   これを下回るE値を持つものだけ集める。

`-max_target_seqs` (500)
:   アラインされた配列をいくつまで保持して報告するか？
    [これを1にして得られる結果がベストヒットとは限らないことに注意。](https://doi.org/10.1093/bioinformatics/bty833)

`-subject`
:   データベース化してないFASTAファイルを検索対象として直接指定。
    ただし `-num_threads` を指定しても並列化できない。

`-num_threads` (1)

`-outfmt` (0)
:   0 = pairwise,\
    1 = query-anchored showing identities,\
    2 = query-anchored no identities,\
    3 = flat query-anchored, show identities,\
    4 = flat query-anchored, no identities,\
    5 = XML Blast output,\
    **6** = tabular,\
    **7** = tabular with comment lines,\
    8 = Text ASN.1,\
    9 = Binary ASN.1\
    **10** = Comma-separated values\
    11 = BLAST archive format (ASN.1)

    TSV/CSVテーブル(6, 7, 10)の場合は出力内容を細かく指定できる。デフォルトは
    `-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"`

    `qseqid`: Query Seq-id\
    `qgi`: Query GI\
    `qacc`: Query accesion\
    `sseqid`: Subject Seq-id\
    `sallseqid`: All subject Seq-id(s), separated by a ';'\
    `sgi`: Subject GI\
    `sallgi`: All subject GIs\
    `sacc`: Subject accession\
    `sallacc`: All subject accessions\
    `qstart`: Start of alignment in query\
    `qend`: End of alignment in query\
    `sstart`: Start of alignment in subject\
    `send`: End of alignment in subject\
    `qseq`: Aligned part of query sequence\
    `sseq`: Aligned part of subject sequence\
    `evalue`: Expect value\
    `bitscore`: Bit score\
    `score`: Raw score\
    `length`: Alignment length\
    `pident`: Percentage of identical matches\
    `nident`: Number of identical matches\
    `mismatch`: Number of mismatches\
    `positive`: Number of positive-scoring matches\
    `gapopen`: Number of gap openings\
    `gaps`: Total number of gap\
    `ppos`: Percentage of positive-scoring matches\
    `frames`: Query and subject frames separated by a '/'\
    `qframe`: Query frame\
    `sframe`: Subject frame\
    `btop`: Blast traceback operations (BTOP)\
    `staxids`: unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)\
    `sscinames`: unique Subject Scientific Name(s), separated by a ';'\
    `scomnames`: unique Subject Common Name(s), separated by a ';'\
    `sblastnames`: unique Subject Blast Name(s), separated by a ';' (in alphabetical order)\
    `sskingdoms`: unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)\
    `stitle`: Subject Title\
    `salltitles`: All Subject Title(s), separated by a '&lt;&gt;'\
    `sstrand`: Subject Strand\
    `qcovs`: Query Coverage Per Subject\
    `qcovhsp`: Query Coverage Per HSP

### `blastn`

Nucleotide query vs Nucleotide subject

<http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.blastn_application_options>

`-word_size` (11, short: 7, mega: 28)
:   最初にexact matchさせる配列の長さ。
    短いほどいろんな開始地点から探せるが、遅くなる。

`-gapopen` (5, mega: 0)

`-gapextend` (2, mega: none)

`-reward` (2, short: 1, mega: 1)

`-penalty` (-3, mega: -2)

`-perc_identity` (0)

`-ungapped`

------------------------------------------------------------------------

`blastp`
:   Protein query vs Protein subject
    <http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.blastp_application_options>

`blastx`
:   Nucleotide query (translated) vs Protein subject
    <http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.blastx_application_options>

`tblastn`
:   Protein query vs Nucleotide subject (translated)
    <http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.tblastn_application_options>

`tblastx`
:   Nucleotide query (translated) vs Nucleotide subject (translated)
    <http://www.ncbi.nlm.nih.gov/books/NBK279675/#appendices.T.tblastx_application_options>
