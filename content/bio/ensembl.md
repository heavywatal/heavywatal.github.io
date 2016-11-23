+++
title = 'Ensembl'
tags = ["genetics", "database"]
[menu.main]
  parent = "bio"
+++

<http://www.ensembl.org/>

よくアノテーションされてて、相同な遺伝子などを探しやすいゲノムデータベース。
[RefSeq](http://www.ncbi.nlm.nih.gov/refseq/) とか [UniProt](http://www.uniprot.org/) に比べると、ゲノム上での位置とかが得意なのかな。
脊椎動物以外は別サイトに分かれてる。

-   <http://metazoa.ensembl.org/>
-   <http://plants.ensembl.org/>
-   <http://fungi.ensembl.org/>
-   <http://protists.ensembl.org/>
-   <http://bacteria.ensembl.org/>

{{%div class="note"%}}
Features の Attributes
でEnsemblやRefSeqのIDを複数同時に指定すると
パラログの対応がちゃんとしてないせいか
ID同士の組み合わせの数だけ結果が膨れてしまう...?
{{%/div%}}

## プログラムからアクセス

[Ensembl API](http://www.ensembl.org/info/docs/api/)
:   Perl 限定なので却下

[REST API](http://beta.rest.ensembl.org/)
:   2012年9月に公開され、Perl以外の言語からでも簡単にアクセスできるようになった。
    例えば Python ラッパーとして
    [pyensemblrest](https://pypi.python.org/pypi/pyensemblrest)
    というパッケージが公開されている。

[biomaRt]({{< relref "rstats/biomart.md" >}})
:   R の [Bioconductor]({{< relref "rstats/bioconductor.md" >}}) 上で利用できる。
    データが少なければこれが一番使い勝手いいかなぁ。
    単なるダウンロードではなくサーバー上の演算が絡むっぽいので
    データが多くなるとかなり遅くなる。

[FTP Download](http://www.ensembl.org/info/data/ftp/index.html)
:   どんなデータが入手可能かブラウザから見てみるのには良い。

[rsync](http://www.ensembl.org/info/data/ftp/rsync.html)
:   複雑なパターンマッチや、既存ファイル無視などができる分、
    実際に落とすときはFTPよりも `rsync` のほうが強力。 cf. [/dev/rsync]({{< relref "dev/rsync.md" >}})

[MySQL](http://www.ensembl.org/info/data/mysql.html)
:   データベースを構成する MySQL にダイレクトにアクセスする。
    下記のようにローカルサーバーにごっそりミラーリングすることも可能。

## ローカルデータベース構築

イトヨの例

1.  MySQL 環境を整える
2.  データをダウンロードする:

        % x=gasterosteus_aculeatus_core_57_1j
        % wget -m ftp://ftp.ensembl.org/pub/current_mysql/$x

3.  落としてきたファイルをデータベースに組み込む:

        % cd ftp.ensembl.org/pub/current_mysql
        % gunzip $x/*.gz
        % mysqladmin5 -u root --password=PASSWORD create $x
        % mysql5 -u root --password=PASSWORD $i < $x/$x.sql
        % mysqlimport5 -u root --password=PASSWORD $i -L $x/*.txt

    全種でやるなら `wget -m ftp://ftp.ensembl.org/pub/current_mysql`
    で全部落として、以下のようなシェルスクリプトを実行する:

        for x in *; do
            echo $x
            gunzip $x/*.gz
            mysqladmin5 -u root --password=PASSWORD create $x
            mysql5 -u root --password=PASSWORD $x < $x/$x.sql
            mysqlimport5 -u root --password=PASSWORD $x -L $x/*.txt
        done

スペースとスピードの観点からすると、本当はいちいち展開したファイルを作らずに
`gunzip -c dna.txt.gz | mysqlimport ...`
ってな感じでパイプで流し込みたい。
["mysqlimport stdin"でググってみる](http://www.google.co.jp/search?q=mysqlimport+stdin)
と、どうやら同じような不満を持つ人がたくさんいるらしいことと、
MySQL側はそれを改善する気がないことと、
Perlでは `MySQL::Slurp` を使ってそれができるということはわかった。
