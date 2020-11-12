+++
title = '遺伝研スパコン'
tags = ["job"]
[menu.main]
  parent = "bio"
+++

https://sc.ddbj.nig.ac.jp/

## 利用開始

1.  [新規ユーザ登録申請](https://sc2.ddbj.nig.ac.jp/index.php/ja-new-application)
1.  責任者にメールが届くので、それに従って誓約書PDFを管理者に送信
1.  アカウント登録証が手元に届く
1.  https://sc2.ddbj.nig.ac.jp/ でログイン
1.  手元のマシンでSSH公開鍵を生成し、
    [サーバーに登録](https://sc2.ddbj.nig.ac.jp/index.php/2014-09-17-05-42-33)。
    ドキュメントの例ではRSAが使われてるけどECDSAのほうが安全で軽い。
    ほんとはed25519のほうが良いけどウェブ登録ではなぜか拒否される。

    ```sh
    ssh-keygen -t ecdsa -b 521 -N '' -C 'heavywatal@nig.ac.jp' -f ~/.ssh/id_ecdsa_nig
    ```

1.  手元の `~/.ssh/config` に設定を追加:

        Host *.ddbj.nig.ac.jp
          User heavywatal
          RequestTTY yes
          IdentityFile ~/.ssh/id_ecdsa_nig

    sshコマンドでユーザ名と `-t` オプションを省略できるようになる。

1.  ゲートウェイノードにSSH接続してログインノードに `qlogin`:

        ssh gw.ddbj.nig.ac.jp qlogin


## 環境整備

- [ハードウェア構成](https://sc.ddbj.nig.ac.jp/ja/guide/hardware)
- [ソフトウェア構成](https://sc.ddbj.nig.ac.jp/ja/guide/software)
  (Phase 3: Red Hat Enterprise Linux 7.5)
    - [Singularity](https://sc.ddbj.nig.ac.jp/ja/guide/software/singularity)
    - [Environment Module](https://sc.ddbj.nig.ac.jp/ja/guide/software/environmental-modules)


## ファイルの送受信

[公式「システムへのファイル転送方法」](https://sc.ddbj.nig.ac.jp/ja/guide/software/file-transfer)
にはsftpかAsperaを使えと書かれてるけど、
[rsync]({{< relref "rsync.md" >}})
を使うのが簡単かつ効率的。

```sh
# send
rsync -auv ~/input/ gw.ddbj.nig.ac.jp:~/input/

# receive
rsync -auv gw.ddbj.nig.ac.jp:~/output/ ~/output/
```

ソースコードは当然[Git]({{< relref "git.md" >}})で管理。


## ジョブ投入、管理

- <https://sc.ddbj.nig.ac.jp/ja/guide/software/univa-grid-engine>
- <https://sc2.ddbj.nig.ac.jp/index.php/ja-howtouse>
- <https://sc2.ddbj.nig.ac.jp/index.php/ja-uge-additional>
- [Univa Grid Engine (UGE)](http://www.univa.com/products/grid-engine)
- <http://gridengine.eu/grid-engine-documentation>

### `qsub`

遺伝研ウェブサイトにはスクリプトを書いてから渡す方法しか書いてないが、
簡単なタスクならコマンド引数として直接渡すのもあり。

```sh
## スクリプトを書いてから渡す
qsub test.sh

## コマンドライン引数として直接渡す
qsub -l short -b y -shell n -cwd -N test "pwd; sleep 5; ls >ls.txt"
```

`-l ***`
:   実行時間や計算ノードなどの要求を伝える。
    管理者が定義したキューから選んで指定する。
    例えば、3日以内に終わる軽いものなら `-l short`、
    2か月かかるものなら `-l epyc`、
    メモリが多めに必要なら `-l medium`、など。
    `qstat -g c` でキューの一覧とそれぞれの負荷を確認できるので空いてるところを探す。
:   1コアあたりのRAM上限(デフォルト8GB)もこのオプションから
    `-l s_vmem=16G -l mem_req=16G` のように変更できる。

`-pe def_slot 8`
:   parallel environment:
    並列処理で使用するCPUコア数を指定。
:   MPIによる並列化の場合はまた違うオプションがある。

`-cwd`
:   カレントディレクトリでジョブ実行。
    デフォルトでは `${HOME}`。

`-N ***`
:   ジョブに名前をつける。
    デフォルトではスクリプト名が採用される。

`-o ***`, `-e ***`
:   標準出力・標準エラー出力の書き出し先。
    デフォルトではワーキングディレクトリ以下に
    `{JOBNAME}.o{JOBID}`, `{JOBNAME}.e{JOBID}`
    という名前で書き出される(空っぽでさえ)。
    不要な出力は `/dev/null` に流し込むべし。

`-S /bin/sh`
:   インタープリタを指定。
    指定しないと `csh` が利用されるせいか、
    標準出力で `Warning: no access to tty` と怒られる。

`-t 1-M`
:   M個のタスクを持つアレイジョブとして投入する。
    タスクごとにパラメータを変えてプログラムを走らせたい場合は、
    スクリプトの中で環境変数 `SGE_TASK_ID` を拾ってどうにかする。
    `qsub` コマンドを生成して何回も呼ぶほうが圧倒的に楽だけど、
    アレイジョブのほうがクラスタ側の負荷が小さいらしい。ほんとかな。

`-tc MAX_RUNNING_TASKS`
:   アレイジョブで同時実行するタスク数の上限を指定。
    システムからユーザーに与えられた上限は `qquota` コマンドで確認できる。
    いまのところ500らしい。

`-v VARIABLE=value`
:   環境変数を定義してジョブに引き継ぐ。
:   大文字 `-V` で全ての環境変数が渡される。

`-b y`
:   計算ノードにバイナリがあるものとしてジョブを投げる。
    これを指定しない場合はスクリプト扱いになり、
    投入ノードから計算ノードへのコピーなど余計なプロセスが挟まるらしい。

`-shell n`
:   環境変数の解決など、
    プログラムの呼び出しにシェルを介す必要がない場合は
    これを指定することで多少コスト削減できる。
    当然 `-b y` のときのみ有効。


### `qsub` スクリプト

`#$` で始まる行は `qsub` へのオプションと見なされる。

ジョブスクリプト内で参照可能な特殊環境変数をプリントしてみるジョブの例:

```sh
#!/bin/sh
#$ -S /bin/sh
#$ -l short
#$ -cwd
#$ -t 1-2
echo HOME: $HOME
echo USER: $USER
echo JOB_ID: $JOB_ID
echo JOB_NAME: $JOB_NAME
echo HOSTNAME: $HOSTNAME
echo SGE_TASK_ID: $SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST: $SGE_TASK_LAST
echo SGE_TASK_STEPSIZE: $SGE_TASK_STEPSIZE
pwd
ls
```

スクリプトはPythonでもいい。
インタープリタを `-S /usr/bin/env python` で指定できないのは残念。

```py
#!/usr/bin/env python
#$ -S $HOME/.pyenv/shims/python
#$ -l short
#$ -cwd
#$ -t 1-2
import os
print("SGE_TASK_ID: " + os.environ["SGE_TASK_ID"])
```

### 補助コマンド

`qstat`
:   現在実行中のジョブ一覧。
    ステータスは1文字に省略されて
    E(rror), r(unning), R(estarted), s(uspended)
    のように表示される。
    詳しくは `man qstat` を参照。

`qstat -g c`
:   クラスタで定義されているキューの一覧と、それぞれの負荷を表示

`qstat -f | less`
:   全ノードの状況をfullに表示

`qstat -u '*' | less`
:   全ユーザのジョブを表示。
    `-s p` でpending中のみに絞ったり、
    `-l medium` でキューの種類を絞ったりできる。

`qstat -j JOBID`
:   ジョブの詳細表示

`qacct -j JOBID`
:   ジョブ実行後にリソース消費を確認

`qdel JOBID`
:   ジョブ削除

`qquota`
:   ユーザーに与えられたリソースを表示


## [Singularity](https://sc.ddbj.nig.ac.jp/ja/guide/software/singularity)

`/usr/local/biotools/` 以下に各種ソフトウェアが用意されている。
[BioContainers](https://biocontainers.pro) のものをほぼそのまま置いているらしい。

利用可能なソフトウェアとバージョンを探す:
```sh
find /usr/local/biotools/ -name 'blast*' | sort
```

イメージとプログラム名を指定して実行:
```sh
singularity exec -e /usr/local/biotools/f/fastp:0.20.0--hdbcaa40_0 fastp --help
```

そこらに落ちてるイメージを拾ってきて使うこともできる。
[例えばTrinityの公式最新版を使いたい場合](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-in-Docker#trinity_singularity):
```sh
find /usr/local/biotools/ -name 'trinity*' | sort
mkdir -p ~/image
wget -P ~/image/ https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/trinityrnaseq.v2.11.0.simg
singularity exec -e ~/image/trinityrnaseq.v2.11.0.simg Trinity --help
```


## misc.

### インストール済みRを使ってみる

https://sc.ddbj.nig.ac.jp/ja/guide/software/r

`module load r/3.5.2` で比較的新しいやつが使えるけど、
古いコンパイラ(おそらく `/usr/bin/gcc` 4.8.5)でビルドされているため
RcppでC++11までしか使えない。
また、各パッケージも同じく古いコンパイラでビルドしなければならない。
`module load gcc` などで新しいgcc/g++がPATH上に乗っていると、
Rcppインストール時などに
<code>/usr/lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found</code>
と怒られる。


### 自前でRをインストールする

`/opt/pkg/r/*/lib64/R/etc/Makeconf`
を参考に新しいコンパイラで自前ビルドを試みる:

```sh
wget -O- https://cran.r-project.org/src/base/R-3/R-3.5.2.tar.gz | tar xz
cd R-3.5.2/
./configure -h | less
./configure --prefix=${HOME}/R --disable-openmp --disable-java '--enable-R-shlib' '--enable-shared' '--with-tcl-config=/usr/lib64/tclConfig.sh' '--with-tk-config=/usr/lib64/tkConfig.sh' 'PKG_CONFIG_PATH=/usr/local/lib64/pkgconfig:/usr/local/lib/pkgconfig:/usr/local/share/pkgconfig:/cm/local/apps/curl/lib/pkgconfig:/usr/lib64/pkgconfig:/usr/lib/pkgconfig:/usr/share/pkgconfig' 'CFLAGS=-I/usr/local/include:/usr/include/X11:/cm/local/apps/curl/include' 'CPPFLAGS=-I/usr/local/include:/usr/include/X11:/cm/local/apps/curl/include' 'CXXFLAGS=-I/usr/local/include:/usr/include/X11:/cm/local/apps/curl/include' 'FFLAGS=-I/usr/local/include:/usr/include/X11:/cm/local/apps/curl/include' 'FCFLAGS=-I/usr/local/include:/usr/include/X11:/cm/local/apps/curl/include' 'LDFLAGS=-L/usr/local/lib64 -L/usr/lib64 -L/usr/lib -L/cm/local/apps/curl/lib'
make -j2
make install
```

`configure: error: libcurl >= 7.22.0 library and headers are required with support for https`
古いcurlは入ってないように見えるのに、なぜ？
