+++
title = '遺伝研スパコン'
tags = ["job"]
[menu.main]
  parent = "bio"
+++

[Homepage](https://sc.ddbj.nig.ac.jp/)

## 利用開始

<https://sc.ddbj.nig.ac.jp/index.php/ja-howtouse>

1.  [新規ユーザ登録申請](https://sc.ddbj.nig.ac.jp/index.php/ja-new-application)
2.  責任者にメールが届くので、それに従って誓約書PDFを管理者に送信
3.  アカウント登録証が手元に届く
4.  Homepage\_ でログイン
5.  [SSH公開鍵登録](https://sc.ddbj.nig.ac.jp/index.php/2014-09-17-05-42-33) (パスワード認証ではSSHできない)
6.  `~/.ssh/config` に設定を追加:

        Host *.ddbj.nig.ac.jp
          User heavywatal
          RequestTTY yes

    sshコマンドでユーザ名と `-t` オプションを省略できるようになる。

7.  システム内ノード間SSHのパスワード認証を省くために
    非公開鍵も送っておく(いいのかな？):

        % rsync -auv ~/.ssh/id_rsa gw2.ddbj.nig.ac.jp:~/.ssh/

8.  ゲートウェイノードにSSH接続してログインノードに `qlogin`:

        % ssh gw2.ddbj.nig.ac.jp qlogin

{{%div class="note"%}}
Phase 1 と Phase 2 という異なるシステムが存在しているが、
基本的にはPhase 2システムを使えば良さそう。
共有されるのはユーザ情報のみで、ホームのファイルシステムも別。
Phase 1システムにデータをコピーするには、
Phase 2ゲートウェイから
`qlogin -l trans`
でデータ移行用ノードにログインし、
`/home_ph1/` 以下に見える自分のホームに `rsync` する。
{{%/div%}}

## 環境整備

<https://sc.ddbj.nig.ac.jp/index.php/system-software-config>

Phase 1 Red Hat Enterprise Linux 6.1\
Phase 2 Red Hat Enterprise Linux 6.4

Linuxbrewで環境を整える。
cf. [/dev/devenv]({{< relref "dev/devenv.md" >}})

{{%div class="note"%}}
mercurialはインストールできてもglibcらへんの関係でうまく動かない。
emacsは `User *** has no home directory` という謎のエラーを吐く。
clang/boostは要調整。
{{%/div%}}

### ログインシェルをzshに変更(しないで対処)

<https://sc.ddbj.nig.ac.jp/index.php/ja-tips>

LDAPで管理されているので `chsh` は効かない:

    % ldapmodify -x -D uid=heavywatal,ou=people,dc=nig,dc=ac,dc=jp -W
    Enter LDAP Password: ********
    dn: uid=heavywatal,ou=people,dc=nig,dc=ac,dc=jp
    changetype:modify
    replace:loginShell
    loginShell:/bin/zsh
    # ctrl-d

反映されるまで少し時間がかかるっぽい。

`/etc/profile.d/*` の読み込みのタイミングのせいか、
zshにすると `ssh -t gw.nig qlogin` で `command not found` になってしまう。
しかも `/bin/zsh` が結構古くて微妙。

`~/.bashrc` にエイリアスを定義して対処したほうが良さそう:

    PATH=${HOME}/.homebrew/bin:${PATH}
    alias zmux='SHELL=$(brew --prefix)/bin/zsh tmux'

ちなみにtmuxセッションの寿命はどうなってるんだろう...

### 最新版Rをインストール

環境が古すぎて手こずった。
依存関係について、何が必要十分なのかまでは検証していないが、
いくつかのライブラリは新しいものをLinuxbrewで入れておく必要がある。
e.g., gcc, binutils, bzip2

```sh
wget -O- https://cran.r-project.org/src/base/R-3/R-3.3.0.tar.gz | tar xz
cd R-3.3.0/
./configure -h | less
./configure LIBS="-lpthread" --prefix=${HOME}/R --disable-openmp --disable-java
make
make install
```

## ジョブ投入、管理

<https://sc.ddbj.nig.ac.jp/index.php/ja-howtouse>

<https://sc.ddbj.nig.ac.jp/index.php/ja-uge-additional>

[Univa Grid Engine (UGE)](http://www.univa.com/products/grid-engine)

<http://gridengine.eu/grid-engine-documentation>

### `qsub`

遺伝研ウェブサイトにはスクリプトを書いてから渡す方法しか書いてないが、
コマンド引数として直接渡すほうが圧倒的に楽チン。

```sh
## スクリプトを書いてから渡す
qsub test.sh

## コマンドライン引数として直接渡す
qsub -l debug -b y -shell n -cwd -N test "pwd; sleep 5; ls >ls.txt"
```

`-l ***`
:   管理者が定義したキューから選んで指定し、
    実行時間や計算ノードなどの要求を伝える。
    例えば、24時間以内に終わるものなら `-l short` や `-l debug`、
    2か月かかるものなら `-l month` を指定。
   `qstat -g c` でキューの一覧とそれぞれの負荷を確認できる。
:   1スレッドあたりのRAM上限(デフォルト4GB)もこのオプションから
    `-l s_vmem=8G -l mem_req=8G` のように変更する。

`-cwd`
:   カレントディレクトリでジョブ実行。
    デフォルトでは `$HOME`。

`-N ***`
:   ジョブに名前をつける。
    デフォルトではスクリプト名。

`-o ***`, `-e ***`
:   標準出力・標準エラー出力の書き出し先。
    デフォルトではワーキングディレクトリ以下に
    `{JOBNAME}.o{JOBID}`, `{JOBNAME}.e{JOBID}`
    という名前で書き出される(空っぽでさえ)。
    不要な出力は `/dev/null` に流し込むべし。

`-b y`
:   計算ノードにバイナリがあるものとしてジョブを投げる。
    これを指定しない場合はスクリプト扱いになり、
    投入ノードから計算ノードへのコピーなど余計なプロセスが挟まるらしい。

`-shell n`
:   環境変数の解決など、
    プログラムの呼び出しにシェルを介す必要がない場合は
    これを指定することで多少コスト削減できる。
    当然 `-b y` のときのみ有効。

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
:   環境変数を定義してジョブに引き継ぐ

`-pe def_slot N`
:   parallel environment:
    並列処理で使用するスレッド数を指定。
    `-l mem_req=*G` の値は1スレッドあたりなので注意。

### `qsub` スクリプト

`#$` で始まる行は `qsub` へのオプションと見なされる。

ジョブスクリプト内で参照可能な特殊環境変数をプリントしてみるジョブの例:

```sh
#!/bin/sh
#$ -S /bin/sh
#$ -l debug
#$ -cwd
#$ -t 1-2
#$ -N test_sh
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
#$ -S $HOME/.virtualenv/py3/bin/python
#$ -l debug
#$ -cwd
#$ -t 1-2
#$ -N test_py
import os
print("SGE_TASK_ID: " + os.environ["SGE_TASK_ID"])
```

### 補助コマンド

`qstat`
:   現在実行中のジョブ一覧

`qstat -g c`
:   クラスタで定義されているキューの一覧と、それぞれの負荷を表示

`qstat -f`
:   全ノードの状況をfullに表示

`qstat -j JOBID`
:   ジョブの詳細表示

`qacct -j JOBID`
:   ジョブ実行後にリソース消費を確認

`qdel JOBID`
:   ジョブ削除

`qquota`
:   ユーザーに与えられたリソースを表示
