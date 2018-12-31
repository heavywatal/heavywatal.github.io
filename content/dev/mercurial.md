+++
title = 'Mercurial'
subtitle = "分散型バージョン管理システム"
tags = ["vcs", "python"]
[menu.main]
  parent = "dev"
+++

-   [Official Website](https://www.mercurial-scm.org/)
-   [Wiki](https://www.mercurial-scm.org/wiki/)
-   [日本語Wiki](https://www.mercurial-scm.org/wiki/JapaneseMercurial)
-   [The Definitive Guide by Bryan O'Sullivan](http://hgbook.red-bean.com/read/)

{{%div class="note"%}}
[GitHub](https://github.com)が確固たる地位を確立し、
今や[BitBucket](https://bitbucket.org)もGitに対応したので、
どうしても既存のMercurialリポジトリを使わなきゃいけない場合を除いて、
基本的にはGitを使うようにしたほうがよさそう。
See [Git]({{< relref "git.md" >}}).
{{%/div%}}

## 基本的な操作

ファイルの操作。`hg add` はWorking directoryにあるファイルをRepositoryに追加する。
そのほか `rm`, `mv`, `cp` などはコマンド直打ちではなく
`hg` を介して行うことにする。
`hg addremove` はWorking directoryに有るファイルを全部 `hg add` して、
無いファイルを全部 `hg rm` する。
これらの変更は `hg commit` するまでRepositoryには反映されない:

    hg add [FILE]...
    hg rm FILE
    hg mv SOURCE... DEST
    hg cp SOURCE... DEST
    hg addremove

Working directoryにおける（まだ `hg commit` されてない）変更を確認:

    hg status

Working directoryにおける変更をRepositoryに反映させる。
引数でファイルを指定できて、省略すると全部。
`-m` オプションでメッセージを指定しない場合は
`$EDITOR` が起動してコメントを求められるので、何か書いて保存、終了:

    hg commit -m "a message that describes the modifications you made"

ssh越しでRepositoryをやり取り。
ディレクトリの指定方法が `scp` とはちょっと違う。
スラッシュ1つだと `~/` 、
スラッシュ2つだとルートからの絶対パス。
`~/.ssh/config` などでちゃんと設定しとけばURLは簡略化可能。
cf. [ssh]({{< relref "ssh.md" >}}):

    hg push ssh://username@example.com//home/username/the_project
    hg push ssh://username@example.com/the_project
    hg pull ssh://username@example.com/the_project

`hg push` や `hg pull` はRepositoryの情報を送受信するだけで、
受け手のWorking directoryを変更しない。
受け取った側が `hg update` した時点で変更が適用される:

    hg update

## プロジェクト開始

<http://mercurial.selenic.com/wiki/QuickStart>

既にあるプロジェクトを取ってくる。`DEST` 省略時は元と同じ名前でディレクトリが作られる:

    hg clone [OPTION]... SOURCE [DEST]
    hg clone http://selenic.com/hg mercurial-repo
    cd mercurial-repo
    hg parents

新しいMercurialプロジェクトを開始する

1.  プロジェクトのルートディレクトリ（無ければ `mkdir` するなどして）に入って初期化:

        cd the_project/
        hg init
        ls -a
        ./ ../ .hg/

1.  プロジェクト固有の設定を `the_project/.hg/hgrc` に記述。
    例えば以下のように書いておけば `pull/push` の対象を省略できる。:

        [paths]
        default = ssh://username@example.com/the_project

1.  一時ファイルやバイナリファイルを無視するように、
    除外設定を `the_project/.hgignore` に記述
1.  除外設定が正しく効いてるか確認:

        hg status

1.  リポジトリにファイルを追加してコミット、確認。
    `add` は個別にファイルを指定できて、省略すると全部。:

        hg add
        hg commit -m "first commit"
        hg parents

## よく使うコマンド

あのrevisionではどんな変更したっけ？:

    hg diff -c 42

あのrevisionから今までにどこが変わった？:

    hg diff -r 42

いろいろやってみたけど今回の変更を全部無かったことにする
(`hg commit` する前):

    hg revert --all --no-backup

直前の `hg commit` を取り消す
(`hg push` する前の1回分のみ有効):

    hg rollback

`hg push` 済みあるいは複数回 `hg commit` してしまった後、変更を取り消す:

    hg backout -r 42
    hg commit

管理対象外のファイルを確認・削除する (`purge` extentionを有効にして):

    hg clean -p
    hg clean

## Merge

1.  `hg heads` でマージすべき2つの頭を確認
1.  `hg merge` でマージ実行
    -   これだけで解決できたら次のステップに
    -   conflictが生じた場合の挙動は設定によって大きく異る。
        以下のような設定にしておくと:

            [ui]
            merge = internal:merge

        Mercurialが出来る限りのマージをして、
        できなかった部分にマークのつけて返してくれるので、
        そのファイルを自分のエディタで編集して、次に進む

    -   `vimdiff` や `ediff` を使う設定になっていると
        `vi` とか `emacs` が起動する。
        以下は `ediff` の説明。
        -   左上(a)は現在の頭、右上(b)が別の頭、下(c)が結果
        -   `?` でヘルプ表示
        -   `n` と `p` で衝突箇所を移動
        -   衝突箇所ごとに `a` と `b`
            を押してどっちの版を採用するか決めていく。
        -   `wc` して結果を書き込む

1.  解決の必要なファイルを確認:

        hg resolve -a

1.  解決済みであることをマークしてコミット:

        hg resolve -m some_source.py
        hg commit

## 設定

<http://www.selenic.com/mercurial/hgrc.5.html>

以下のリストの上から順に探して読んでって、どんどん上書きしていく。すなわち下のやつほど優先順位が高い:

    <install-root>/etc/mercurial/hgrc.d/*.rc
    <install-root>/etc/mercurial/hgrc
    /etc/mercurial/hgrc.d/*.rc
    /etc/mercurial/hgrc
    ~/.hgrc
    <repo>/.hg/hgrc

ユーザー設定は `~/.hgrc` で。
`username` は `commit` するときに使われる:

    [ui]
    username = Jean Sibelius <username@example.com>
    ignore = ~/.hgignore

`push` や `pull` の受け取り側の `.hg/hgrc` に以下が記述されていると、
転送完了のあと自動的に `update` される:

    [hooks]
    changegroup = hg update >&2

シェルのようにエイリアス設定も可能:

    [alias]
    ll = glog --stat --limit 6
    rep = !$HG locate --print0 | xargs -0 grep `$@

便利な [Extension](http://mercurial.selenic.com/wiki/UsingExtensions)
もここで設定:

    [extensions]
    color =
    fetch =
    graphlog =
    pager =
    purge =
    schemes =

    [pager]
    pager = LESS='-R' less
    attend = help, diff, log, glog, annotate

除外設定はWorking directory直下の `.hgignore` に記述。
あるいは上記のようにユーザーレベルでも指定できる。
正規表現とグロブの2つの表記法がある。
cf. <http://mercurial.selenic.com/wiki/.hgignore> :

    syntax: regexp
    ._
    .DS_Store
    \.out$`
    \.o$
    \.pyc$
    ~$

    syntax: glob
    ._*
    .DS_Store
    *.out
    *.o
    *.pyc
    *~

## インストール

PyPIに登録されてるPythonパッケージなので
[pip]({{< relref "pip.md" >}}) でインストールできる:

    pip install mercurial

でもPythonから `import` して使うことは無いので、
Macなら [Homebrew]({{< relref "homebrew.md" >}}) で入れちゃうほうが管理が楽チン:

    brew install mercurial

Linuxでソースからインストールしたい場合は
`python-devel` 的なパッケージを入れた上で:

    wget -O- http://mercurial.selenic.com/release/mercurial-3.0.tar.gz | tar xz
    cd mercurial-2.9.2/
    sudo make install-bin

`apt-get` や `yum` でもインストールできるが
たいてい公式リポジトリのやつはバージョンが古すぎてダメ。

## Gitに移行する

<https://github.com/frej/fast-export>

1.  `fast-export` をダウンロード:

        git clone https://github.com/frej/fast-export.git

1.  移行先のローカルリポジトリを作成:

        mkdir dst_git
        cd dst_git
        git init

1.  実行:

        path/to/fast-export/hg-fast-export.sh -r path/to/src_hg

1.  作業ディレクトリに反映:

        git status
        git checkout master
