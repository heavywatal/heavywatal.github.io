+++
title = 'rsync'
tags = ["communication"]
[menu.main]
  parent = "dev"
+++

ファイルをコピーし、2つのディレクトリを同期する。
更新があったものだけをコピーする、
ひとつのsshセッションで複数のファイルを送受信する、
という使い方が可能なので `cp` や `scp` よりも便利な場合が多い。

## 基本

オプションについては後述するとして、基本は:

    rsync -auv SRC/ DST/

SRC側の末尾のスラッシュの有無によって結果が大きく異なることに注意。
DST側は付けても付けなくても同じ。

```sh
## SRC/ 以下のファイルが DST/ 以下に入る
rsync -auv SRC/ DST

## ディレクトリ DST/SRC が作られる
rsync -auv SRC DST

## 結果は同じだが、下のほうがより明示的
rsync -auv SRC/DIR DST
rsync -auv SRC/DIR/ DST/DIR/
```

[sshの設定]({{< relref "ssh.md" >}})をしておけばリモートホストへの転送も可能。
その場合は宛先を `remote_machine:~/DST` のようにコロンで指定する。

## Options

`-a, --archive`
:   バックアップ用途のオプション一式 `-rlptgoD`

`-u, --update`
:   受け手の方が新しいファイルをスキップ

`-v, --verbose`
:   冗長なメッセージ表示

`-n, --dry-run`
:   実際に送受信を行わず試してみる

`-z, --compress`
:   通信量を減らしたければ

`--delete`
:   SRC側に存在しないものがDST側にあったとき削除

`--delete-excluded`
:   除外設定されているファイルが受け手側にあったら削除（危険！）

`--ignore-existing`
:   受信側に存在していたら無視

## Exclude and include

`--include=<PATTERN>`
:   マッチするファイル・ディレクトリを除外しない

`--exclude=<PATTERN>`
:   マッチするファイル・ディレクトリを除外

`--exclude-from=<FILE>`
:   ファイルに記述した除外パターンを読む

先に記述したものほど優先される。

1.  いつでも除外したいものを `$HOME/.rsync/exclude` の中に記述する:

        ._*
        .DS_Store
        .Trash
        .Trashes
        .Spotlight-*
        .hidden
        .vol
        .localized
        *~
        *.o
        *.out
        *.pyc
        *.zhistory
        known_hosts

2.  `--exclude-from` オプションでそのファイルを読ませる。
    例えば `.zshrc` にこう書く:

        alias rsync='rsync --exclude-from=$HOME/.rsync/exclude'

3.  そのほかで除外したいものは `--exclude` オプションで個別に指定
