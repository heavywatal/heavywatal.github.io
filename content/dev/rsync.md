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
:   圧縮・展開のCPUコストはかかるけど通信量は減る

`--delete`
:   SRC側に存在しないものがDST側にあったとき削除

`--delete-excluded`
:   除外設定されているファイルが受け手側にあったら削除（危険！）

`--ignore-existing`
:   受信側に存在していたら無視

`--iconv`
:   文字コードを変換する。
    異なるOS/FS間でウムラウトや日本語を含むファイルを送受信するときに使う。
    リモート側だけ指定すれば十分だけどローカルも明示的に指定できる。
    その場合、pushかpullかによらず `={LOCAL},{REMOTE}` の順番で。
    例えばローカルのLinuxマシンとリモートの古いMacでやり取りする場合は
    `--iconv=utf8-mac` or `--iconv=utf8,utf8-mac` 。
    新しいMacのAPFSはLinuxと同じ `utf8` と見なしておけば良さそう...?


## Exclude and include

`--include=<PATTERN>`
:   マッチするファイル・ディレクトリを除外しない

`--exclude=<PATTERN>`
:   マッチするファイル・ディレクトリを除外

`--exclude-from=<FILE>`
:   ファイルに記述した除外パターンを読む

先に記述したものほど優先される。

1.  いつでも除外したいものを適当なファイル(`~/.config/rsync-exclude`)に記述する:

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

1.  `--exclude-from` オプションでそのファイルを読ませる。
    例えば `.zshrc` にこう書く:

        alias rsync='rsync --exclude-from=${HOME}/.config/rsync-exclude'

1.  そのほかで除外したいものは `--exclude` オプションで個別に指定


## SSH越しの送受信

```sh
rsync -auv user@example.com:~/dir/ ~/dir/
```

- [SSH公開鍵を設定]({{< relref "ssh.md" >}})してパスワード無しでログインできるようにしておく。
- `.ssh/config` でユーザー名とかも登録しておくとさらに楽。
  `RequestTTY yes` を付けてると当然ながら怒られる。
- リモート側の `.bashrc` とかで標準出力に何かを表示するようにしてあるとコケる
