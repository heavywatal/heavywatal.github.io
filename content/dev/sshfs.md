+++
title = 'sshfs'
tags = ["communication"]
[menu.main]
  parent = "dev"
+++

SSH File Transfer Protocol

リモートのファイルシステムを [ssh]({{< relref "ssh.md" >}}) 経由でローカルにマウントする。
必要なライブラリや導入方法はOS毎に異なる。
マウントにはFUSEというライブラリを利用している。

[<http://fuse.sourceforge.net/sshfs.html>](http://fuse.sourceforge.net/sshfs.html)

## GUI: Mac

1.  [OSXFUSE](http://osxfuse.github.com/) をインストール
2.  [Macfusion](http://macfusionapp.org/) をインストール
3.

    Macfusion を起動
    :   1.  左下の[+]からSSHFSを選択
        2.  [ssh]({{< relref "ssh.md" >}}) のconfigをちゃんとしておけばここでの設定はとりあえずホスト名とかユーザ名だけよい
        3.  圧縮、シムリンク、Mac固有の不可視ファイルなどをどうするかは好みで
        4.  Mount

## CUI: Linux, Mac

1.  [ssh]({{< relref "ssh.md" >}}) の設定を済ます
    （公開鍵によるパスワード無し認証、configによるユーザ名やホスト名などの設定）
2.  sshfs をインストール。 (Synaptic でもよい):

        % sudo apt-get install sshfs
        or
        % sudo port install sshfs

3.  マウントポイントにする適当なディレクトリを作る。(`$ mkdir /mnt/meme`)
4.  マウントする。[ssh]({{< relref "ssh.md" >}}) のconfigをしちゃんとしとけば下段のように省略できる。:

        % sshfs USER@HOST:DIR MOUNTPOINT
        or
        % sshfs meme: /mnt/meme

5.  アンマウントする:

        % fusermount -u /mnt/meme

## Windows

FUSEのWindows版に相当するのが [Dokan](http://dokan-dev.net/download/) 。
なんか、俺の環境だと遅くてダメだった。

1.  Microsoft Visual C++ 2005 SP1をインストール。(2008ではダメだった)
2.  Dokanライブラリをインストール。
3.  Dokan SSHFSをインストール。
4.  起動して、アドレスや鍵を指定してマウントする。
5.  タスクバーのアイコンからUnmountしたりExitしたりする。
