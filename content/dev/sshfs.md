+++
title = 'sshfs'
tags = ["communication"]
[menu.main]
  parent = "dev"
+++

[SSH File Transfer Protocol](https://github.com/libfuse/sshfs)

リモートのファイルシステムをssh経由でマウントし、
Finderとかで普通のディレクトリのように扱えるようにする。
まずは、公開鍵によるパスワード無し認証、ユーザ名やホスト名のconfigなど、
[sshの設定]({{< relref "ssh.md" >}}) を済ませておく。

## GUI: Mac

1.  [OSXFUSE](https://osxfuse.github.io/) をインストール
1.  [Macfusion](http://macfusionapp.org/) をインストール
1.  Macfusion を起動
    1.  左下の[+]からSSHFSを選択
    2.  [ssh]({{< relref "ssh.md" >}}) のconfigをちゃんとしておけばここでの設定はとりあえずホスト名とかユーザ名だけよい
    3.  圧縮、シムリンク、Mac固有の不可視ファイルなどをどうするかは好みで
    4.  Mount

## CUI: Linux, Mac

1.  sshfs をインストール:
    ```sh
    % sudo apt-get install sshfs
    or
    % brew cask install sshfs
    ```

1.  マウントポイントにする適当なディレクトリを作る `mkdir ~/tmp/mnt`

1.  マウントする:
    ```sh
    % sshfs USER@HOST:DIR MOUNTPOINT
    % sshfs meme: ~/tmp/mnt/
    ```

1.  アンマウントする:
    ```sh
    % fusermount -u ~/tmp/mnt/
    or
    % umount ~/tmp/mnt
    ```
