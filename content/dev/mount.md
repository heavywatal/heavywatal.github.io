+++
title = 'mount'
tags = ["communication"]
[menu.main]
  parent = "dev"
+++

<http://manpages.ubuntu.com/manpages/precise/man8/mount.8.html>

## 基本的な使い方

予め作っておいた空のディレクトリ `dir` に `device` をマウントする:

    % [sudo] mount [-t type] [-o option[,option]...] device dir

マウントを解除する:

    % [sudo] umount dir

既にマウントされているものを列挙:

    % mount -l

### 主な `-o` オプション

`defaults`
:   `rw,suid,dev,exec,auto,nouser,async` と同じ

`rw`
:   読み書きモード。逆は `ro`

`suid`
:   SUIDとSGIDを有効に

`dev`
:   マウントしたファイルシステム上にあるデバイスを使えるように

`exec`
:   プログラム実行を許可

`auto`
:   `mount -a` で一緒にマウントされるように

`nouser`
:   `root` 以外のユーザでマウントできないように

`async`
:   とりあえずメモリに置いたらいいことにして処理を進め、裏でディスクに書き込む。逆は `sync`

`nounix`
:   Unix拡張機能を無効に (cifs)

`iocharset=utf8`
:   文字コードの設定。デフォルトは `iso8859-1`

`uid`, `gid`
:   ユーザID、グループIDを指定

## cifs/samba

### Ubuntu 12.04 から cifs マウント

<http://manpages.ubuntu.com/manpages/precise/man8/mount.cifs.8.html>

1.  cifs-utils をインストール:

        % sudo apt-get install cifs-utils

2.  パスワードをコマンド履歴や `/etc/fstab` に残さなくて済むように
    `~/.cifs` のようなファイルを作っておく:

        username=iwasaki
        password=******

3.  `mount` コマンドでマウント:

        % sudo mount -t cifs -o defaults,iocharset=utf8,nounix,uid=$(id -u),gid=$(id -g),credentials=$HOME/.cifs //ADDRESS/VOLUME ~/mnt

4.  起動時に自動でマウントさせるには `/etc/fstab` に追記:

        //ADDRESS/VOLUME /home/iwasaki/mnt cifs credentials=/home/iwasaki/.cifs,uid=iwasaki,gid=iwasaki,nounix,iocharset=utf8,defaults 0 0

### Ubuntu 12.04 のホームディレクトリを cifs/smb マウント出来るようにする

1.  samba をインストール:

        % sudo apt-get install samba

2.  `/etc/samba/smb.conf` の一部を編集:

        [homes]
        comment = Home Directories
        browseable = no

        create mask = 0644

        directory mask = 0755

        valid users = %S

3.  サービスを再起動:

        % sudo service smbd restart

## afp

### Ubuntu 12.04 のコマンドラインから afp でマウント

1.  afpfs-ng-utils をダウンロードしてインストール:

        % wget http://launchpadlibrarian.net/90192653/afpfs-ng-utils_0.8.1-2_amd64.deb
        % sudo dpkg -i afpfs-ng-utils_0.8.1-2_amd64.deb

2.  以下のようなコマンドでマウント。できなかった:

        % mount_afp 'afp://user:password@address/volume/' ~/mnt

### Ubuntu 12.04 の Nautilus から afp でマウント

1.  Nautilusをアクティブにして `control + l`
    (あるいはメニューバーから `Go --> Location...`）
2.  Location に `afp://***.***.***.***` という形でIPアドレスを入力してConnect

### Mac の Finder からマウント

1.  Finderをアクティブにして `command + k`
    (あるいはメニューバーから `Go --> Connect to Server...`)
2.  Server Address に `afp://***.***.***.***` という形でIPアドレスを入力してConnect


## sshfs

https://github.com/libfuse/sshfs

リモートのファイルシステムを[ssh]({{< relref "ssh.md" >}})経由でマウントし、
Finderとかで普通のディレクトリのように扱えるようにする。
マウントされる側はsshでログインできさえすればいいので、
手元のマシンに必要な準備をする。

1.  Macなら[FUSE for macOS](https://osxfuse.github.io/)を、
    Linuxなら[libfuse](https://github.com/libfuse/libfuse)をインストール

1.  sshfsをインストール:
    ```sh
    % brew install sshfs
    ```

1.  マウントポイントにする適当なディレクトリを作る。
    e.g., `mkdir ~/mnt`

1.  マウントする:
    ```sh
    % sshfs watal@example.com:/home/watal ~/mnt
    ```

1.  アンマウントする:
    ```sh
    % umount ~/mnt
    ```
