+++
title = '/etc'
tags = ["shell"]
[menu.main]
  parent = "dev"
+++

## `/etc/default/`

### `/etc/default/ntpdate`

時刻合わせ。起動時に /usr/sbin/ntpdate-debian が動くことになっている。
その設定ファイルは `/etc/default/ntpdate`:

    NTPDATE_USE_NTP_CONF=no
    NTPSERVERS="130.34.11.117 130.34.48.32 ntp.nict.jp ntp.ubuntu.com"

## `/etc/hosts`

IPアドレスとホスト名のマッピングを行うファイル。
コマンドラインでIPを直打ちしなくて済むようになったり、
DNSサーバに問い合わせなくて済むようになったりする。

    [IP_ADDRESS] [FULL_NAME] [ALIAS] [ALIAS] ...

    127.0.0.1       localhost
    ::1             localhost

    130.34.107.135 meme.biology.tohoku.ac.jp meme

## `/etc/hosts.allow`, `/etc/hosts.deny`

外部からのアクセスの許可・拒否を設定する。
基本的な記述は「デーモン : ホスト」で、スペースかカンマで区切って複数指定することも可能。
ワイルドカード、ネットマスク、例外なども指定できる。
アクセスは以下の順序で評価され、該当した時点で評価終了。

1.  `/etc/hosts.allow` の記述に該当していれば許可:

        ALL : localhost
        sshd : 130.34.107.

2.  `/etc/hosts.deny` の記述に該当していれば拒否:

        ALL : ALL

3.  どちらにも該当しなかったアクセスはすべて許可
