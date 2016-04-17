+++
title = 'MySQL'
tags = ["database"]
[menu.main]
  parent = "dev"
+++

[<http://dev.mysql.com/>](http://dev.mysql.com/)

The MySQL™ software delivers a very fast, multi-threaded, multi-user, and robust SQL (Structured Query Language) database server.

## Operation

ログイン:

    % mysql -u root -p
    Enter password:

文末にはセミコロン:

    mysql> show databases;
    mysql> show status;

`like` とパーセント記号によるワイルドカードで絞り込める:

    mysql> show variables like "%log%"

## Administration

ユーザーの指定方法。ホストのアドレスにはワイルドカードが使える。パスワード・ユーザー名・ホストは省略可能:

    USERNAME@'HOST' identified by 'PASSWORD'
    meme@'130.34.107.135'
    klabo@'130.34.107.%'
    biol@'%.biology.tohoku.ac.jp'
    root@''
    ''@localhost

ユーザーの一覧:

    mysql> select User,Host,Password from mysql.user;

ユーザーの追加・削除:

    mysql> create user [USER];
    mysql> drop user [USER];

パスワードの設定・変更:

    mysql> set password for USERNAME@HOST=password("PASSWORD");

権限を与える:

    mysql> show grants for [USER];
    mysql> grant select on *.* to [USER];

設定を反映:

    mysql> flush privileges;

ログを削除:

    mysql> show binary logs;
    mysql> purge master logs before now();
