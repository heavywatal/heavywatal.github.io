+++
title = 'Mailman'
[menu.main]
  parent = "misc"
+++

Webベースのメーリングリスト管理システム

公式サイト: <http://www.gnu.org/software/mailman/>\
日本語情報: <http://docs.python.jp/contrib/mailman/>

## Install on Mac OS 10.8 Mountain Lion Server

新しい Mac サーバーでは標準搭載されていないらしいので、
ソースコードから自分でセットアップする必要がある。

やり方は基本的に
<http://www.livetime.com/mountain-lion-mailman-mailing-list/>
を参考にした。

1.  App Store から Xcode をインストール
2.  Xcode を起動して、設定画面から
    Command Line Tools をインストール
3.  ここからは Terminal.app 上での作業となる。
    インストール先のディレクトリを準備。
    グループはアンダースコアが付いて `_mailman` となるが気にしない:

        % sudo mkdir /usr/local/mailman
        % cd /usr/local/mailman
        % sudo chgrp mailman .
        % sudo chmod a+rx,g+ws .

4.  Mailman をダウンロードして展開:

        % cd /tmp/
        % wget -O- http://ftp.gnu.org/gnu/mailman/mailman-2.1.15.tgz | tar xz
        % cd mailman-2.1.15/

5.  `/usr/bin/python` が使われるように注意して
    `configure` & `make`。
    初めから `sudo` してやるほうがいいかも:

        % sudo ./configure
        % sudo make
        % sudo make install

6.  インストール先に行ってアクセス権を修復:

        % cd /usr/local/mailman/
        % sudo ./bin/check_perms -f

7.  `/etc/apache2/extra/httpd_mailman.conf` を作る

    ```apache
    ScriptAlias /mailman/ "/usr/local/mailman/cgi-bin/"
    Alias /pipermail/ "/usr/local/mailman/archives/public/"
    Alias /icons/ "/usr/local/mailman/icons/"
    <Directory "/usr/local/mailman/archives/public/">
        Options FollowSymLinks MultiViews Indexes
        AllowOverride None
        Order allow,deny
        Allow from all
    </Directory>
    ```

8.  `/Library/Server/Web/Config/apache2/httpd_server_app.conf` を編集する。
    550行目らへん、SSL/TLSのコメントのすぐ上あたりに下記の行を追加:

        Include /private/etc/apache2/extra/httpd-mailman.conf

9.  Apache ウェブサーバを再起動:

        % sudo apachectl graceful

10. ウェブブラウザから `http:://{SERVER_ADDRESS}/mailman/admin` にアクセスしてみる。
    エラーが発生していたら
    `/usr/local/mailman/logs/` 以下あるいは
    `/var/log/apache2/` 以下のログを見てみる。
11. `/usr/local/mailman/Mailman/mm_cfg.py` に以下を追記:

        MTA = `Postfix`

12. エイリアスを作ってアクセス権を調整:

        % sudo ./bin/genaliases
        % sudo chown mailman:mailman data/aliases*
        % sudo chmod g+w data/aliases*

13. `/Library/Server/Mail/Config/postfix/main.cf` に以下を追記。
    当方の環境では `hash:/etc/aliases` を書くと失敗した:

        alias_maps = hash:/usr/local/mailman/data/aliases

14. Postfix メールサーバの設定を再読み込み:

        % sudo postfix reload

15. Mailman のサイトパスワードを設定:

        % sudo /usr/local/mailman/bin/mmsitepass XXXXXXXX

16. 新しいメーリスを作ってみる(必要かどうかは不明):

        % sudo ./bin/newlist mailman
        % sudo ./bin/config_list -i data/sitelist.cfg mailman

17. qrunner を起動:

        % sudo ./bin/mailmanctl start

18. `/System/Library/LaunchDaemons/org.list.mailmanctl.plist`
    を作り、システム起動後に自動的に始まるようにする

    ```xml
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
    <plist version="1.0">
    <dict>
        <key>Debug</key>
        <true/>
        <key>Disabled</key>
        <false/>
        <key>Label</key>
        <string>org.list.mailmanctl</string>
        <key>Program</key>
        <string>/usr/local/mailman/bin/mailmanctl</string>
        <key>ProgramArguments</key>
        <array>
            <string>mailmanctl</string>
            <string>-s</string>
            <string>start</string>
        </array>
        <key>KeepAlive</key>
        <false/>
        <key>RunAtLoad</key>
        <true/>
        <key>AbandonProcessGroup</key>
        <true/>
    </dict>
    </plist>
    ```

19. リストを引き継ぐには `mailman/` ディレクトリの
    `lists/` と `archives/` をコピーすれば良い。
    古いサーバーから移行してきた場合、それらは
    `/Library/Server/Migrated/var/mailman` にある:

        % sudo rsync -auv /Library/Server/Migrated/var/mailman/lists/ /usr/local/mailman/lists/
        % sudo rsync -auv /Library/Server/Migrated/var/mailman/archives/ /usr/local/mailman/archives/

20. エイリアスを作り直し、アクセス権を修復し、再起動:

        % sudo ./bin/genaliases
        % sudo ./bin/check_perms -f
        % sudo ./bin/mailmanctl restart

21. テストメールを送って確認。
