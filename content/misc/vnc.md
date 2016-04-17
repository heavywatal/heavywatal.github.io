+++
title = 'VNCによる画面共有'
[menu.main]
  parent = "misc"
+++

Virtual Network Computing

## サーバー側

### Ubuntu 12.04 LTS Precise Pangolin

1.  vnc4server をインストール:

        % sudo apt-get install vnc4server

2.  一旦起動してパスワードを設定:

        % vncserver
        Password:
        Verify:
        % vncserver -kill :1

3.  設定ファイル `$HOME/.vnc/xstartup` を書き換える:

        #!/bin/sh

        # Uncomment the following two lines for normal desktop:
        unset SESSION_MANAGER
        # exec /etc/X11/xinit/xinitrc
        gnome-session --session=ubuntu-2d &

        [ -x /etc/vnc/xstartup ] && exec /etc/vnc/xstartup
        [ -r $HOME/.Xresources ] && xrdb $HOME/.Xresources
        xsetroot -solid grey
        vncconfig -iconic &
        #x-terminal-emulator -geometry 80x24+10+10 -ls -title "$VNCDESKTOP Desktop" &
        #x-window-manager &

4.  オプションを与えつつ起動 (省略時は `:1 -geometry 1024x768`):

        % vncserver :6 -geometry 1536x1024

    コロンの後の数字はディスプレイ番号。
    アクセスするときのポートが **5900 + ディスプレイ番号** になる。

5.  リソース節約のため、使い終わったらディスプレイ番号を指定して止める:

        % vncserver -kill :6

### Mac

1.  `System Preferences... --> Sharing --> Screen Sharing`
    をオンにする。
2.  `Computer Settings...` で両方にチェックし、パスワード設定

## クライアント側

### Linux

1.  Remote Desktop Client を起動
2.  Server に `{SERVER}:{PORT}` を指定
3.  ユーザー名やパスワードを適宜入力して接続

### Mac

1.  Finder をアクティブにして `command + k`
    あるいは `Go --> Connect to Server...`
2.  Server Address: に以下のように入力して接続:

        vnc://charles:5901

### SSH port forwarding

ただのVNCだと暗号化されていない情報がダダ漏れなので、
SSHに乗せてセキュアな通信を行うべきである。
公開鍵によるパスワード無し認証やユーザ名の省略など、
事前に [SSHの設定]({{< relref "dev/ssh.md" >}}) をしておくとよい。

1.  転送するポートを指定しつつリモートホストにSSH接続:

        ssh -L {port}:{host}:{hostport} {host}

    クライアント側のポートは何番でもよい。例えば:

        ssh -L 9999:charles:5906 charles

2.  リモートでVNCサーバを起動:

        % vncserver :6 -geometry 1536x1024

3.  ローカルホスト宛てにVNC接続するとSSHが転送してくれる:

        vnc://localhost:9999
