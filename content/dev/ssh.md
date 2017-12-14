+++
title = 'ssh'
tags = ["shell", "communication"]
[menu.main]
  parent = "dev"
+++

<http://www.openssh.com/>

## 公開鍵によるパスワード無しの認証

1.  パスフレーズ無しで鍵ペア (秘密鍵 `id_rsa` と公開鍵 `id_rsa.pub`) を作る:

        % ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
        Generating public/private rsa key pair.
        Enter file in which to save the key ($HOME/.ssh/id_rsa): ＃そのままリターン
        Created directory '$HOME/.ssh'.
        Enter passphrase (empty for no passphrase): ＃そのままリターン
        Enter same passphrase again: ＃そのままリターン
        Your identification has been saved in $HOME/.ssh/id_rsa.
        Your public key has been saved in $HOME/.ssh/id_rsa.pub.

2.  できあがった公開鍵を `authorized_keys` にリネームまたは追加。
    公開されるのはこの `authorized_keys`:

        % cd ~/.ssh/
        % mv id_rsa.pub authorized_keys

    既に存在する `authorized_keys` に追記する場合は↓:

        % cat id_rsa.pub >> authorized_keys
        % rm id_rsa.pub

3.  パーミッションを設定:

        % chmod 700 ~/.ssh
        % chmod 600 ~/.ssh/id_rsa
        % chmod 600 ~/.ssh/authorized_keys

4.  リモートホスト・ローカルホストの `~/.ssh/` に鍵を適宜配置する。
    秘密鍵を持っているホストから公開鍵を持っているホストにパスワード無しで
    `ssh`, `scp` できるようになる。

{{%div class="note"%}}
古い環境を考慮しなくていい場合はRSAより高性能なEd25519を使い、
鍵ファイルの名前も適宜読み替える。
{{%/div%}}


## Configuration

http://man.openbsd.org/OpenBSD-current/man5/ssh_config.5

### `~/.ssh/config`

`ssh` する側のユーザー毎の設定。
ユーザー名やアドレスをセットにしてニックネームをつけることで、入力を省略できる。
先に書いたものが優先されるので、一般設定は最後に:

    Host finlandia
      Hostname 192.168.0.26
    Host tapiola
      Hostname 192.168.0.112
    Host beast
      Hostname 192.168.6.66
      User ironmaiden
    Host *
      User sibelius
      Protocol 2
      GSSAPIAuthentication no
      StrictHostKeyChecking no

上の設定を `~/.ssh/config` に書いておけば以下の２つは等価:

    % ssh -2 sibelius@192.168.0.112
    % ssh tapiola

「ローカル - 中継機 - 計算機」でローカルから計算機へ一気に多段 `ssh`
をキメる(`ssh keisanki1`)には以下のように書いておく。
でも一旦中継機にログインしてから計算機に何かしたいとき（少なくない）には余計な感じになってしまう。
「中継機にいないときは」みたいな条件文を書ければ最高なんだけど:

    Host keisanki1 keisanki2
      ProxyCommand nohup ssh chuukeiki nc %h %p

### `/etc/ssh_config`:

`ssh` する側のマシン全体の設定。認証の試行を絞ることで少しでも高速にしたい:

    Protocol 2
    ChallengeResponseAuthentication no
    GSSAPIAuthentication no

### `/etc/sshd_config`:

`ssh` される側のデーモンの設定。余計な入り口を塞ぐべし:

    Protocol 2
    PermitRootLogin without-password
    ChallengeResponseAuthentication no
    PasswordAuthentication no
    KerberosAuthentication no
    GSSAPIAuthentication no
    UsePAM no

## Environmental variables

    % echo $SSH_CONNECTION
    [client IP] [client port] [server IP] [server port]

どこから `ssh` したか、を取得するには:

    % echo $SSH_CONNECTION | awk '{print $1}'

## ファイル転送

ファイルひとつなら`scp`でもいいけどそれ以上なら
[rsync]({{< relref "rsync.md" >}})を使ったほうがよい。
