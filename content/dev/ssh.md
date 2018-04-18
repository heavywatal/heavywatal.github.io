+++
title = 'ssh'
tags = ["shell", "communication"]
[menu.main]
  parent = "dev"
+++

<http://www.openssh.com/>

## 公開鍵によるパスワード無しの認証

秘密鍵 `~/.ssh/id_ed25519` を持つホストから、
公開鍵 `~/.ssh/authorized_keys` を持つホストに、
パスワード無しでログインできるようにする。

1.  ローカルホストで鍵ペア (秘密鍵 `id_ed25519` と公開鍵 `id_ed25519.pub`) を作成:
    ```
    % ssh-keygen -t ed25519 -N ''
    Generating public/private ed25519 key pair.
    Enter file in which to save the key (~/.ssh/id_ed25519): ＃そのままリターン
    Your identification has been saved in ~/.ssh/id_ed25519.
    Your public key has been saved in ~/.ssh/id_ed25519.pub.
    ```
    Ed25519を使えない古い環境ではRSA (`-t rsa -b 4096`)
    やECDSA (`-t ecdsa -b 521`) を使う。

1.  できあがった公開鍵をリモートホストの `~/.ssh/authorized_keys` に追加:
    ```sh
    % ssh-copy-id -i ~/.ssh/id_ed25519.pub watal@example.com
      # or
    % ssh watal@example.com "mkdir ~/.ssh; touch ~/.ssh/authorized_keys"
    % cat id_ed25519.pub | ssh watal@example.com "cat >> ~/.ssh/authorized_keys"
    ```

1.  ユーザー本人だけが読み書きできるパーミッションに設定されていることを確認:
    ```sh
    % ls -l ~/.ssh
    % ssh watal@example.com "ls -l ~/.ssh"
      # if necessary
    % chmod 700 ~/.ssh
    % chmod 600 ~/.ssh/id_ed25519
    % chmod 600 ~/.ssh/authorized_keys
    ```


## 設定

http://man.openbsd.org/OpenBSD-current/man5/ssh_config.5

### `~/.ssh/config`

`ssh` する側のユーザー毎の設定。
ユーザー名やアドレスをセットにしてニックネームをつけることで、入力を省略できる。
先に書いたものが優先されるので、一般設定は最後に:

```
Host beast
  Hostname 192.168.6.66
  User eddie
Host *.ddbj.nig.ac.jp
  User heavywatal
  RequestTTY yes
Host *
  Protocol 2
  GSSAPIAuthentication no
  StrictHostKeyChecking no
  VisualHostKey yes
```

例えば上のような設定を `~/.ssh/config` に書いておけば以下の２つは等価:
```
% ssh -2 eddie@192.168.6.66
% ssh beast
```

### `/etc/ssh_config`:

`ssh` する側のマシン全体の設定。認証の試行を絞ることで少しでも高速にしたい:

```
Protocol 2
ChallengeResponseAuthentication no
GSSAPIAuthentication no
```

### `/etc/sshd_config`:

`ssh` される側のデーモンの設定。余計な入り口を塞ぐべし:
```
Protocol 2
PermitRootLogin without-password
ChallengeResponseAuthentication no
PasswordAuthentication no
KerberosAuthentication no
GSSAPIAuthentication no
UsePAM no
```


## ファイル転送

ファイルひとつなら`scp`でもいいけどそれ以上なら
[rsync]({{< relref "rsync.md" >}}) を使ったほうがよい。
あるいは[sshfsでマウント]({{< relref "mount.md#sshfs" >}})してしまうのも楽ちん。


## 環境変数

```
% echo $SSH_CONNECTION
[client IP] [client port] [server IP] [server port]
```

どこから `ssh` したか、を取得するには
`echo $SSH_CONNECTION | awk '{print $1}'`
