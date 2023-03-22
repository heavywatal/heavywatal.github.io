+++
title = 'ssh'
tags = ["shell", "communication"]
[menu.main]
  parent = "dev"
+++

<https://www.openssh.com/>

## 公開鍵によるパスワード無しの認証

秘密鍵 `~/.ssh/id_ed25519` を持つホストから、
公開鍵 `~/.ssh/authorized_keys` を持つホストにログインできるようにする。

手元のコンピュータ1台につき1ペア生成して、公開鍵をリモートホストに登録するのが基本。
リモートホストごとに違う鍵を使うことも可能。

1.  手元のコンピュータで鍵ペア (秘密鍵 `id_ed25519` と公開鍵 `id_ed25519.pub`) を作成:
    ```
    ssh-keygen -t ed25519
    # Generating public/private ed25519 key pair.
    # Enter file in which to save the key (~/.ssh/id_ed25519): <return>
    # Enter passphrase (empty for no passphrase):              パスフレーズ<return>
    # Enter same passphrase again:                             パスフレーズ<return>
    # Your identification has been saved in ~/.ssh/id_ed25519
    # Your public key has been saved in ~/.ssh/id_ed25519.pub
    ```
    - Ed25519を使えない古い環境ではRSA (`-t rsa -b 4096`)
      やECDSA (`-t ecdsa -b 521`) を使う。
    - 鍵ファイルの名前はデフォルトのまま<kbd>return</kbd>が簡単。
      特定のリモートサーバー用に別の鍵を用意する場合だけ変える。
      使用中のものを上書きしてしまわないように注意。
    - コメントは「この鍵がどこで作られたか」の人間用メモ。
      ホスト名がまともならデフォルト `-C username@hostname` の形でいい。
    - ここで入力するパスフレーズはログインパスワードではなく鍵の暗号化に使われる。
      パスフレーズ無しのほうが断然楽ちんなのでそうされることが多いが、
      ssh-agentやKeyChainを設定しておけば入力をほとんど省略できるので、
      簡単なものでも設定しておいたほうが鍵ファイル流出に強くなる。

1.  できあがった鍵ペアのうち公開鍵のほう `~/.ssh/id_ed25519.pub`
    をリモートホストの `~/.ssh/authorized_keys` に追加する。
    やり方は場合によって異なる。

    - ウェブブラウザからアップロードする方法。
      (e.g., [GitHub > Settings](https://github.com/settings/keys))。
    - ターミナルからコマンドで送り込む方法
      (パスワードなど別の認証でログインできる場合):
      ```sh
      ssh-copy-id -i ~/.ssh/id_ed25519.pub watal@example.com
      # or
      ssh watal@example.com "mkdir ~/.ssh; touch ~/.ssh/authorized_keys"
      cat id_ed25519.pub | ssh watal@example.com "cat >> ~/.ssh/authorized_keys"
      ```

1.  `~/.ssh/` ディレクトリ内のパーミッションを確認:
    ```sh
    ls -al ~/.ssh
    # drwx------  .
    # -rw-------  authorized_keys
    # -rw-------  config
    # -rw-------  id_ed25519
    # -rw-r--r--  id_ed25519.pub
    # -rw-------  known_hosts
    ```
    外に出していい公開鍵 `.pub` 以外はユーザー本人だけが読み書きできるように。
    必要があれば次のようなコマンドで修正:
    ```sh
    chmod 700 ~/.ssh
    chmod 600 ~/.ssh/id_ed25519
    chmod 600 ~/.ssh/authorized_keys
    ```


## 設定

### `~/.ssh/config`

`ssh` する側のユーザー毎の設定。
設定項目は `ssh_config` と同じ。
先に書いたものが優先されるので、一般設定は最後に:

```
Host beast
  Hostname 192.168.6.66
  User eddie
Host *.ddbj.nig.ac.jp
  User heavywatal
  RequestTTY yes
Host *
  PasswordAuthentication no
  KbdInteractiveAuthentication no
  GSSAPIAuthentication no
  StrictHostKeyChecking no
  VisualHostKey yes
  AddKeysToAgent yes
  UseKeychain yes
  IdentityFile ~/.ssh/id_ed25519
```

ユーザー名やアドレスをセットにしてニックネームをつけることで、入力を省略できる。
例えば上のような設定を `~/.ssh/config` に書いておけば以下の２つは等価:
```sh
ssh eddie@192.168.6.66
ssh beast
```

macOSで `UseKeychain yes` にしておけばパスフレーズ入力は最初の1回だけで済む。
パスフレーズ無しのような手軽さで秘密鍵の流出に強くなる。
`AddKeysToAgent yes` を設定せず
`ssh-add -l` が `The agent has no identities.` を返す状態でも行ける。


### `/etc/ssh_config`

[`man ssh_config`](https://manpages.ubuntu.com/manpages/man5/ssh_config.5.html)

`ssh` する側のコンピュータ全体の設定。
普通は `~/.ssh/config` で事足りるのでいじらない。
サーバー間の通信を調整する管理者向け。

### `/etc/sshd_config`

[`man sshd_config`](https://manpages.ubuntu.com/manpages/man5/sshd_config.5.html)

`ssh` される側のデーモンの設定。余計な入り口を塞ぐべし:
```
PermitRootLogin no
PasswordAuthentication no
KbdInteractiveAuthentication no
KerberosAuthentication no
GSSAPIAuthentication no
UsePAM no
```

`PermitRootLogin`
: `without-password` は `prohibit-password` に取って代わられた。

`KbdInteractiveAuthentication`
: `ChallengeResponseAuthentication` は取って代わられた。


## ファイル転送

ファイルひとつなら`scp`でもいいけどそれ以上なら
[rsync]({{< relref "rsync.md" >}}) を使ったほうがよい。
あるいは[sshfsでマウント]({{< relref "mount.md#sshfs" >}})してしまうのも楽ちん。


## 環境変数

```
echo $SSH_CONNECTION
[client IP] [client port] [server IP] [server port]
```

どこから `ssh` したか、を取得するには
`echo $SSH_CONNECTION | awk '{print $1}'`
