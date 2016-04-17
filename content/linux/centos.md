+++
title = 'CentOS 6.5'
tags = ["linux"]
[menu.main]
  parent = "linux"
+++

Red Hat系

## 初期設定

### `sudo` 許可

1.  `root` になる:

        % su -

2.  `visudo` 実行して `wheel` グループについてのコメントアウト削除:

        %wheel ALL=(ALL) ALL

    {{%div class="note"%}}
`x` で1文字削除して `:wq` で保存終了
    {{%/div%}}

3.  ユーザーを `wheel` に加える:

        # usermod -G wheel watal
        # id watal

4.  `root` もユーザーも一旦ログアウト:

        # exit
        % exit

### `yum` 関連

1.  国内ミラーだけ使うように
    `/etc/yum/pluginconf.d/fastestmirror.conf` を編集:

        include_only=.jp

2.  EPEL (Extra Packages for Enerprise Linux) リポジトリを追加:

        % su -c 'rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm'

3.  [RepoForge](http://repoforge.org/use/) (旧RPMforge) リポジトリを追加:

        % wget http://pkgs.repoforge.org/rpmforge-release/rpmforge-release-0.5.3-1.el6.rf.x86_64.rpm
        % sudo rpm -Uvh rpmforge-release-0.5.3-1.el6.rf.x86_64.rpm

4.  とりあえずアップデート:

        % sudo yum update

### 開発環境

-   `gcc`, `make`, `git` などをまとめて入れる:

        % yum groupinfo "Development Tools"
        % sudo yum groupinstall "Development Tools"

-   ログインシェルを `zsh` に。 cf. [/dev/zsh]({{< relref "dev/zsh.md" >}})
-   最新の Python を入れる。cf. [/python/install]({{< relref "python/install.md" >}})
-   ソースから Mercurial を入れる。cf. [/dev/mercurial]({{< relref "dev/mercurial.md" >}})
-   EPELリポジトリを追加して `sudo yum install R`
-   `git` はRepoForgeでも古い1.7しか入らないのでソースから:

        % sudo yum install curl-devel expat-devel gettext-devel openssl-devel zlib-devel
        % wget -O- https://github.com/git/git/archive/v2.3.1.tar.gz | tar xz
        % cd git-2.3.1/
        % make prefix=${HOME}/local all
        % make prefix=${HOME}/local install

-   `emacs` も23だと Cask が使えないのでソースから？
