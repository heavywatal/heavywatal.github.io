+++
title = 'apt/dpkg'
tags = ["linux", "package"]
[menu.main]
  parent = "linux"
+++

-   Applications &gt; Ubuntu Software Center
-   System &gt; Administration &gt; Synaptic Package Manager, Update Manager
-   Applications &gt; Accessories &gt; Terminal で `apt-get`, `apt-cache`, `dpkg`

## apt-get

<https://help.ubuntu.com/community/AptGet/Howto>

大概のことは [Synaptic](https://help.ubuntu.com/community/SynapticHowto) で済むけど、コマンドのほうが早い場合もある。

-   パッケージリストをアップデートし、アップグレード可能なものがあるか確認。したらアップグレード:

        % sudo apt-get update && apt-get -s upgrade
        % sudo apt-get -u upgrade

-   不要なものを自動判別して削除:

        % sudo apt-get autoclean
        % sudo apt-get autoremove

-   パッケージの設定をインストール時の状態からやり直す:

        % dpkg-reconfigure [options] packages

-   インストール済みのパッケージ一覧。適当に `grep` して使うべし:

        % dpkg -l linux-image* | grep ii
        ii  linux-image-2.6.32-22-generic         2.6.32-22.36                                        Linux kernel image for version 2.6.32 on x86
        ii  linux-image-2.6.32-23-generic         2.6.32-23.37                                        Linux kernel image for version 2.6.32 on x86
        ii  linux-image-generic                   2.6.32.23.24                                        Generic Linux kernel image

## Repository

[Medibuntu](http://medibuntu.org/):

    sudo wget --output-document=/etc/apt/sources.list.d/medibuntu.list http://www.medibuntu.org/sources.list.d/$(lsb_release -cs).list
    sudo apt-get --quiet update
    sudo apt-get --yes --quiet --allow-unauthenticated install medibuntu-keyring
    sudo apt-get --quiet update

[Dropbox](http://dropbox.com/):

    sudo apt-key adv --keyserver pgp.mit.edu --recv-keys 5044912E
    sudo add-apt-repository "deb http://linux.dropbox.com/ubuntu $(lsb_release -sc) main"

[Personal Package Archives for Ubuntu](https://launchpad.net/ubuntu/+ppas)
公式版のパッケージは安定性のために新しさを犠牲にしている場合が多い。か
と言って、ソースをダウンロードして自分でビルド・インストールするのは難しい。
誰かが非公式にビルドした最新版がこのPPAで提供されていれば、
そのリポジトリを追加することで簡単にインストールできる。
以下のように、その提供者の\`ppa:&lt;user-name&gt;/&lt;ppa-name&gt;さえ分かれば
コマンド１行でリポジトリと鍵を登録できる。
このコマンドは python-software-properties というパッケージで提供されている。

::

    sudo add-apt-repository ppa:&lt;user-name&gt;/&lt;ppa-name&gt;
    sudo apt-get update

Toolchain test builds &lt;<https://launchpad.net/~ubuntu-toolchain-r/+archive/test>&gt;\_: ppa:ubuntu-toolchain-r/test\`
GCCなど開発ツールの最新版を自分でビルドしたくない場合はこれを。
ただし `base-files` パッケージのインストールにより
`/etc/issue` や `/etc/lsb-release` が書き換わったりして
困ったことになる場合があるのでそうなったら手で書きなおす。
[ref.](http://askubuntu.com/questions/126498/ubuntu-12-04-reports-itself-as-quantal)

[PPA for Ubuntu Japanese Team](https://launchpad.net/~japaneseteam/+archive/ppa/): `ppa:japaneseteam/ppa`
:   日本語環境のサポート

[PPA for Japanese packages for testers](https://launchpad.net/~japanese-testers/+archive/ppa): `ppa:japanese-testers/ppa`
:   日本語入力 (ibus) らへんのアップデート
