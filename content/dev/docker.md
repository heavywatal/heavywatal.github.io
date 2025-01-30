+++
title = "Docker"
tags = ["package"]

[menu.main]
  parent = "dev"
+++

<https://docs.docker.com/>

ホストOSから隔離された環境でアプリケーションを動かすためのプラットフォーム。
ゲストOSを丸ごと動かす仮想マシンとは異なり、ホストLinuxのカーネルを使うので軽い。
LinuxじゃないホストOSでは仮想化が挟まるぶんのコストが結局かかる。

管理者権限を持っている状況でサービスを開発・運用するのに向いている。
HPCにおけるCLIツール導入や一般ユーザーからの利用が目的なら [Apptainer]({{< relref "apptainer.md" >}}) のほうが便利。

## Get Started

### Glossary

container
: ホストOSから隔離されたプロセスやファイルシステムを持つ実行環境。
  目的のアプリケーションを動かすのに必要な構成要素のうちkernelより上をすべて含む。
  imageを読み込んで起動したインスタンス。

image
: containerの構成要素を集めて固めたファイル。
: Docker独自の形式をベースに[OCI](https://github.com/opencontainers/image-spec)として標準化されている。

dockerfile
: imageのレシピ。


### Linux

<https://docs.docker.com/engine/install/ubuntu/>

apt repository を追加してインストール。

`sudo` なしで実行できるようにするには:
```sh
sudo gpasswd --add $(whoami) docker
sudo chgrp docker /var/run/docker.sock
sudo systemctl restart docker
exit
# login again
docker --version
```


### Mac

Docker Desktop より [Orbstack](https://orbstack.dev/) のほうが軽くて使いやすい。
いずれにせよ [Homebrew]({{< relref "homebrew.md" >}}) で入れるのが簡単:
```sh
brew install --cask orbstack
open -a OrbStack
```

OrbStack.app を起動するとシェルの設定を勝手にいじってパスを通してくれる:
```sh
orb --help
docker --version
docker-compose --version
```

Silicon Mac で x86_64/amd64 のイメージを動かすには
`--platform linux/amd64`
のようなオプションを明示的に与えてRosettaを介す。
<https://docs.orbstack.dev/docker/#intel-x86-emulation>


### Hello world

```sh
docker image pull hello-world
docker image ls
docker container run hello-world
docker container ls --all
```

`run` で手元にimageが見つからなければ勝手に `pull` される。
containerは実行終了後も残る。


## CLI

操作対象がcontainerかimageかなどによってサブコマンドを使う。
前は `docker pull`, `docker run` のようにフラットなコマンド体系だった。

### `docker container` subcommands

<https://docs.docker.com/engine/reference/commandline/container/>

`docker container ls`
:   実行中のcontainerを表示。
:   `-a`: 終了後のものも表示。
:   `-q`: IDのみ表示。

`docker container rm CONTAINER`
:   containerを削除。
:   停止中のものを一括削除するなら `docker container prune`

`docker container run [OPTIONS] IMAGE [COMMAND]`
:   新しいcontainerを走らせてコマンドを実行。
:   `-i, --interactive`: 標準入力を開けておく。
:   `-t, --tty`: 標準出力を開けておく。
:   `-d, --detach`: バックグラウンドで実行させておく。
    うっかりフォアグラウンドで起動してしまっても
    <kbd>^p</kbd><kbd>^q</kbd> でデタッチできる。
:   `--rm`: 終了後に削除。これを指定しないと `ls -a` のリストに残る。
:   `--name string`: コンテナに名前を付ける。 `exec` で指定するときに便利。
:   `--mount`: containerから外のファイルシステムにアクセスできるように割り当てる。
    - `type=volume,src=<VOLUME-NAME>,dst=<CONTAINER-PATH>`:
      Docker管理下に場所を確保する。
      container再起動やcontainer間共有のための永続化にはこっちを使う。
      Docker外でも使うデータを読み書きするのには向かない。
    - `type=bind,src=<HOST-PATH>,dst=<CONTAINER-PATH>`: ホスト側の場所を直接指定する。
      古い `-v` オプションによる指定は非推奨。
    -	`--mount type=bind,src="$PWD",dst="$PWD"`\
      `--workdir "$PWD"`\
      `--user "$(id -u):$(id -g)"`\
      のようにするとcontainer内のプログラムにカレント以下を渡せる。
      Apptainerならこのへんの設定をデフォルトでやってくれる。

`docker container exec [OPTIONS] IMAGE COMMAND`
:   既に走ってるcontainerでコマンドを実行。

`docker container stats`


### `docker image` subcommands

<https://docs.docker.com/engine/reference/commandline/image/>

`docker image ls`
:   手元にあるimageを一覧表示。

`docker image pull NAME:TAG`
:   Docker Hub などからimageをダウンロード。


### misc.

```sh
docker system info
docker system df
```

## Registry

[BioContainers](https://biocontainers.pro/registry/)
: [BioConda](https://bioconda.github.io/) recipes を使っているらしい。
: 実際のregistry機能をホストしているのは他所のサーバー:
  - <https://hub.docker.com/u/biocontainers> OCI
  - <https://quay.io/organization/biocontainers> OCI
  - <https://depot.galaxyproject.org/singularity/> SIF.
    [singularity-build-bot](https://github.com/BioContainers/singularity-build-bot)
    がQuayからOCIを取得して変換。
    通信が遅く、~500MB以上のイメージをpullしようとするとタイムアウトになって使えない。
    一旦 `wget` か何かで落とすか、自分でOCIから変換するほうが早くて確実。

[Docker Hub](https://hub.docker.com/)
: 100 pulls / 6 hours の制限がある。
  無料のPersonalアカウントで認証していれば200。

[Quay.io](https://quay.io/)
: 無料アカウントは無い。pull制限は?

[Google Container Registry](https://gcr.io)
: シェルさえ含まず軽量セキュアな [distroless](https://github.com/GoogleContainerTools/distroless) を提供。

[GitHub Container Registry (GitHub Packages)](https://ghcr.io)

[GitLab Container Registry](https://docs.gitlab.com/ee/user/packages/container_registry/)
