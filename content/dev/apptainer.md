+++
title = "Apptainer"
tags = ["package"]

[menu.main]
  parent = "dev"
+++

<https://apptainer.org/>

- [Docker]({{< relref "docker.md" >}}) のようなもの。
- 管理者権限なしで実行できるので一般ユーザーに使わせやすい。
- ホスト機のファイル読み書きもデフォルトでやりやすい。
- イメージはOCI形式ではなく[SIF](https://github.com/apptainer/sif)形式。
  拡張子は `.sif` が標準的で `.simg` も無印も見かける。

Singularityだったものが[Linux Foundationへの移管](https://apptainer.org/news/community-announcement-20211130/)に伴って改名。
紛らわしいことに[Sylabs社がSingularity CEと称しているfork](https://github.com/sylabs/singularity)はとりあえず無視。


## Admin

<https://apptainer.org/docs/admin/latest/>


### Installation

<https://apptainer.org/docs/admin/latest/installation.html>

#### Linux

```sh
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

setuid版というのはかなり古いkernel向けっぽいのでとりあえず無視。


#### Mac

kernelが違うのでネイティブには動かない。
今のところDocker Desktop的なものも無いので自分で仮想環境を用意する必要があるらしい。


### Configuration

<https://apptainer.org/docs/admin/latest/configfiles.html>

```sh
apptainer help
apptainer version
apptainer buildcfg
apptainer remote list
```

`/etc/singularity` や `/usr/local/etc/singularity` が残ってると怒られるので消す。
`singularity` コマンドは互換性のためのエイリアスとしてしばらく提供されるっぽいけど。


## User

<https://apptainer.org/docs/user/latest/>

### CLI

<https://apptainer.org/docs/user/latest/cli.html>

[`pull`](https://apptainer.org/docs/user/latest/cli/apptainer_pull.html)
: `apptainer pull [pull options...] [output file] <URI>`
: SIFイメージをダウンロードする。
  OCIイメージだったら `build` を呼んでSIFに変換する。
  e.g., `apptainer pull docker://alpine` すると `alpine_latest.sif` ができる。

[`exec`](https://apptainer.org/docs/user/latest/cli/apptainer_exec.html)
: `apptainer exec [exec options...] <container> <command>`
: コマンドを実行する。
: container引数にはSIFファイルへのパスを渡せる。
: `docker exec` との違い:
  - containerが先に走っている必要はない。 `--rm` も不要。
  - `--mount` 無しでも `/home/$USER`, `/tmp`, `$PWD` がbindされる。
  - `-it` 無しでも自然な入出力。

[`run`](https://apptainer.org/docs/user/latest/cli/apptainer_run.html)
: container内のrunscriptを実行する。
  runscriptは `/apptainer` (`/singularity`) に置かれたシェルスクリプトで、
  `exec` と同様にコマンドを実行できるようになっていることが多そう。

[`shell`](https://apptainer.org/docs/user/latest/cli/apptainer_shell.html)
: container内でシェルを起動する。
  `exec` や `run` で `bash` するのとどう違うか？


## Repositories

`docker://`
: [Docker向けregistry]({{< relref "docker.md#registry" >}})からpullしてSIFに変換。
  まともなSIF registryが見当たらない現状では結局これが楽ちん。

`shub://` [Singularity Hub](https://singularityhub.github.io/singularityhub-docs/)
: 2021-04 に凍結。
  [DataLad](https://datasets.datalad.org/?dir=/shub/)
  に引き継がれた既存イメージはpullできるらしい。

`library://` [Sylabs Singularity library](https://cloud.sylabs.io/library)
: Apptainerデフォルトではオフ。
