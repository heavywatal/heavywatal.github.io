+++
title = 'Pythonインストール'
toc = true
tags = ["python"]
[menu.main]
  parent = "python"
  weight = -99
+++

## 標準・公式

MacやLinuxならシステムの一部として
`/usr/bin/python3` が既にインストールされているので、
何もしなくても使えっちゃ使える。
でも大概そこに入ってるのはちょっと古いバージョンなので、
[python.jp/環境構築ガイド](https://www.python.jp/install/install.html)
に従って最新版を入れるのがよい。


## uv

<https://docs.astral.sh/uv/>

プロジェクトの環境構築を一切合切面倒見てくれる管理ツール。
ruffや[rye](#rye)と同じチームが開発していて、同じくrust製。
次のワンライナーでプログラム本体が `~/.cargo/bin/uv` に配置される:
```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

rust/cargoを使わない人は `source ~/.cargo/bin/env`
みたいな文をシェルの設定ファイルに追加してPATHを通す必要があるかも。

任意のバージョンのPythonを入れるためのツールとして、
つまり[pyenv](#pyenv)的な位置付けでも使える。
しかもビルド済みのPythonを
<https://github.com/indygreg/python-build-standalone>
から取ってくるので、自前ビルド環境に左右されずCPUも使わず簡単。

### Pythonインストーラーとして使う

<https://docs.astral.sh/uv/reference/cli/#uv-python>

```sh
# バージョン一覧
uv python list

# インストール
uv python install
```

バージョンを省略すると適当な最新版。
`3` とか `3.12` みたいな指定でもその中での最新版を入れられる。

`~/.local/share/uv/python/` 以下に配置される。
`uv run` や `uv venv` 越しに使う前提ならここにPATHを通す必要はない。

### 設定

<https://docs.astral.sh/uv/configuration/environment/>

設定ファイルは `~/.config/uv/` に置くらしいけど、
今のところ環境変数を使うのが主流。

ほかのツールで入れたPythonやシステム標準の `/usr/bin/python`
まで探しに行って報告しようとしてくれる。
uv自身で入れたPythonだけに専念してもらうと少し早くなる:
```sh
export UV_PYTHON_PREFERENCE=only-managed
```

[PEP 668](https://peps.python.org/pep-0668/) `EXTERNALLY-MANAGED` が有効なので
[`uv venv`](https://docs.astral.sh/uv/reference/cli/#uv-venv)
で仮想環境を作って使う。

[`uv pip`](https://docs.astral.sh/uv/reference/cli/#uv-pip)
は普通の[pip]({{< relref "pip.md" >}})と比べて圧倒的に速い。

ほかにも[サブコマンド](https://docs.astral.sh/uv/reference/cli/#uv)がたくさん。


## rye

<https://rye.astral.sh/>

プロジェクトの環境構築を一切合切面倒見てくれる管理ツール。
ruffや[uv](#uv)と同じチームが開発していて、同じくrust製。
次のワンライナーでプログラム本体や設定ファイルなどが `~/.rye/` に配置される:
```sh
curl -sSf https://rye.astral.sh/get | bash
```

内部的にuvを使っていて、役割がややかぶっている。
uvの開発が急速に進む一方、こちらはやや遅れてきて存在意義が揺らいできている...?

### Pythonインストーラーとしてだけ使う

パッケージ運用はまだ普通にpipとかでいいかなと思う人の使い方。

```sh
# インストール済みバージョン一覧
rye toolchain list

# 利用可能バージョン一覧
rye toolchain list --include-downloadable

# インストール
rye toolchain fetch cpython@3.12
```

今のところグローバルに `pip3 install` しても怒られない。
`global-python = true`
の設定でshimsにPATHが通っていればryeの管理下にあるPythonをプロジェクト外でも使えるが、逆に
`pyproject.toml` が存在するディレクトリでのみそれができないという問題
[#1121](https://github.com/astral-sh/rye/issues/1121)
がある。
shimsに頼らず自分でPATHを設定するworkaround:
```sh
if [ -d "${RYE_HOME:=${HOME}/.rye}" ]; then
  py_versions=($(ls "${RYE_HOME}/py" | sort -V))
  export PY_PREFIX=${RYE_HOME}/py/${py_versions[@]: -1}
  PATH=${PY_PREFIX}/bin:${PATH}:${RYE_HOME}/shims
  unset py_versions
fi
```


## pyenv

<https://github.com/pyenv/pyenv>

管理者権限なしでホーム以下にインストールできる。
ソースコードを取ってきて自前ビルドするという点で上記[uv](#uv)や[rye](#rye)と異なる。
共有ライブラリやフレームワークなどを有効にしたい場合に使いやすい。

1.  [Homebrew]({{< relref "homebrew.md" >}}) か
    [Git]({{< relref "git.md" >}}) を使ってpyenvをインストール:

    ```sh
    brew install pyenv
    # or
    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    mkdir -p ~/.pyenv/cache
    ```

1.  <https://github.com/pyenv/pyenv/wiki>
    を参考に依存ライブラリをインストールしておくとビルドが少し軽くなる。

1.  Pythonのインストール先を決める環境変数
    `PYENV_ROOT` を公式推奨の `~/.pyenv` に設定し、
    ついでにPATHも追加しておく。
    シェルの設定ファイル (e.g., `~/.bashrc`) に次のように追記:

    ```sh
    if [ -d "${PYENV_ROOT:=${HOME}/.pyenv}" ]; then
      py_versions=($(ls "${PYENV_ROOT}/versions" | sort -V))
      export PY_PREFIX=${PYENV_ROOT}/versions/${py_versions[@]: -1}
      PATH=${PY_PREFIX}/bin:$PATH
      unset py_versions
    fi
    ```

    `pyenv shell` や `pyenv local`
    を使ってPythonのバージョンを頻繁に切り替える場合は、
    [公式の説明](https://github.com/pyenv/pyenv#installation)どおりに
    `eval "$(pyenv init --path)"` や
    `eval "$(pyenv init -)"`
    を設定してshimsを使う方法のほうがいいかもしれないけど、
    そうでなければこのようにPATHだけ設定するほうが単純で、
    起動時間も短くなる。

1.  シェルを再起動して設定を反映し、
    目当てのバージョンを探してインストール:

    ```sh
    exec $SHELL -l
    pyenv install -l | less
    pyenv install 3.12
    exec $SHELL -l
    ```

    R から [`reticulate`](https://rstudio.github.io/reticulate/)
    越しに呼ぶ場合は共有ライブラリを有効にしてビルドする:
    ```sh
    env PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.10
    # or
    Rscript -e 'reticulate::install_python("3.10")'
    ```

1.  [pip]({{< relref "pip.md" >}}) のパスを確認し、パッケージを入れる:

    ```sh
    which pip3
    pip3 install -U setuptools pip
    pip3 install -r /path/to/requirements.txt
    ```

    よく使うパッケージは
    [`requirements.txt`](https://github.com/heavywatal/dotfiles/blob/master/.config/python/requirements.txt)
    の形でまとめておくと楽。


<https://github.com/pyenv/pyenv/wiki/Common-build-problems>



## Anaconda

Scientificな用途で使いたい場合は
[Numpy/Scipy]({{< relref "scipy.md" >}})
などの主要パッケージもまとめて面倒みてくれる
[Anaconda](https://docs.continuum.io/anaconda/)
を使うという選択肢もある。
私は使わない。
GUIのインストーラでもいいし、Homebrewでも入れられる。

ただし`PATH`上でシステムコマンドを上書きしちゃうヤンチャな面もあるので、
それが気になる人はpyenv越しに入れることで汚染をある程度防げる。
全部入りに抵抗がある場合は
`pyenv install miniconda3-latest`
から小さくスタートすることも可能。
パッケージ管理では `pip` の代わりに非公式の `conda` を使うことになる。


## 環境変数

https://docs.python.org/3/using/cmdline.html#environment-variables

シェルの設定ファイル(`~/.zshrc` とか)で `export` しておく。
一時的に無効したいときは `python -E` で起動。

### `PYTHONPATH`

`import` の探索パス (`sys.path`) の先頭付近に場所を追加できる。
例えば自分で書いたモジュールやパッケージの置き場所を指定しておけば、
いつでも優先的に `import` できるようになる。


### `PYTHONUSERBASE`

[`pip install`]({{< relref "pip.md" >}}) や `setup.py install` における
`--user` オプションの目的地を指定できる。
デフォルトでは `${HOME}/.local` 。
MacのFramework buildでは `${HOME}/Library/Python/2.7` とかになる。

現在の設定は
[`site`](https://docs.python.org/3/library/site.html)
の `USER_BASE` で確認できる (`python -m site`)。
`USER_SITE` はその下の `{BASE}/lib/python*.*/site-packages` に配置され、
`sys.path` に含まれる。
除外したいときは `PYTHONNOUSERSITE` をセットするか `python -s` で起動。


### `PYTHONSTARTUP`

インタラクティブモードで起動するときに読み込むファイルを指定できる。
