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
`/usr/bin/python` が既にインストールされているので、
何もしなくても使えっちゃ使える。
でも大概そこに入ってるのは古い2.7とかなので、
ちゃんと使える3.xを使いたければ
[python.jp/環境構築ガイド](https://www.python.jp/install/install.html)
に従って最新版を入れるのが簡単。


## rye

<https://rye.astral.sh/>

プロジェクトの環境構築を一切合切面倒見てくれる管理ツール。
ruffやuvと同じチームが開発していて、同じくrust製。
次のワンライナーでプログラム本体や設定ファイルなどが `~/.rye/` に配置される:
```sh
curl -sSf https://rye.astral.sh/get | bash
```

任意のバージョンのPythonを入れるためのツールとして、
つまり[pyenv](#pyenv)的な位置付けでも使える。
しかもビルド済みのPythonを
<https://github.com/indygreg/python-build-standalone>
から取ってくるので、自前ビルド環境に左右されずCPUも使わず簡単。

### Pythonインストーラーとしてだけ使う

パッケージ運用はまだ普通にpipとかでいいかなと思うので、今のところ私はこの使い方。

```sh
# インストール済みバージョン一覧
rye toolchain list

# 利用可能バージョン一覧
rye toolchain list --include-downloadable

# インストール
rye toolchain fetch cpython@3.12.3
```

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
ソースコードを取ってきて自前ビルドするのが上記[rye](#rye)と比べたときのメリットでありデメリット。

1.  [Homebrew]({{< relref "homebrew.md" >}}) か
    [Git]({{< relref "git.md" >}}) を使ってpyenvをインストール:

    ```sh
    brew install pyenv
    # or
    git clone https://github.com/pyenv/pyenv.git ~/.pyenv
    mkdir -p ~/.pyenv/cache
    ```

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
    pyenv install 3.12.3
    exec $SHELL -l
    ```

    R から [`reticulate`](https://rstudio.github.io/reticulate/)
    越しに呼ぶ場合は共有ライブラリを有効にしてビルドする:
    ```sh
    env PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv install 3.7.13
    # or
    Rscript -e 'reticulate::install_python("3.7.13")'
    ```


1.  [pip]({{< relref "pip.md" >}}) のパスを確認し、パッケージを入れる:

    ```sh
    which pip3
    pip3 install -U setuptools pip wheel
    pip3 install -r /path/to/requirements.txt
    ```

    よく使うパッケージは
    [`requirements.txt`](https://github.com/heavywatal/dotfiles/blob/master/.config/python/requirements.txt)
    の形でまとめておくと楽。


### 既知の問題

<https://github.com/pyenv/pyenv/wiki/Common-build-problems>

-   3.1.0より古い[matplotlib]({{< relref "matplotlib.md" >}}) で
    `macosx` backend を使いたい場合などは環境変数
    `PYTHON_CONFIGURE_OPTS="--enable-framework"`
    をセットしてFramework型でビルドする。

-   [Mojaveで "zlib not available" と怒られる問題](https://github.com/pyenv/pyenv/issues/1219)は
    `CFLAGS="-I$(xcrun --show-sdk-path)/usr/include"` を定義して回避。

-   [PEP 394](https://www.python.org/dev/peps/pep-0394/)
    になるべく沿うように、
    `python` では `/usr/bin/python` が呼び出される状態を維持しつつ、
    `python3`, `pip3` を明示的に使うようにしたい。
    が、いまのところできない？


## Anaconda

Scientificな用途で使いたい場合は
[Numpy/Scipy]({{< relref "scipy.md" >}})
などの主要パッケージもまとめて面倒みてくれる
[Anaconda](https://docs.continuum.io/anaconda/)
を使うという選択肢もある。
私は使わない。
GUIのインストーラでもいいし、Homebrewでも入れられる:

```sh
brew install anaconda
export PATH=/usr/local/anaconda3/bin:"$PATH"
```

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
