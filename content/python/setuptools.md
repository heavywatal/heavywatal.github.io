+++
title = 'setuptools'
subtitle = "Pythonパッケージ作成"
tags = ["python", "package"]
[menu.main]
  parent = "python"
+++

- https://packaging.python.org/
- https://setuptools.readthedocs.io/

## ファイル構成

GitHubから `pip` で直接インストールできるようにしたい。

```sh
% pip install git+https://github.com/heavywatal/pywtl.git
% python -c 'import wtl; wtl.hello()'
```

```sh
pywtl/
├── LICENSE
├── README.md
├── entry_points.cfg
├── setup.cfg
├── setup.py
└── wtl/
    ├── __init__.py
    └── hello.py
```

リポジトリ名(`pywtl`)とパッケージ名(`wtl`)は必ずしも一致してなくてもよい。


### `setup.py`

新しいsetuptoolsでは `setup.cfg` の設定を読み込んでくれるので、
`setup()` の引数にいろいろ渡す必要は無くなった。

```py
from setuptools import setup
setup()
```

### `setup.cfg`

https://setuptools.readthedocs.io/en/latest/setuptools.html#configuring-setup-using-setup-cfg-files

```ini
[metadata]
name = wtl
version = 0.1
url = https://github.com/heavywatal/pywtl
author = Watal M. Iwasaki
author_email = heavy.watalあgmail.com
license = MIT
description = wtl: Personal Python package
long_description = file: README.md

[options]
zip_safe = False
packages = find:
entry_points = file: entry_points.cfg
```

`license = file: LICENSE` のように外部ファイルを参照することも可能。

### `entry_points`

- https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
- https://packaging.python.org/distributing/#entry-points

用途はいろいろあるけど
`${prefix}/bin/` に実行可能スクリプトを配置するのによく使われる。
設定は下記のように別ファイル `entry_points.cfg` として
`setup.cfg` から読ませるのが楽チン。

```ini
[console_scripts]
hello.py = wtl.hello:main
```

引数を取らない関数のみ利用可能。
コマンドライン引数を処理したい場合は `argparse` が有効。


### ソースコード

`wtl/__init__.py`
: このディレクトリがひとつのパッケージであることを示すファイル。
  空でもいいし、初期化処理やオブジェクトを記述してもよい。

`wtl/hello.py`
```py
def main():
    print('Hello, world!')
```

### `MANIFEST.in`

Pythonモジュール以外で配布物に含めたいもの、除外したいものがあれば
`include`, `recursive-include`, `exclude`
で指定する。


## setuptools

[Command Reference](https://setuptools.readthedocs.io/en/latest/setuptools.html#command-reference)

```sh
python setup.py --help
python setup.py --help-commands
python setup.py ^i  # zsh completion

# 開発中に使う
python setup.py build
python setup.py clean
python setup.py check

# 配布できる形に固める
python setup.py sdist
python setup.py bdist_wheel
```