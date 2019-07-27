+++
title = 'setuptools'
subtitle = "Pythonパッケージ作成"
tags = ["python", "package"]
[menu.main]
  parent = "python"
+++

- https://docs.python.org/3/tutorial/modules.html
- https://docs.python.org/3/reference/import.html
- https://packaging.python.org/
- https://setuptools.readthedocs.io/

## ファイル構成

GitHubやローカルの開発環境から `pip` で直接インストールできる形。

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
ローカルからのインストールは `-e,--editable`
オプションを付けることで余計なコピーを減らせるっぽい。

```sh
pip install -v -e ~/git/pywtl/
pip install -v git+https://github.com/heavywatal/pywtl.git
python -m wtl.hello
```


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
version = attr: wtl.__version__
url = https://github.com/heavywatal/pywtl
author = Watal M. Iwasaki
author_email = heavy.watal@gmail.com
license_file = LICENSE
description = wtl: Personal Python package
long_description = file: README.md

[options]
zip_safe = False
packages = find:
install_requires =
  psutil
  requests
entry_points = file: entry_points.cfg
```

`version` を書くべき場所が複数あって悩ましかったけど
`attr:` や `file:` がサポートされて一括管理可能になった。
(Thanks, [@sgryjp](https://twitter.com/sgryjp)-san!)

`license = MIT` のように文字列で済ませることも可能。

依存パッケージを列挙する `install_requires` はセミコロン区切りで1行に収めることも可能。
それらは `pip install` で自動的にインストールされ、
その後アンインストールするとそれを使わないモジュールでさえ
`distribution was not found and is required`
などと怒って動いてくれなくなるくらい強制力が強い。

`requirements.txt` はsetuptoolsではなくpipの機能で、
能動的に `pip install -r requirements.txt`
を打たなきゃインストールされないし、
そこに列挙されたパッケージが本当に必要になるまで怒られない。


### `entry_points`

- https://setuptools.readthedocs.io/en/latest/setuptools.html#automatic-script-creation
- https://packaging.python.org/distributing/#entry-points

用途はいろいろあるけど
`${prefix}/bin/` に実行可能ファイルを配置するのによく使われる。
設定は下記のように別ファイル `entry_points.cfg` として
`setup.cfg` から読ませるのが楽チン。

```ini
[console_scripts]
hello.py = wtl.hello:main
```

引数を取らない関数のみ利用可能。
コマンドライン引数を受け取りたい場合はその関数の中で標準の
[`argparse`](https://docs.python.org/3/library/argparse.html)
を使って処理する。


### ソースコード

`wtl/__init__.py`
: このディレクトリがひとつのパッケージであることを示すファイル。
  空でもいいし、初期化処理やオブジェクトを記述してもよい。
  文字列変数 `__version__ = '0.1.2'` を定義して
  `wtl.__verison__` のように参照できるようにしておくのが慣例。
  それを `setup.cfg` から参照することで冗長性を抑えられる。

`wtl/hello.py`
```py
"""Sample module
"""


def main():
    import getpass
    print('Hello, ' + getpass.getuser() + '!')


if __name__ == '__main__':
    main()
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
