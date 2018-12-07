+++
title = 'dadi'
tags = ["genetics", "python"]
[menu.main]
  parent = "bio"
+++

https://bitbucket.org/gutenkunstlab/dadi

## Installation

https://bitbucket.org/gutenkunstlab/dadi/wiki/Installation

1.  macOS 標準 Python 2.7 に
    [pip と virtualenv]({{< relref "pip.md" >}}) を入れる:

        % /usr/bin/python -m ensurepip -v --user
        % ~/Library/Python/2.7/bin/pip install --user -U setuptools pip virtualenv

1.  dadi専用のvirtualenvを作って、その中に依存パッケージをインストール:

        % ~/Library/Python/2.7/bin/virtualenv ~/.virtualenv/dadi
        % source ~/.virtualenv/dadi/bin/activate
        (dadi) % pip install -U setuptools pip flake8
        (dadi) % pip install -U numpy scipy matplotlib

1.  リポジトリのクローンを取得し、リリース版をチェックアウトしてインストール:

        (dadi) % git clone https://bitbucket.org/gutenkunstlab/dadi.git
        (dadi) % cd dadi/
        (dadi) % git tag
        (dadi) % git checkout 1.7.0
        (dadi) % python setup.py install

1.  テスト:

        (dadi) % cd tests/
        (dadi) % python run_tests.py
        (dadi) % cd ..

1.  `doc/` 以下に結構詳しいドキュメントが入ってるので読む (これもウェブで公開してくれたらいいのに)


## 入力データ

https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats

## Tips

- $\theta = 4N_0\mu$ は ms と同じだが、
  時間とmigration rateは $2N_0$ 単位。
- 時間とmigration rateは大きすぎると遅くなるので
  upper bound はそれぞれ 10, 20 を推奨
- population size は逆に小さいとき遅くなるので
  lower bound `1e-3` 程度を推奨

