+++
title = 'dadi'
[menu.main]
  parent = "bio"
+++

<https://bitbucket.org/RyanGutenkunst/dadi>

## Installation

<https://bitbucket.org/RyanGutenkunst/dadi/wiki/Installation>

1.  Python 2.7 で `dadi` 用の `virtualenv` を作って、中に入る:

        % virtualenv ~/.virtualenv/dadi
        % source ~/.virtualenv/dadi/bin/activate

2.  依存パッケージをインストール:

        % pip install -U numpy==1.9.2 scipy matplotlib

    {{%div class="note"%}}
`numpy==1.10.0` だと
`ValueError: Cannot operate with a folded Spectrum and an unfolded one.`
というエラーが出て使えなかったのでとりあえず
`numpy==1.9.2` を使うべし。
    {{%/div%}}

3.  `git` でリポジトリのクローンを取得:

        % git clone git@bitbucket.org:RyanGutenkunst/dadi.git

4.  リリース版をチェックアウト:

        % cd dadi/
        % git tag
        % git checkout 1.7.0

5.  インストール:

        % python setup.py build
        % python setup.py install

6.  テスト:

        % cd tests/
        % python run_tests.py
        % cd ..

7.  `doc/` 以下に結構詳しいドキュメントが入ってるので読む (これもウェブで公開してくれたらいいのに)

## 入力データ

<https://bitbucket.org/RyanGutenkunst/dadi/wiki/DataFormats>
