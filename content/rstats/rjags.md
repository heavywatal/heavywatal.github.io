+++
title = 'rjags'
subtitle = "JAGSを使ってベイズ統計"
[menu.main]
  parent = "rstats"
+++

-   <http://mcmc-jags.sourceforge.net/>
-   <http://www.rdocumentation.org/packages/rjags>

JAGS (Just Another Gibbs Sampler)
はギブズサンプリングをやってくれるCプログラム。
それをRから動かすためのパッケージが `rjags` 。
データを R で整えて JAGS に渡し、
出て来た結果を R で受け取って解析する。

## Installation

### 普通にコンパイル済みのやつを入れる

1.  <http://sourceforge.net/projects/mcmc-jags/files/>
    からOSに合ったインストールファイルをダウンロード。
2.  ダウンロードしたファイルを実行し、「はい」とか「次へ」で普通にインストール。
    これでJAGS本体のインストール完了。
3.  Rを立ち上げて `rjags` というパッケージをインストール。
    画面上部のメニューからインストーラを立ち上げるか、以下のコマンドをコピペ

    ```r
    install.packages("rjags", dependencies=TRUE)
    ```

    Rのバージョンによって相性があるようなので、問題があったらそのへんをいじってみるといいかも。

### Homebrew で入れる

1.  [Homebrew]({{< relref "mac/homebrew.md" >}}) をインストール
2.  JAGS本体をインストール

    ```sh
    % brew install jags
    ```

    ソースコードからビルドしてもよい

    ```sh
    % wget -O- http://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/Source/JAGS-3.4.0.tar.gz | tar xz
    % cd JAGS-3.4.0/
    % ./configure --help
    % ./configure --prefix=/usr/local/jags
    % make
    % sudo make install
    ```

3.  R を立ちあげ、JAGSの場所を指定して
    `rjags` をソースコードからインストール

    ```r
    prefix = system('brew --prefix', TRUE)
    template = "--with-jags-include=%s/include/JAGS --with-jags-lib=%s/lib"
    args = sprintf(template, prefix, prefix)
    install.packages("rjags", type="source", configure.args=args)
    ```
