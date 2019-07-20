+++
title = '開発環境'
tags = ["package"]
[menu.main]
  parent = "dev"
+++

https://github.com/heavywatal/dotfiles

## Linux

### ユーザ追加

ユーザーを作って管理権限を与える:

    # adduser USERNAME
    # gpassed -a USERNAME sudo

### すぐ入れるパッケージ

`sudo apt-get install`:

    build-essential
    zsh
    tmux
    git

これらシステム標準のものが古すぎたり、
管理者権限がなくて自由にインストールできない場合は
[Homebrew]({{< relref "homebrew.md" >}}) を利用してユーザのホームに入れる。

## Mac

### Command Line Tools

コンパイラや `make` などはOSに付いてこないので別途インストールが必要。
<https://developer.apple.com/downloads/> からダウンロードするか、
ターミナルから以下のコマンドを実行:

```sh
xcode-select -p
xcode-select --install
```

インストールされているバージョンなどを確認するには:
```sh
pkgutil --pkg-info=com.apple.pkg.CLTools_Executables
clang -v
```

総合開発環境 Xcode をインストールしたければ、App Store から [Xcode](https://itunes.apple.com/jp/app/xcode/id497799835) を選択。

### パッケージ管理ツール

-   [/mac/homebrew]({{< relref "homebrew.md" >}})
-   [/mac/macports]({{< relref "macports.md" >}})

### その他のプログラム

-   [MenuMeters](https://member.ipmu.jp/yuji.tachikawa/MenuMetersElCapitan/)
-   [QuickLook plugins]({{< relref "quicklook.md" >}})
-   [`defaults`コマンドで各種設定]({{< relref "command.md#defaults" >}})


## 共通

### Python

- [install]({{< relref "/python/install.md" >}})
- [pip]({{< relref "pip.md" >}})
- `Pillow` をインストールする前に:

      sudo apt-get install libtiff5-dev libwebp-dev libfreetype6-dev liblcms2-dev libopenjpeg-dev

### C++

- [boost]({{< relref "boost.md" >}})
- [SFMT]({{< relref "random.md" >}})

### R

[/rstats/config]({{< relref "/rstats/config.md" >}})

### エディタ

- [atom]({{< relref "atom.md" >}})
- [emacs]({{< relref "emacs.md" >}})
- [vim]({{< relref "vi.md" >}})
- [nano]({{< relref "nano.md" >}})

### Trash

`rm` はゴミ箱を経由せず削除してしまうので、
間違って消してしまっても基本的には元に戻せない。
以下のような対策によりその危険が少しは減るかも。

`alias rmi='rm -i'`
:   ホントに消していいかどうか確認してくれるようなオプションつきのエイリアスを
    `.zshrc` に設定しておく。

        rmi .DS_Store
        rm: remove regular file '.DS_Store'?

    エイリアス名を `rm` そのものにしてしまうと、
    結局ろくすっぽ確認せず `y` を押す癖や、
    いちいち確認されないように `rm -rf` する癖がつくので逆に危険。
    普段は `rmi` を使う癖をつけ、必要なときたまに `rm` を使い、
    `-f` はよほどのことが無い限り使わないようにする。

`trash-cli`
:   Python製なので [pip]({{< relref "pip.md" >}}) で
    `pip install trash-cli` して入れる。
    すると以下のようなコマンドがインストールされる:

        restore-trash
        trash
        trash-empty
        trash-list
        trash-put
        trash-rm

    ちなみにLinuxのゴミ箱は `~/.local/share/Trash`

`rmtrash`
:   Macでも `trash-cli` を使えないことはないが、
    ゴミ箱のパスがMac標準の `~/.Trash` ではなく
    Linuxのものになってしまうので
    [Homebrew]({{< relref "homebrew.md" >}}) で `rmtrash` を入れる:

        brew install rmtrash
