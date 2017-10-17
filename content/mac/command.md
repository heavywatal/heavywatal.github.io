+++
title = 'Mac固有コマンド'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

## 簡単・便利

### `open`

Finder でダブルクリックするのと同じように、
Terminal から一発でファイルを開くことができる。

```sh
# 関連付けられたデフォルトのアプリケーションで開く
% open hudson1992g.pdf

# アプリケーションを指定して開く
% open -a Skim hudson1992g.pdf

# カレントディレクトリをFinderで開く
% open .
```

### `pbcopy`, `pbpaste`

ターミナルからクリップボード(pasteboard)を操作する。

```sh
% pwd | pbcopy
% pbpaste | wc
```

[tmux内で使うにはreattach-to-user-namespaceが必要。]({{< relref "tmux.md" >}})


### `say` 音声読み上げ

```sh
% say hello world
% pbpaste | say
```

声は環境設定の Accessibility > Speech で変更可能。


### `killall`

Finder や Dock など、GUIから終了させにくいアプリケーションを再起動する。
アプリケーションの動作が不安定になったとき、設定変更を反映させたいとき、
メモリを開放したいときなどに使える。:

```sh
% killall Finder
% killall Dock
% killall Kotoeri
```


## 設定関連

### `defaults`

`/Library/Preferences/` 以下にある各種設定ファイル.plistを編集する。
`true` に設定した項目を元に戻すには、
項目自体を `delete` するか、`false` に設定する。

```sh
defaults [write/delete] DOMAIN KEY -TYPE VALUE

# Finderのタイトルバーにフルパスを表示
defaults write com.apple.finder _FXShowPosixPathInTitle -boolean true

# iTerm2のタブの横幅を広くする
defaults write com.googlecode.iterm2 OptimumTabWidth -int 360

# Launchpadの並び順をリセット
defaults write com.apple.dock ResetLaunchPad -bool true; killall Dock

# 使わないXcodeのための重いindexingを切っておく
defaults write com.apple.dt.Xcode IDEIndexDisable 1

# Quicklook上でコピペできるようにする
defaults write com.apple.finder QLEnableTextSelection -bool true
```

[Onyx](http://www.titanium.free.fr) や
[Tinkertool](http://www.bresink.com/osx/TinkerTool.html)
などのGUIアプリを使うほうが簡単で安心かも


### `launchctl`

サービスの起動、終了


### `lsregister`

Open with で表示されるアプリケーションが重複しまくったときなど、
ファイルとアプリケーションの関連付けに関する古い情報を消して再構築。
まっさらに戻るわけではない。:

    % /System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister -kill -r -domain local -domain system -domain user

設定ファイルは `~/Library/Preferences/com.apple.LaunchServices.plist`

### インストールなど

`.dmg` のマウント、`.pkg` からのインストール、
システムのソフトウェア・アップデートなどを
`ssh` 越しにやらねばならぬときもある:

    % hdiutil mount SomeDiskImage.dmg
    % sudo installer -pkg SomePackage.pkg -target /
    % sudo softwareupdate -i -a

## Obsolete

### アカウント管理: `niutil`, `nidump`

Leopard以前

項目をリストアップ:

    % niutil -list . /
    % niutil -list . /users
    % niutil -list . /groups

中身を見る:

    % niutil -read . /users/iwasaki
    % niutil -read . /groups/admin

一覧で一気に:

    % nidump passwd . /
    % nidump group . /

新規ユーザーの追加:

    % niutil -create / /users/hoge
    % niutil -createprop / /users/hoge shell /bin/zsh
    % niutil -createprop / /users/hoge uid 1050
    % niutil -createprop / /users/hoge gid 20
    % niutil -createprop / /users/hoge home /Users/hoge
    % niutil -createprop / /users/hoge _shadow_passwd
    % passwd hoge
    % niutil -appendprop / /groups/staff users hoge
