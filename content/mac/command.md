+++
title = 'Mac Command'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

## `open`

Finder でダブルクリックするのと同じように、
Terminal から一発でファイルを開くことができる。

```sh
#### 関連付けられたデフォルトのアプリケーションで開く
% open hudson1992g.pdf

#### アプリケーションを指定して開く
% open -a Skim hudson1992g.pdf

#### TextEdit.app で開く
% open -e README

#### デフォルトのテキストエディタで開く
% open -t README

#### カレントディレクトリをFinderで開く
% open .
```

## `defaults`

`/Library/Preferences/` 以下にある各種.plistを編集し、
いろいろな設定を変更するコマンド。
`true` に設定した項目を元に戻すには、
項目自体を `delete` するか、`false` に設定する。

    % defaults [write/delete] DOMAIN KEY -TYPE VALUE

    ### Finderで隠しファイルを見えるようにする。
    % defaults write com.apple.finder AppleShowAllFiles -bool true

    ### 元に戻すには下のいずれか
    % defaults write com.apple.finder AppleShowAllFiles -bool false
    % defaults delete com.apple.finder AppleShowAllFiles

とは言え
[Onyx](http://www.titanium.free.fr) や
[Tinkertool](http://www.bresink.com/osx/TinkerTool.html)
などのGUIアプリを使うほうが簡単で安心かも

iTerm2のタブの横幅を広くする
:   `defaults write com.googlecode.iterm2 OptimumTabWidth -int 360`

## `killall`

Finder や Dock など、GUIから終了させにくいアプリケーションを再起動する。
アプリケーションの動作が不安定になったとき、設定変更を反映させたいとき、
メモリを開放したいときなどに使える。:

    % killall Finder
    % killall Dock
    % killall Kotoeri

## `lsregister`

Open with で表示されるアプリケーションが重複しまくったときなど、
ファイルとアプリケーションの関連付けに関する古い情報を消して再構築。
まっさらに戻るわけではない。:

    % /System/Library/Frameworks/CoreServices.framework/Frameworks/LaunchServices.framework/Support/lsregister -kill -r -domain local -domain system -domain user

設定ファイルは `~/Library/Preferences/com.apple.LaunchServices.plist`

## インストール関連

`.dmg` のマウント、`.pkg` からのインストール、
システムのソフトウェア・アップデートなどを
`ssh` 越しにやらねばならぬときもある:

    % hdiutil mount SomeDiskImage.dmg
    % sudo installer -pkg SomePackage.pkg -target /
    % sudo softwareupdate -i -a

## アカウント管理: `niutil`, `nidump`

{{%div class="note"%}}
Leopard以前の古いOSでしか使えない
{{%/div%}}

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
