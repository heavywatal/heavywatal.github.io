+++
title = 'Keyboard'
tags = ["mac", "writing"]
[menu.main]
  parent = "mac"
+++

## Default Shortcuts

`shift + command + 3`
:   スクリーンショット（全画面）

`shift + command + 4`
:   スクリーンショット（範囲選択）

`shift + command + 4 + space`
:   スクリーンショット（ウィンドウ選択）

`control + eject`
:   Dialog (Restart/Sleep/Shut Down)

`shift + control + eject`
:   Sleep Display

`option + command + eject`
:   Sleep

`control + command + eject`
:   Restart

`control + option + command + eject`
:   Shut Down

## Customize

US配列のほうがいろいろ合理的だし、余計な印字もなくてカッコイイ。
JIS配列のほうが優れている点
（左手小指Aの隣はControlの指定席、
スペース左右の「英数」「かな」キーは便利、
右手だけでForward Deleteできる）は、
以下のように変更すればUS配列でも実現できる。

`System Preferences --> Keyboard --> Modifier Keys...`
: `Caps Lock` `^ Control`

そのほか、Cocoaアプリのキーバインドは以下のファイルで定義できる。

-   `${HOME}/Library/KeyBindings/DefaultKeyBinding.dict`

### [Karabiner](https://pqrs.org/osx/karabiner/)

旧KeyRemap4MacBook

-   Change `Delete` Key
    - `Option+Delete` to Forward Delete
      (あるいは、`Option_R` を `fn` にしてしまってもいいかも)

-   For Japanese
    -   左右のコマンドキーを「英数/かな」としても使う


## ウムラウト、アクセント符号

`option` キーを使って修飾つきアルファベットを入力できる

Umlaut (ä ï ü ë ö ÿ)
:   `option-u` + 母音

Acute accent (á í ú é ó)
:   `option-e` + 母音

Grave accent (à ì ù è ò)
:   `option-`\` + 母音

Circumflex (â î û ê ô)
:   `option-i` + 母音

Tilde (ã õ ñ)
:   `option-n` + 母音

## Entering Unicode Text and Symbols

1.  `System Preferences --> Language & Text --> Input Sources` で Unicode Hex Input にチェック
2.  `⌘command + space` で Unicode Hex Input を選択
3.  `⌥option` を押しつつ番号を入力 (e.g. `⌥option + 2318` で ⌘ が入力される)
