+++
title = 'Keyboard'
tags = ["mac", "writing"]
[menu.main]
  parent = "mac"
+++

## Default Shortcuts

スクリーンショット（全画面）
:   <kbd>shift-command-3</kbd>

スクリーンショット（マウスで範囲選択）
:   <kbd>shift-command-4</kbd>

スクリーンショット（撮りたいウィンドウをクリック）
:   <kbd>shift-command-4</kbd> <kbd>space</kbd>

Dialog (Restart/Sleep/Shut Down)
:   <kbd>control-eject</kbd>

Sleep Display
:   <kbd>shift-control-eject</kbd>

Sleep
:   <kbd>option-command-eject</kbd>

Restart
:   <kbd>control-command-eject</kbd>

Shut Down
:   <kbd>control-option-command-eject</kbd>


## Customize

US配列のほうがいろいろ合理的だし、余計な印字もなくてカッコイイ。
JIS配列のほうが優れている点
（左手小指Aの隣はControlの指定席、
スペース左右の「英数」「かな」キーは便利、
右手だけでForward Deleteできる）は、
以下のように変更すればUS配列でも実現できる。

System Preferences > Keyboard > Modifier Keys...
: <kbd>caps lock</kbd> to <kbd>^ control</kbd>

そのほか、Cocoaアプリのキーバインドは以下のファイルで定義できる。

-   `~/Library/KeyBindings/DefaultKeyBinding.dict`

### [Karabiner](https://pqrs.org/osx/karabiner/)

旧KeyRemap4MacBook

-   Change <kbd>Delete</kbd> Key
    - <kbd>Option+Delete</kbd> to Forward Delete
      (あるいは、<kbd>Option_R</kbd> を <kbd>fn</kbd> にしてしまってもいいかも)

-   For Japanese
    -   左右のコマンドキーを「英数/かな」としても使う


## ウムラウト、アクセント符号

<kbd>option</kbd> キーを使って修飾つきアルファベットを入力できる

Umlaut (ä ï ü ë ö ÿ)
:   <kbd>option-u</kbd> + 母音

Acute accent (á í ú é ó)
:   <kbd>option-e</kbd> + 母音

Grave accent (à ì ù è ò)
:   <kbd>option-`</kbd> + 母音

Circumflex (â î û ê ô)
:   <kbd>option-i</kbd> + 母音

Tilde (ã õ ñ)
:   <kbd>option-n</kbd> + 母音

## Entering Unicode Text and Symbols

1.  `System Preferences --> Language & Text --> Input Sources` で Unicode Hex Input にチェック
2.  <kbd>⌘command-space</kbd> で Unicode Hex Input を選択
3.  <kbd>⌥option</kbd> を押しつつ番号を入力 (e.g. <kbd>⌥option + 2318</kbd> で ⌘ が入力される)
