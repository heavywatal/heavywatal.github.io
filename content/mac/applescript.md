+++
title = 'AppleScript'
tags = ["mac"]
[menu.main]
  parent = "mac"
+++

## メニューバーから実行可能にする

- `~/Library/Scripts/` にスクリプトを置く
- `Script Editor.app` を起動
- <kbd>cmd</kbd><kbd>,</kbd> Preferences...
- ✅ Show Script menu in menu bar

## Music Equalizer

イコライザを設定するスクリプト。
`band 1` が低音で `band 10` が高音。
下記の設定は Equal-loudness contour(等ラウドネス曲線) を考慮して
小さい音量でもシャカシャカせずそれなりのバランスで聴こえるようにした例。
耳に届きやすい 4kHz (band 8) らへんを落とすのが肝。
機器や音量や好みによって調整し、別々のスクリプトとして保存しておくとよい:

```applescript
tell application "Music"
    tell EQ preset "Manual"
        set band 1 to 1.5
        set band 2 to 0.7
        set band 3 to 0.3
        set band 4 to 0.1
        set band 5 to -0.1
        set band 6 to -0.5
        set band 7 to -1
        set band 8 to -2
        set band 9 to 0.3
        set band 10 to 0.7
        set preamp to 0
    end tell
    set current EQ preset to EQ preset "Manual"
end tell
```

## Finder: Get file path

```applescript
tell application "Finder"
    set the clipboard to POSIX path of (the selection as alias)
end tell
```

## osascript: スクリプトをターミナルで実行

スクリプトを置いたディレクトリにパスを通せれば楽なんだけど。できるのかな？:

```sh
osascript ${HOME}/Library/Scripts/Applications/Bibdesk/exportOOXML.scpt
```

## Bibdesk からGroupとTemplateを指定してExport

出力先、グループ名、テンプレートを最初の３行で指定し、Exportしたいライブラリを開いた状態で実行:

```applescript
set theOutFile to POSIX file "/Users/Iwasaki/Presentations/20100118Thesis/items1.xml"
set theGroup to "20100118 Thesis"
set theTemplateName to "OOXML template"

using terms from application "BibDesk"
    tell application "BibDesk"
        set theDoc to first document
        set thePubs to publication of static group theGroup of theDoc
        tell theDoc
            export for thePubs using template theTemplateName to theOutFile
        end tell
    end tell
end using terms from
```
