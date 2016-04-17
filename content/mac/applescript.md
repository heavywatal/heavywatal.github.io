+++
title = 'AppleScript'
[menu.main]
  parent = "mac"
+++

## Finder: Get file path

    tell application "Finder"
        set the clipboard to POSIX path of (the selection as alias)
    end tell

## iTunes Equalizer

イコライザを設定するスクリプト。
`band 1` が低音で `band 10` が高音。
下記の設定は Equal-loudness contour(等ラウドネス曲線) を考慮して
小さい音量でもシャカシャカせずそれなりのバランスで聴こえるようにした値。
音量やスピーカや好みによって調整し、別々のスクリプトとして保存しておくとよい:

    tell application "iTunes"
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
    end tell

## osascript: スクリプトをターミナルで実行

スクリプトを置いたディレクトリにパスを通せれば楽なんだけど。できるのかな？:

    % osascript ${HOME}/Library/Scripts/Applications/Bibdesk/exportOOXML.scpt

## mi で選択中の文字列をターミナルで実行

    tell application "mi"
            tell document 1
                    set CommandLines to selection
            end tell
    end tell

    tell application "Terminal"
            activate
            ignoring application responses
                    do script with command CommandLines in window 1
            end ignoring
    end tell

## Bibdesk からGroupとTemplateを指定してExport

出力先、グループ名、テンプレートを最初の３行で指定し、Exportしたいライブラリを開いた状態で実行:

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
