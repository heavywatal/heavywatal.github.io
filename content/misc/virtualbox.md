+++
title = 'VirtualBox'
[menu.main]
  parent = "misc"
+++

## Free virtualization software

<http://www.virtualbox.org/>

現在使っている条件

-   ホストOS: Mac OS X 10.7 Lion
-   ゲストOS: Windows 7 Ultimate 64-bit

<http://www.virtualbox.org/manual/>

## Installation

1.  VirtualBoxをダウンロードしてインストール
2.  VirtualBoxを立ち上げ、環境設定 `command + ,`
    -   `General --> Default Machine Folder`
        を適当なところに (e.g. `~/Library/VirtualBox`)
    -   `Input --> Host Key` を適当なキーに (e.g. `右option`)

3.  ツールバーのNewからゲストOSの設定を作成

    Name
    :   Win7u64

    Memory
    :   2048

    Virtual Hard Disk
    :   Create new

    -   Dynamically expanding storage
    -   Location: 適当なところに置き、Time Machineバックアップから除外する
    -   Size: デフォルト20GBでも足りなくはないけど、OSアップデートやソフトウェアの追加を考えるともうちょいあったほうが安心かも

4.  ゲストOSをインストールする前にいくつかの設定を変更
    -   `System > Processor`: 2
    -   `Display > Video > Video Memory`: 64MBくらい？

5.  設定したVirtual Machineをスタートし、インストールディスクから普通にOSインストール
6.  続けて、下記の手順でゲストOSに Guest Additions\_ をインストールする
7.  ゲストOSのセキュリティアップデートやライセンス認証などを済ます

## Guest Additions

<http://www.virtualbox.org/manual/ch04.html>

ゲストOS側にこれをインストールすることでホストOSとの親和性が高まり、
以下のように便利になるのでとりあえず入れよう

-   Mouse Integrationが有効になり、ホストOSとゲストOSのマウスの動きが連動するようになる
-   ホストOSとゲストOSの共有フォルダをつくれるようになる
-   クリップボードの共有が可能になる
-   Seamless modeを使えるようになる `Host + L`
-   ゲストOSの時刻がホストOSに合う

インストール手順

1.  ゲストOSをセーフモードで起動（VM起動直後に `F8` キー連打）
2.  ホストOSメニューバーの `Devices > Install Guest Additions...`
3.  ゲストOSのCD/DVDドライブに読み込まれるインストーラを実行

{{%div class="note"%}}
CentOSの場合はまず
`gcc`, `make`, `kernel-devel`
をインストールして再起動してから。
{{%/div%}}

## Network

以下のようなコマンドでIPマスカレードを設定する。
`~/Library/VirtualBox/***/***.vbox` というXMLファイルに書き込まれる:

    % VBoxManage setextradata "***" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/***/Protocol" ***

### ssh

    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/ssh/Protocol" TCP
    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/ssh/GuestPort" 22
    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/ssh/HostPort" 60022

    % ssh -p 60022 user@localhost

### VNC

    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/vnc/Protocol" TCP
    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/vnc/GuestPort" 5906
    % VBoxManage setextradata "LinuxMint" "VBoxInternal/Devices/pcnet/0/LUN#0/Config/vnc/HostPort" 5906

    % vnc4server :6

    % vnc://localhost:5906
