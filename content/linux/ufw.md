+++
title = 'ufw'
[menu.main]
  parent = "linux"
+++

Linux のファイアウォール `iptables` をもっと簡単に設定できるようにしたもの

<https://help.ubuntu.com/community/UFW>

## Usage

-   オン、オフ:

        % sudo ufw enable
        % sudo ufw disable

-   現状確認:

        % sudo ufw status verbose
        % sudo ufw show raw

-   設定例:

        # eth1
        sudo ufw allow in on eth1
        #sudo ufw allow from 192.168.0.0/24

        # ssh
        sudo ufw allow from 130.34.107.0/25 to any port ssh
        sudo ufw limit ssh

        # cifs/samba
        #sudo ufw allow from 130.34.107.0/25 to any app samba
        sudo ufw allow proto udp from 130.34.107.0/25 to any port 137
        sudo ufw allow proto udp from 130.34.107.0/25 to any port 138
        sudo ufw allow proto tcp from 130.34.107.0/25 to any port 139
        sudo ufw allow proto tcp from 130.34.107.0/25 to any port 445

        # avahi/bonjour/zeroconf
        sudo ufw allow proto udp to any port 5353 from 130.34.107.0/25

        # http
        sudo ufw allow 80/tcp

        # gema 626
        # iMac 8612
        sudo ufw allow proto udp from 130.34.107.0/25

-   `status` で上に表示される(=番号が若い)ほど優先順位が高いので指定順に注意。
    並べ替えはできないので `reset` で全消しして始めからやり直すか、
    `delete` や `insert` で調整:

        % sudo ufw status numbered
        % sudo ufw delete 1
        % sudo ufw insert 3 allow 80/tcp
        % sudo ufw reset

## Installation

    % sudo apt-get install ufw
