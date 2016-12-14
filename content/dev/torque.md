+++
title = 'TORQUE'
tags = ["job", "concurrent"]
[menu.main]
  parent = "dev"
+++

<http://www.adaptivecomputing.com/products/open-source/torque/>

環境
: Ubuntu 12.04 LTS Precise Pangolin
: Torque 2.4.16


## Usage

[Commands Overview](http://www.adaptivecomputing.com/resources/docs/torque/2-5-9/a.acommands.php)

Pythonのvirtualenvは計算ノードに引き継がれないので、
メインのコマンドの前に明示的にactivateする必要あり:

    echo ". ~/.virtualenv/project/bin/activate && ./a.out" | /usr/bin/qsub [options]

## Installation

<http://tech.ckme.co.jp/torque.shtml>\
<http://www.clusterresources.com/torquedocs/>\
<http://www.adaptivecomputing.com/resources/docs/torque/2-5-9/index.php>

### Compute node

1.  `torque-mom` と `torque-client` をインストール:

        % sudo apt-get install torque-mom torque-client

2.  `/etc/hosts` を確認。
    ヘッドノードのホスト名が解決できるように:

        192.168.0.254 charles

    ヘッドノードが複数のインターフェイスを持つ場合は
    計算ノードから見えるIPを記述した行が一番上になるようにする。
    そうしないと `qsub` でキューを出しても実行されず、
    以下のような警告が静かに書き出される:

        % sudo momctl -d3
        WARNING: invalid attempt to connect from server 192.168.0.254:1022  (server not authorized)

        % less /var/spool/torque/mom_logs/2012xxxx
        pbs_mom;Job;process_request;request type QueueJob from host charles rejected (host not authorized)
        pbs_mom;Req;req_reject;Reject reply code=15008(Access from host not allowed, or unknown host MSG=request not authorized), aux=0, type=QueueJob, from PBS_Server@charles

3.  ヘッドノードのホスト名を `/etc/torque/server_name` に設定: `charles`
4.  `/var/spool/torque/mom_priv/config` に設定記述 (ref. [MOM Configuration](http://www.adaptivecomputing.com/resources/docs/torque/2-5-9/a.cmomconfig.php))

    `/home` を NFSマウントしてるので `scp` しないように:

        $usecp *:/home /home

    ログファイルの末尾に計算ノードの名前を追加:

        $log_file_suffix %h

5.  `pbs_mom` を再起動:

        % sudo service torque-mom restart

### Head node

1.  `torque-server` をインストール:

        % sudo apt-get install torque-server

2.  `pbs_server` が動いてしまっていたら `kill`:

        % ps aux | grep pbs

3.  `/etc/hosts` を確認。
    計算ノードのホスト名が解決できるように:

        192.168.0.1 node01
        192.168.0.2 node02
        192.168.0.3 node03

4.  ヘッドノードのホスト名を `/etc/torque/server_name` に設定: `charles`
5.  計算ノードの情報を `/var/spool/torque/server_priv/nodes` に書き込む:

        node01 np=4
        node02 np=4
        node03 np=4

    {{%div class="note"%}}
`np` はコア数
    {{%/div%}}

6.  `/var/spool/torque/spool` のパーミションを設定:

        % sudo chmod 777 /var/spool/torque/spool/
        % sudo chmod o+t /var/spool/torque/spool/

7.  `pbs_server` を新規作成:

        % sudo pbs_server -t create

8.  `pbs_server` の設定 (ref. [Server Parameters](http://www.adaptivecomputing.com/resources/docs/torque/2-5-9/a.bserverparameters.php)):

        % sudo qmgr -c "create queue batch queue_type=execution"
        % sudo qmgr -c "set queue batch priority=0"
        % sudo qmgr -c "set queue batch resources_default.nodes=1"
        % sudo qmgr -c "set queue batch resources_default.ncpus=1"
        % sudo qmgr -c "set queue batch resources_default.walltime=96:00:00"
        % sudo qmgr -c "set queue batch keep_completed=3600"
        % sudo qmgr -c "set queue batch enabled=true"
        % sudo qmgr -c "set queue batch started=true"

        % sudo qmgr -c "create queue low queue_type=execution"
        % sudo qmgr -c "set queue low priority=-255"
        % sudo qmgr -c "set queue low resources_default.nodes=1"
        % sudo qmgr -c "set queue low resources_default.ncpus=1"
        % sudo qmgr -c "set queue low resources_default.walltime=96:00:00"
        % sudo qmgr -c "set queue low keep_completed=3600"
        % sudo qmgr -c "set queue low enabled=true"
        % sudo qmgr -c "set queue low started=true"

        % sudo qmgr -c "create queue high queue_type=execution"
        % sudo qmgr -c "set queue high priority=255"
        % sudo qmgr -c "set queue high resources_default.nodes=1"
        % sudo qmgr -c "set queue high resources_default.ncpus=1"
        % sudo qmgr -c "set queue high resources_default.walltime=96:00:00"
        % sudo qmgr -c "set queue high keep_completed=3600"
        % sudo qmgr -c "set queue high enabled=true"
        % sudo qmgr -c "set queue high started=true"

        % sudo qmgr -c "set server default_queue=batch"
        % sudo qmgr -c "set server scheduling=true"
        % sudo qmgr -c "set server allow_node_submit = True"

9.  `pbs_server` を再起動:

        % sudo service torque-server restart
        % sudo service torque-scheduler restart

10. 計算ノードの認識を確認 (`state = free`):

        % pbsnodes -a

11. 試しにジョブを投入:

        % echo "sleep 30" | qsub

12. 実行されているか確認 (`S` の下が `R` なら走ってる):

        % qstat

## Maintenance / Trouble Shooting

### Head node

設定確認:

    % qmgr -c 'p s'

計算ノードの認識を確認:

    % pbsnodes -a

投入したジョブの状態を確認:

    % qstat
    % qstat -q

ログを見る:

    % less /var/spool/torque/server_logs/2012xxxx
    % less /var/spool/torque/sched_logs/2012xxxx

サービス再起動:

    % sudo service torque-server restart
    % sudo service torque-sched restart

計算ノードが落ちたあと `qstat` リストに残ってしまうジョブを消す:

    % sudo qdel -p

### Compute node

状態確認:

    % sudo momctl -d3

ログを見る:

    % less /var/spool/torque/mom_logs/2012xxxx

サービス再起動:

    % sudo service torque-mom restart
