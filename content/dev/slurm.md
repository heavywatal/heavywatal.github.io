+++
title = "Slurm"
tags = ["job", "concurrent"]
aliases = ["torque.html"]
[menu.main]
  parent = "dev"
+++

Migrated from [TORQUE](https://wiki.archlinux.org/title/TORQUE) and [OpenPBS](https://github.com/openpbs/openpbs).

- Ubuntu 24.04 Noble Numbat
- [Slurm 23.11.4](https://launchpad.net/ubuntu/noble/+package/slurmd)


## Installation

- <https://slurm.schedmd.com/quickstart_admin.html>
- <https://memo.fifthdimension.jp/configure-slurm-workload-manager-on-ubuntu>

Common for both controler and compute nodes:
```sh
sudo apt updata && apt list --upgradable
sudo apt upgrade
sudo apt install munge slurmd
```

Set NTP server and time zone:
```sh
sudo timedatectl set-timezone Asia/Tokyo
sudo vim /etc/systemd/timesyncd.conf
```
```
[Time]
NTP=YOUR_NTP_SERVER_HERE
```
```sh
sudo systemctl restart systemd-timesyncd
systemctl status systemd-timesyncd
timedatectl timesync-status
```

### Controller node

```sh
sudo apt install slurm-wlm slurmctld
```

Create `/etc/munge/munge.key`:
```sh
sudo -u munge /usr/sbin/mungekey -v
```

Create [`/etc/slurm/cgroup.conf`](https://slurm.schedmd.com/cgroup.conf.html).
It can be empty, but you may want to try resource constraints:
```text
ConstrainCores=yes
ConstrainDevices=yes
ConstrainRAMSpace=yes
```

Use `/usr/share/doc/slurmctld/slurm-wlm-configurator.easy.html` to generate
an almost-default [`/etc/slurm/slurm.conf`](https://slurm.schedmd.com/slurm.conf.html).

- `ClusterName`
- `SlurmctldHost`
- `NodeName`
- `PartitionName`
- `CPUs`: Leave it blank to get `Sockets * CoresPerSocket * ThreadsPerCore`.
- `Sockets`
- `CoresPerSocket`
- `ThreadsPerCore`
- `RealMemory`: Total memory in MB minus some margin for OS.

Prepare directories for log files and state files:
```sh
sudo mkdir -p /var/lib/slurm/slurmctld
sudo mkdir -p /var/lib/slurm/slurmd
sudo mkdir -p /var/log/slurm
sudo chown -R slurm:slurm /var/lib/slurm
sudo chown -R slurm:slurm /var/log/slurm
```

Set permissions for configuration files:
```sh
sudo chown -R slurm:slurm /etc/slurm
sudo chmod 644 /etc/slurm/slurm.conf
sudo chmod 644 /etc/slurm/cgroup.conf
```

Enable and start services:
```sh
sudo systemctl enable --now munge slurmd slurmctld
systemctl status munge slurmd slurmctld
```

Copy `munge.key`, `slurm.conf`, and `cgroup.conf` to some NFS-shared directory, e.g., `/home/admin/etc/slurm/`, to distribute them to compute nodes later.


### Compute node

Prepare directories for log files and state files:
```sh
sudo mkdir -p /var/lib/slurm/slurmd
sudo mkdir -p /var/log/slurm
sudo chown -R slurm:slurm /var/lib/slurm
sudo chown -R slurm:slurm /var/log/slurm
```

Copy `munge.key`, `slurm.conf`, and `cgroup.conf` from the controler node:
```sh
cd /home/admin/etc/slurm/

sudo cp munge.key /etc/munge/
sudo chmod 400 /etc/munge/munge.key
sudo chown munge:munge /etc/munge/munge.key

sudo cp slurm.conf /etc/slurm/
sudo cp cgroup.conf /etc/slurm/
sudo chmod 644 /etc/slurm/*.conf
```

Enable and start services:
```sh
sudo systemctl enable --now munge slurmd
systemctl status munge slurmd
```


### Check

```sh
sinfo
scontrol show nodes
systemctl status munge slurmd slurmctld
sudo less /var/log/slurm/slurmctld.log
sudo less /var/log/slurm/slurmd.log
```

You may need to update the state manually:
```sh
sudo scontrol reconfigure
sudo scontrol update nodename=compute01 state=resume
```
