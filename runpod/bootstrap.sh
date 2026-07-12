#!/bin/bash
# Pod bootstrap: sshd for result retrieval + kick off the QE batch.
apt-get update -qq && apt-get install -y -qq openssh-server git rsync >/dev/null 2>&1
mkdir -p /run/sshd /root/.ssh
echo 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIIWsQNj8wVSSCkStELR4hMyiP8xm49qjTopL+IdmG9Yw runpod-access' >> /root/.ssh/authorized_keys
echo 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAICARC8tdRQcP2eIR0zcjnyMEmLO4PxMQF8k+mE0f0RY7 santosh@agnostiq.ai' >> /root/.ssh/authorized_keys
chmod 700 /root/.ssh && chmod 600 /root/.ssh/authorized_keys
/usr/sbin/sshd
git clone https://github.com/santoshkumarradha/NaxCoO2-Superconductivity /workspace/repo
cd /workspace/repo/runpod && nohup bash run_all.sh > run_all.log 2>&1 &
sleep infinity
