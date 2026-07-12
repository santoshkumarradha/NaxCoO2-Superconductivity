#!/bin/bash
# Pod bootstrap: sshd for result retrieval + kick off the QE batch.
apt-get update -qq && apt-get install -y -qq openssh-server git rsync python3-numpy >/dev/null 2>&1
mkdir -p /run/sshd /root/.ssh
echo 'ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQDSfSnxVADxon3kKkE8ah+zowVaEArPF4SA5MmkRrZw+TRjXzAqWrpYcPFJHxhU9+f5f4YXjteHIRbAIpI+xKcJ212jfKypOmx8eQWj3yfzR7cMpyzDeO5zj8zwyi4NtLTZy61odXDU5ip3e9Fr+hlv6UaLMDURYY5KSbPYvyS2tlg6BCGQyi6VziJiddiWxdHlBT3+74zpq7YwGW6DXN0o3Q3HY2ECgULrSDJeAxs7HgZQ3h4nZPNwZuXWrBPrN/l3YXzf782Kz8WFqHz7WTZecRiJUw/zmKUFfBnX0hGCBZb8gqf232opI0X5h38M/Tjap86MWxqIEsiXaFqkWBKh runpodctl-ssh-key' >> /root/.ssh/authorized_keys
[ -n "$PUBLIC_KEY" ] && echo "$PUBLIC_KEY" >> /root/.ssh/authorized_keys
usermod -p '*' root
echo 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIIWsQNj8wVSSCkStELR4hMyiP8xm49qjTopL+IdmG9Yw runpod-access' >> /root/.ssh/authorized_keys
echo 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAICARC8tdRQcP2eIR0zcjnyMEmLO4PxMQF8k+mE0f0RY7 santosh@agnostiq.ai' >> /root/.ssh/authorized_keys
chmod 700 /root/.ssh && chmod 600 /root/.ssh/authorized_keys
/usr/sbin/sshd
git clone https://github.com/santoshkumarradha/NaxCoO2-Superconductivity /workspace/repo
cd /workspace/repo/runpod && nohup bash run_bands.sh > run_bands.log 2>&1 &
sleep infinity
