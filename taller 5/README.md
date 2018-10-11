# MPICH2 in Ubuntu 14
## Install
Sure Updated last apt list
```bash
$ sudo apt-get update -Force::IPv4=true
```
verify system is alredy up to date (upgraded), remember if want to upgrade but dont change version ubuntu use dist-upgrade
```bash
$ sudo apt-get dist-upgrade -o Force::IPv4=true
```

now can install mpi, in ubuntu 14 is mpich2 (in others is mpich)
```bash
$ sudo apt-get install mpich2
```

this automatically install all dependencies (build eseencials and others) and set in enviroment variables mpicc, mpiCC, mpicxx, mpic++. This doit in all machines to work (Header and workers)

**Note:** All devices should have the same username and installed on this user

## Configure comunication
for correct execution in workers machines from Header machine, need to enable comunication first on all machines setup name machine (for simplicity Header, Worker1, Worker2... Workern) in ubuntu the name of device is configured in */etc/hostname*

```bash
$ sudo nano /etc/hostname
```

On this file is only the name of device.
### Header Device
```bash
Header
```
### Worker1 Device
```bash
Worker1
```

if name of device change the "localhost" modifiy config file in */etc/hosts*
### Header Device
```bash
127.0.0.1   localhost
172.0.1.1   Header
```
### Worker Device
```bash
127.0.0.1   localhost
127.0.1.1   Worker1
```

## compiling