Linux troubleshooting credit from [Guide](https://support.arduino.cc/hc/en-us/articles/360016495679-Fix-port-access-on-Linux):

When programming fails of matlab setup step, fix user permissions:

To see the group of permissions of dev:
```bash
ls -l /dev/ttyACM*
```

Usually is `dialout` so:

```bash
sudo usermod -a -G dialout $USER
`` 