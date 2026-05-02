Linux troubleshooting:

To see the group of permissions of dev:
```bash
ls -l /dev/ttyACM*
```

Usually is `dialout` so:

```bash
sudo usermod -a -G dialout $USER
`` 