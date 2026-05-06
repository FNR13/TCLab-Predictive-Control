Linux troubleshooting credit from [Guide](https://support.arduino.cc/hc/en-us/articles/360016495679-Fix-port-access-on-Linux):

When programming fails of matlab setup step, fix user permissions:

To see the group of permissions of dev:
```bash
ls -l /dev/ttyACM*
```

Usually is `dialout` so:

```bash
sudo usermod -a -G dialout $USER
```

Parameter | Value
Transistors (heaters) | BJT TIP31C in TO-220 package
Maximum | heater power10 W
PWM | discretization levels28
Thermistors (sensors) | TMP36GZ
Operating range | −40◦ C to 150◦ C
Sensor accuracy at room temperature | (25◦ C)±1◦ C
General sensor accuracy | ±2◦ C
Heater shut-off temperature | 100◦ C