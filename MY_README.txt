to use as a non-root user follow this procedure:
- connect drs4 board to PC
- open terminal
- lsusb (to list devices connect)


- sudo nano /etc/udev/rules.d/60-drs4board.rules


Put SUBSYSTEM=="usb", ATTRS{idVendor}=="VID", ATTRS{idProduct}=="PID", MODE="0666"
example: SUBSYSTEM=="usb", ATTRS{idVendor}=="04b4", ATTRS{idProduct}=="1175", MODE="0666"


    VID is the USB-IF-assigned Vendor ID of the device in question *
    PID is the Vendor-assigned Product ID of the device in question *

    0666 gives universal read/write access to whatever matches this line

    * $ lsusb to see all attached USB devices and their ID's.

In /etc/udev/rules.d/xx-my-rule.rules (may need root/sudo permissions)

    xx is any number > 50 (the defaults are in 50, and higher numbers take priority)
    my-rule is whatever you want to call it
    must end in .rules

Then udevadm control --reload-rules (may also need root/sudo permissions), and it should "just work" for that specific VID/PID pair.
