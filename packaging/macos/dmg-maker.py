#! /usr/bin/env python
"""
Creates a nice disk image, with background and /Applications symlink
for the app.

Usage: dmg_maker.py path/to/Regina.app

Thanks to Nathan Dunfield for this original version of this script.
It has been modified since, and so any bugs are most likely mine. - B.B.

One issue here is that Snow Leopard uses a different (undocumented, of
course) format for the .DS_Store files than earlier versions, which makes
disk images created on it not work correctly on those systems.   Thus this "solution" uses a .DS_Store file created on Leopard as follows:

(1) Use Disk Utility to create a r/w DMG large enough to store everything and open it.

(2) Copy over the application and add a symlink to /Applications.

(3) Create a subdirectory ".background" containing the file "background.png".

(4) Open the disk image in the Finder and do View->Hide Tool Bar and then View->Show View Options.  To add the background picture inside the hidden directory, use cmd-shift-g in the file dialog.  Adjust everything to suit, close window and open it.   Then copy the .DS_Store file to dotDS_store.  

"""
import os, sys, re, commands
from math import ceil

name = "Regina"
version = "4.93"

allow64 = commands.getstatusoutput('sysctl -n hw.optional.x86_64')[1]
if allow64 == '1':
    arch = 'x86_64'
else:
    arch = 'i386'

kernel = commands.getstatusoutput('uname -r | cut -d. -f1')[1]
if kernel == '9':
    osver = 'Leopard'
elif kernel == '10':
    osver = 'SnowLeopard'
elif kernel == '11':
    osver = 'Lion'
else:
    print 'Unknown MacOS kernel version:', kernel
    sys.exit(1)

dmg_real = name + "-" + version + "_" + osver + "-" + arch + ".dmg";
dmg_tmp = name + "-" + version + "-tmp.dmg";
dist_dir = "dist"

def main():
    # Make sure we are running from the right location.
    if not (os.path.exists('dotDS_Store') and os.path.exists('background.png')):
        print "Please run this script from the directory that contains it."
        sys.exit(1)

    # Make sure the argument is valid.
    if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " path/to/Regina.app"
        sys.exit(1)
    if not (os.path.exists(sys.argv[1] + "/Contents/MacOS/Regina")):
        print "The argument " + sys.argv[1] + " looks invalid."
        print "Usage: " + sys.argv[0] + " path/to/Regina.app"
        sys.exit(1)

    # Strip any trailing slash on the argument, which messes up the
    # recursive copy.
    app = sys.argv[1]
    if app[-1:] == '/':
        app = app[:-1]

    # If there is already a sandbox, be a coward and back out.
    if os.path.exists(dist_dir):
        print "There is already a sandbox at " + dist_dir + "."
        print "Please either delete it or move it out of the way."
        sys.exit(1)

    # Make sure the dmg isn't currently mounted, or this won't work.  
    mount_name = "/Volumes/" + name
    while os.path.exists(mount_name):
        print("Trying to eject " + mount_name)
        os.system('hdiutil detach "%s"' % mount_name)

    # Remove old dmg if there is one
    while os.path.exists(dmg_real):
        os.remove(dmg_real)
    while os.path.exists(dmg_tmp):
        os.remove(dmg_tmp)

    # Build a new sandbox.
    os.mkdir(dist_dir)
    os.mkdir(dist_dir + '/.background')
    os.symlink("/Applications", dist_dir + "/Applications")
    os.system('cp background.png "%s/.background"' % dist_dir)
    os.system('cp dotDS_Store "%s/.DS_Store"' % dist_dir)
    os.system('cp -pr "%s" "%s"' % (app, dist_dir))
    
    # figure out the needed size:
    raw_size = os.popen("du -sh " + dist_dir).read()
    size, units = re.search("([0-9.]+)([KMG])", raw_size).groups()
    new_size = "%d" % ceil(1.2 * float(size)) + units

    # Run the main script:
    print "Building an initial DMG..."
    os.system('hdiutil makehybrid -hfs -hfs-volume-name "%s" -hfs-openfolder "%s" "%s" -o "%s"' % (name, dist_dir, dist_dir, dmg_tmp))
    print "Converting the DMG format..."
    dmg_format = 'UDZO' # 'UDRW' if you want to edit the .DS_Store
    os.system('hdiutil convert -format %s "%s" -o "%s"' % (dmg_format, dmg_tmp, dmg_real))
    print "Cleaning up..."
    os.remove(dmg_tmp)
    print "Done."
              
    
    
if __name__ == "__main__":
    main()



