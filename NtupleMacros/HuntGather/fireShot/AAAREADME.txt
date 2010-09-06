# Create initial directory structure
mkdir picks shots xwd

# Make frame buffer on :1 and always
# use it
export DISPLAY=:1
Xvfb $DISPLAY -screen 0 1600x1200x24 -fbdir /store/disk00/jribnik/fireShot &

# Substitute with whatever cmsShow v
# is supposed to work
./cmsShow35c-2010-03-31/cmsShow -c xvfb_1600x1200_cmsShow35c.fwc --play 5 --port 9999 ./cmsShow35c-2010-03-31/data.root &

# This script runs forever and feeds
# the cmsShow, takes screenshots and
# and sends them to the uaf
./fireShot
