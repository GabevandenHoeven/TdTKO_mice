import sys
import time


def animate():
    sys.stdout.write('Loading ')
    for _ in range(4):
        time.sleep(.5)
        sys.stdout.write('.')
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write('\r')
