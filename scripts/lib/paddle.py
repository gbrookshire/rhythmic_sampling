"""
Use the Electa-provided MEG-compatible LED paddle buttons in psychopy.
For an example of usage, see the test function at the end of the file.
"""
# Might be able to use psychopy.iohub.devices.daq

try:
    from psychopy.hardware.labjacks import U3
except ImportError:
    print 'Labjack not found'
    U3 = 'NOT FOUND'

from psychopy import event

# Registers of the buttons
LEFT_BUTTON = 16  # Works better
RIGHT_BUTTON = 17

# Constants to show button position
UP = 0
DOWN = 1

# LabJack instance (global to allow access to multiple buttons)
LJ = None

def init_labjack(dummy_run=False):
    global LJ
    if not dummy_run:
        if isinstance(LJ, U3):
            print 'LabJack is already initialized. Doing nothing.'
        LJ = U3()
    else:
        print 'Dry run: Press the space bar instead of the LabJack button.'
        LJ = DummyU3()

def close_labjack():
    if LJ is not None:
        LJ.close()
    else:
        print 'LabJack was not initialized'

class DummyU3(object):
    """ For testing on a computer that doesn't have a LabJack device.
        Check whether the space bar is pressed
    """
    def __init__(self):
        event.clearEvents()
    def getDIOState(self, register):
        resp = event.getKeys(['space'])
        if resp:
            r = 1
        else:
            r = 0
        event.clearEvents()
        return r

class Paddle(object):
    """ Class to collect responses on the Electra button box
        connected w/ LabJack.

        Arguments
        register: A LabJack register
        direction: To collect key presses choose paddle.DOWN, for releases
            choose paddle.UP
        clock: An instance of psychopy.core.Clock for computing RTs.
            Note that temporal resolution will only be as good as the rate
            at which you call check_for_events()
    """

    def __init__(self, register, clock, direction=DOWN):
        self.register = register
        self.direction=direction
        self.clock = clock
        self.clear_events()

    def clear_events(self):
        self.prev_state = None
        self.responses = []

    def poll(self):
        return LJ.getDIOState(self.register)

    def check_for_events(self):
        curr_state = self.poll()
        if (curr_state != self.prev_state) and (curr_state == self.direction):
            rt = self.clock.getTime()
            self.responses.append(rt)
        self.prev_state = curr_state

    def last_resp(self):
        if self.responses == []:
            return None
        else:
            return self.responses[-1]


if __name__ == '__main__':
    # Run a test of the paddle
    from psychopy import core, visual
    rt_clock = core.Clock()
    init_labjack(dummy_run=False)
    p = Paddle(register=LEFT_BUTTON,
               direction=UP,
               clock=rt_clock)
    win = visual.Window(size=(600, 600))
    for n in range(3):
        rt_clock.reset()
        p.clear_events()
        print 'Trial %i: Press the button a couple times' % (n + 1)
        while rt_clock.getTime() < 5:
            p.check_for_events()
            core.wait(0.1)
        print '    RTs:', p.responses
    close_labjack()
