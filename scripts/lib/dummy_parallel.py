"""
Simple class that does nothing, to test PsychoPy scripts on
a computer without parallel ports.
"""

class ParallelPort:
    print_triggers = False  # Set this to True to print the trigger
    def __init__(self, address):
        pass
    def setData(self, x):
        if self.print_triggers:
            print "parallel:", x
    def setPin(self, *args):
        if self.print_triggers:
            print "parallel:", args
