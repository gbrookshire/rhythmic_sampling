"""
Rhythmic sampling experiment.
- Project name: 'Do brain waves help you to perceive?'
- Project code: 2018-289
- Grant code: GMGM.ROKZ19867
Fall 2018
G.Brookshire@bham.ac.uk

v2
One (or zero) target per trial
Target ends when participant responds or at 4 s, whichever comes sooner

Distances
Eyes to screen: 147.5 cm
Screen width: 70.3 cm (25.5 deg)
Screen height: 39.5 cm (14.3 deg)

Functions to compute visual angle
pix2cm = lambda x: x * 70.3 / 960
pix2deg = lambda x: np.rad2deg(np.arctan(pix2cm(x) / 147.5))
"""

from __future__ import division

# Change to False for testing on office desktop
IN_MEG_LAB = True

import os
import random
import datetime
import numpy as np
from psychopy import visual, core, data, event, monitors

import refcheck # To check the screen refresh rate
import flicker as fl # Flickering stims on the Propixx projector
import paddle # MEG-compatible LED response box over LabJack
import eye_wrapper # Wrapper around pylink to interface with eye tracker

if IN_MEG_LAB:
    from psychopy import parallel
else:
    import dummy_parallel as parallel

############
# Settings #
############

# Clocks
START_TIME = datetime.datetime.now().strftime('%Y-%m-%d-%H%M')
RT_CLOCK = core.Clock() # response times

port = parallel.ParallelPort(address=0xBFF8)
port.setData(0) # Reset the parallel port

LOG_DIR = '../logfiles/'
STIM_DIR = '../stimuli/'
assert os.path.exists(LOG_DIR)
assert os.path.exists(STIM_DIR)

# Load instructions
with open('lib/instruct.txt') as f:
    instruct_text = f.readlines()

# Parallel port triggers - Send triggers on individual channels
triggers = {'trial begin': 3,
            'target': 4,
            'response': 5}

flicker_freqs = [63., 78.]

stim_size = 2 ** 8 # Diameter of faces in pixels
eccent = int(stim_size * 6/10) # Stimulus eccentricity
frame_center = (fl.FRAME_CENTER[0], fl.FRAME_CENTER[1] + 80) # Shift stimuli up
stim_loc = [(frame_center[0]-eccent, frame_center[1]-eccent),
            (frame_center[0]+eccent, frame_center[1]-eccent)]
stim_orientation = [45, -45] # faces are rotated to point toward fixation

trial_length = 4 # Seconds
target_duration = 0.033
target_interval = [1.5, trial_length - 1]
early_target_interval = [0.5, 1.5]
fix_dur = [1.0, 2.0] # Select from a uniform distribution
iti = 1.0 # Blank screen between trials

n_targets_per_block = 48
n_early_targ_per_block = 8
n_no_targ_per_block = 0

n_blocks_thresh = 4 # N blocks in the thresholding portion of exp
n_blocks_main = 6 # N trials in the main portion of the exp

n_trials = n_targets_per_block + n_early_targ_per_block + n_no_targ_per_block
n_trials_thresh = n_trials * n_blocks_thresh
n_trials_main = n_trials * n_blocks_main

CS = 'rgb'  # ColorSpace
WHITE = [1., 1., 1.]
GREY = [0., 0., 0.] # the background of face pictures
BLACK = [-1., -1., -1.]

# Initialize and calibrate eye-tracker #
if IN_MEG_LAB:
    fl.close() # Make sure projector is in normal mode before starting
    el = eye_wrapper.SimpleEyelink(fl.FULL_RES)
    el.startup()

# Initialize the Propixx projector
fl.init(use_propixx=IN_MEG_LAB)

# Initialize paddle buttonbox
paddle.init_labjack(dummy_run=(not IN_MEG_LAB))
buttonbox = paddle.Paddle(register=paddle.LEFT_BUTTON,
                          direction=paddle.DOWN,
                          clock=RT_CLOCK)

# Keep track of responses within each
RESP_COUNTER = {'hit':0, 'miss':0, 'fa':0}


######################
# Window and Stimuli #
######################

win = visual.Window(fl.FULL_RES, monitor='testMonitor',
                    fullscr=IN_MEG_LAB,
                    color=GREY, colorSpace=CS,
                    allowGUI=False, units='pix')

# Common parameters used across stimuli
stim_params = {'win': win, 'units': 'pix'}

text_stim = fl.QuadStim(visual.TextStim, pos=frame_center, # For instructions
                           color=WHITE, colorSpace=CS,
                           height=30, **stim_params)

fixation = fl.QuadStim(visual.Circle, radius=3,
                          fillColor=WHITE, lineColor=WHITE,
                          fillColorSpace=CS, lineColorSpace=CS,
                          pos=frame_center, **stim_params)

# Patch for measuring temporal characteristics with photo diode
photodiode_size = 40
photodiode_pos = (frame_center[0] - frame_center[0] - photodiode_size/2,
                  frame_center[1] - frame_center[1] + photodiode_size/2)
photodiode_patch = fl.BrightnessFlickerStim(visual.ImageStim,
                                  image=np.array([[1.,1.],[1.,1.]]),
                                  colorSpace=CS,
                                  pos=photodiode_pos,
                                  size=(photodiode_size, photodiode_size),
                                  **stim_params)

# Target stimuli
face_angular_width = np.arctan(1/2)
face_angular_center = -3 * np.pi / 4
face_angles = [face_angular_center + face_angular_width,
               face_angular_center - face_angular_width]
targ_ang = np.linspace(*face_angles, num=12)
# Exclude targets overlapping the nose
targ_ang = np.delete(targ_ang,
                     np.where(np.abs(targ_ang + (3 * np.pi / 4)) < 0.15))
n_targ_locs = len(targ_ang)
# Get target locations for the other side of the screen
targ_ang = np.concatenate([targ_ang, targ_ang + np.pi / 2])
# Exclude angles that overlap with the eyes or mouth
targ_loc = (stim_size * 0.95) * np.exp(1j * targ_ang) # Location in polar coords
targ_loc = zip(np.real(targ_loc), np.imag(targ_loc)) # Polar to cartesian
targets = [] # List of target stimuli
n_stim = -1 # Keep track of target stimuli to exclude overlap with eyes/mouth
for targ_x, targ_y in targ_loc:
    n_stim +=1
    pos = (targ_x + frame_center[0],
           targ_y + frame_center[1])
    s = fl.QuadStim(visual.GratingStim,
                       tex=None,
                       size=(60, 60),
                       mask='gauss',
                       maskParams={'sd': 5},
                       color=BLACK,
                       pos=pos, **stim_params)
    targets.append(s)

# Build the list of test stimuli
images = ['{}/{}.jpg'.format(STIM_DIR, n) for n in range(1,7)]
stimuli = []
for fname in images:
    s = fl.BrightnessFlickerStim(visual.ImageStim, image=fname,
                        colorSpace=CS,
                        pos=frame_center,
                        size=(stim_size, stim_size),
                        mask='circle',
                        **stim_params)
    stimuli.append(s)


#############################
# Build the trial structure #
#############################

cond_combo = [{'target_side': side, 'target_side_freq':freq}
                for side in ('left', 'right') for freq in flicker_freqs]

assert n_targets_per_block % len(cond_combo) == 0
assert n_no_targ_per_block % len(cond_combo) == 0
assert n_early_targ_per_block % len(cond_combo) == 0

def freq_sides(target_side, target_side_freq):
    """ Get the frequency of flickers at each side
    given the side target is on, and the flicker freq of that side
    """
    if target_side == 'left':
        if target_side_freq == flicker_freqs[0]:
            return flicker_freqs
        else:
            return flicker_freqs[::-1]
    elif target_side == 'right':
        if target_side_freq == flicker_freqs[0]:
            return flicker_freqs[::-1]
        else:
            return flicker_freqs

trial_info = []
for i_block in range(n_blocks_thresh + n_blocks_main):
    # Put together the main trials
    # Get all the counterbalanced conditions
    if i_block < n_blocks_thresh:
        quest_thresholding = True
    else:
        quest_thresholding = False
    conds = cond_combo * int(n_trials / len(cond_combo))
    block_trials = []
    for i_trial in range(n_trials):
        d = {}
        d['quest'] = quest_thresholding
        d['fix_dur'] = random.uniform(*fix_dur)
        d['stim_left'], d['stim_right'] = random.sample(range(6), 2)
        # Assign flicker freq to each side
        d['target_side'] = conds[i_trial]['target_side']
        d['target_side_freq'] = conds[i_trial]['target_side_freq']
        d['freq_left'], d['freq_right'] = freq_sides(d['target_side'],
                                                     d['target_side_freq'])
        # Main target trials
        if i_trial < n_targets_per_block:
            d['target_t'] = random.uniform(*target_interval)
        # Catch trials with early targets
        elif i_trial < n_targets_per_block + n_early_targ_per_block:
            d['target_t'] = random.uniform(*early_target_interval)
        # Catch trials with no target
        else:
            d['target_t'] = None
        # Set up the target position
        if d['target_t'] is None:
            d['target_pos'] = None
        else:
            d['target_pos'] = random.choice(xrange(n_targ_locs))
        # If target is on the right  side, make sure it gets the correct inx
        # Right side targets are in the second half of the list
        if d['target_t'] is not None and d['target_side'] == 'right':
            d['target_pos'] += n_targ_locs

        block_trials.append(d)

    random.shuffle(block_trials)
    trial_info.extend(block_trials)

trials = data.TrialHandler(trial_info, nReps=1, method='sequential')

# QUEST Adaptive procedure for titrating target opacity
guess_prob = 0.01 # P(trials) on which subject just guesses blindly
def quest():
    """ Convenience function to help make multiple QUEST staircases
    """
    q = data.QuestHandler(startVal=0.25, startValSd=0.225,
                          pThreshold=0.5,
                          method='mean', # quantile, mean, or mode
                          beta=3.5, # Steepness of psychometric func
                          delta=guess_prob, # P(trials) S presses blindly
                          gamma=guess_prob/2, # P(corr) when intensity=-inf
                          minVal=0, maxVal=1, range=1)
    return q

staircases = {side: {freq: quest() for freq in flicker_freqs}
                for side in ('right', 'left')}


###################################
# Functions to run the experiment #
###################################

def send_trigger(trig):
    """ Send triggers to the MEG acquisition computer and the EyeLink computer.
    """
    port.setPin(trig, 1)
    el.trigger(trig)

def s2f(x):
    """ Convert seconds to frame flips. """
    return int(x * fl.OUTPUT_FRAME_RATE)

def f2s(x):
    """ Convert frame flips to seconds. """
    return x / fl.OUTPUT_FRAME_RATE

def check_accuracy(t_target, t_response, threshold):
    if t_target is None: # No target was shown
        hit = None
        false_alarm = bool(t_response)
    else: # There was a target shown
        if t_response is None:
            hit = False
            false_alarm = False
        elif t_response < t_target:
            hit = None
            false_alarm = True
        elif t_response < t_target + threshold:
            hit = True
            false_alarm = False
        elif t_response >= t_target + threshold:
            hit = False
            false_alarm = True
        else:
            raise Exception('Failed to classify hits and false alarms')
    return hit, false_alarm

def show_text(text):
    """ Show text at the center of the screen.
    """
    text_stim.set('text', text)
    text_stim.draw()
    win.flip()

def instructions(text):
    """ Show instructions and go on after pressing space
    """
    show_text(text)
    event.waitKeys(keyList=['space'])
    win.flip(clearBuffer=True) # clear the screen
    core.wait(0.2)

def example_screen():
    """ Show an example of what a stimulus display will look like.
    """
    fixation.draw()
    for n_stim in range(2):
        stimuli[n_stim].set_pos(stim_loc[n_stim])
        stimuli[n_stim].set('ori', stim_orientation[n_stim])
        stimuli[n_stim].draw()
    instructions('The images will look like this:\n \n \n \n \n')
    fixation.draw()
    for n_stim in range(2):
        stimuli[n_stim].set_pos(stim_loc[n_stim])
        stimuli[n_stim].set('ori', stim_orientation[n_stim])
        stimuli[n_stim].draw()
    targets[4].draw()
    instructions('The target will look like this:\n \n \n \n \n')
    instructions('Any questions before beginning?')

def take_a_break():
    global N_BLOCKS_REMAINING
    global RESP_COUNTER
    N_BLOCKS_REMAINING -= 1
    msg = "Take a break\nBlocks remaining: {}\n\nH: {}, M: {}, F: {}"
    if (N_BLOCKS_REMAINING % 2) == 0:
        msg += '\n\n(Start next recording)'
    show_text(msg.format(N_BLOCKS_REMAINING,
                         RESP_COUNTER['hit'],
                         RESP_COUNTER['miss'],
                         RESP_COUNTER['fa']))
    RESP_COUNTER = {'hit':0, 'miss':0, 'fa':0} # Reset response counter
    # Press ESCAPE to quit early
    r = event.waitKeys(keyList=['space', 'escape'])
    if 'escape' in r:
        show_text('Press Q to exit the study, or SPACE to go on.')
        r = event.waitKeys(keyList=['q', 'space'])
        if 'q' in r:
            return 'quit'

def show_trial(trial):
    """ Present a trial, given a dictionary of trial attributes
    """
    port.setData(0) # set all pins low

    # Brief blank screen before each trial
    win.flip(clearBuffer=True)
    core.wait(iti)

    # Variable fixation before each trial
    fixation.draw()
    win.flip()
    core.wait(trial['fix_dur'])

    # Set the display parameters of the stimuli
    for n_side, side in enumerate(['left', 'right']):
        stimuli[trial['stim_'+side]].flicker(trial['freq_'+side])
        stimuli[trial['stim_'+side]].set_pos(stim_loc[n_side])
        stimuli[trial['stim_'+side]].set('ori', stim_orientation[n_side])

    # Make sure the photodiode patch is set to the same as the right side stim
    photodiode_patch.flicker(trial['freq_right'])

    # Get the frame numbers on which targets should be drawn
    if trial['target_t'] is None:
        target_on_frames = []
    else:
        target_on_frames = range(s2f(trial['target_t']),
                                 s2f(trial['target_t'] + target_duration))

    # Show the stimuli
    inter_frame_int = [] # Keep track of frame timing
    target_on = False
    buttonbox.clear_events()
    event.clearEvents()
    RT_CLOCK.reset()
    for n_frame in xrange(s2f(trial_length)):
        # Should the target be displayed during this frame?
        if n_frame in target_on_frames:
            target_on = True
        else:
            if target_on:
                target_on = False

        # Draw all the stimuli
        fixation.draw()
        photodiode_patch.draw()
        stimuli[trial['stim_left']].draw()
        stimuli[trial['stim_right']].draw()
        if target_on:
            targets[trial['target_pos']].draw()
        win.flip()

        # Send triggers right after the window flips
        if target_on: # Trigger if the target is being displayed
            send_trigger(triggers['target'])
        elif n_frame == 0: # Trigger for beginning of the trial
            send_trigger(triggers['trial begin'])
        else:
            port.setData(0)

        # End the trial if there are any responses
        buttonbox.check_for_events()
        if buttonbox.responses:
            send_trigger(triggers['response'])
            break
        inter_frame_int.append(RT_CLOCK.getTime())

    # delta = np.diff(inter_frame_int)
    # print 'Frame period:', np.mean(delta), '+/-', np.std(delta)

    # Clear the screen
    win.flip(clearBuffer=True)
    core.wait(0.1)

    # Check for responses
    rt = buttonbox.last_resp()
    hit, false_alarm = check_accuracy(trial['target_t'], rt, threshold=1.0)

    # print 'targ: {}, rt: {}, hit: {}, fa: {}'.format(
    #     trial['target_t'], rt, hit, false_alarm)

    # Keep track of responses for feedback between blocks
    if hit == True:
        RESP_COUNTER['hit'] += 1
    elif hit == False:
        RESP_COUNTER['miss'] += 1
    elif false_alarm == True:
        RESP_COUNTER['fa'] += 1

    return rt, hit, false_alarm

def set_target_opacity(trial):
    """ Set the opacity of the targets for this trial
    """
    side = trial['target_side']
    freq = trial['target_side_freq']
    if trial['quest'] and \
        trial['target_t'] and \
        trial['target_t'] > target_interval[0]:
            opacity = staircases[side][freq].next()
    else:
        opacity = staircases[side][freq].mean()
    for targ in targets:
        targ.set('opacity', opacity)
    return opacity

def update_staircase(trial, hit):
    """ Update the correct QUEST staircase given the trial dict, and whether
    this trial was a Hit
    """
    side = trial['target_side']
    freq = trial['freq_' + side]
    if trial['target_t'] and trial['target_t'] > target_interval[0]:
        if hit is not None:
            staircases[side][freq].addResponse(int(hit)) # Hit=1, Miss=0

def data_collection_phase():
    """ Main data collection segment of the experiment
    """
    global N_BLOCKS_REMAINING # For showing how many blocks of trials remain
    N_BLOCKS_REMAINING = n_blocks_main + n_blocks_thresh

    # Show the trials
    for n, trial in enumerate(trials):

        # Take a break between blocks
        if (n > 0) and (n % n_trials == 0):
            if take_a_break() == 'quit':
                break

        opacity = set_target_opacity(trial)
        rt, hit, false_alarm = show_trial(trial)
        trials.addData('rt', rt)
        trials.addData('hit', hit)
        trials.addData('false_alarm', false_alarm)
        trials.addData('target_opacity', opacity)

        if trial['quest']:
            update_staircase(trial, hit)

    # Save the data
    trials.saveAsWideText('{}/{}.csv'.format(LOG_DIR, START_TIME),
                          delim=',',
                          fileCollisionMethod='rename')

    # Save the thresholds and save the staircases as a pickled object
    for side in ('left', 'right'):
        for freq in flicker_freqs:
            staircases[side][freq].saved_thresh = staircases[side][freq].mean()
            fname = '{}/{}_{}_{}.csv'.format(LOG_DIR, START_TIME, side, round(freq))
            # Important objects in pickled file are d.data and d.intensities
            staircases[side][freq].saveAsPickle(fname[:-4])
            staircases[side][freq].saveAsJson(fname[:-4]+'.json')

def stim_test():
    """ Test the stimuli and target positions on the screen
    """
    # Show all targets
    fixation.draw()
    for n in [0, 1]:
        stimuli[n].set_pos(stim_loc[n])
        stimuli[n].set('ori', stim_orientation[n])
        stimuli[n].draw()
    for targ in targets:
        targ.draw()
    win.flip()
    event.waitKeys(keyList=['space'])

    # Show just one target
    fixation.draw()
    for n in [0, 1]:
        stimuli[n].draw()
    targets[5].draw()
    win.flip()
    event.waitKeys(keyList=['space'])

def run_exp():
    """ Run the different sections of the experiment
    """
    refcheck.check_refresh_rate(win, fl.OUTPUT_FRAME_RATE)

    # Show the instructions
    for text in instruct_text:
        instructions(text)
    example_screen()

    # Main experiment
    instructions('Ready to begin the study?')
    core.wait(1)
    data_collection_phase()

    show_text('That was it -- thanks!')
    event.waitKeys(keyList=['escape'])

    # Close everything down
    win.close()
    fl.close()
    paddle.close_labjack()
    if IN_MEG_LAB:
        el.shutdown()
    core.quit()
