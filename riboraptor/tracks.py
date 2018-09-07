import glob
import os
from textwrap import indent, dedent
import textwrap
import sys
import trackhub

INDENT_SPACES = '\t'

def create_track(hub_name,
                 short_label,
                 long_label,
                 genome,
                 glob_pattern,
                 email='skchoudh@usc.edu',
                 host='localhost',
                 remote_dir='/staging/as/skchoudh/riboraptor_trackhub'):
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=hub_name,
        short_label=short_label,
        long_label=long_label,
        genome=genome,
        email=email)
    for bigwig in glob.glob(glob_pattern):
        sample_name = os.path.dirname(bigwig).split(os.path.sep)[-2]
        fragment_length = os.path.dirname(bigwig).split(os.path.sep)[-1]
        name = '{}_{}_{}'.format(
            sample_name, fragment_length,
            trackhub.helpers.sanitize(os.path.basename(bigwig)))
        negate_values = 'off'
        if 'neg' in name:
            negate_values = 'on'
        track = trackhub.Track(
            name=name,  # track names can't have any spaces or special chars.
            source=bigwig,  # filename to build this track from
            visibility='full',  # shows the full signal
            color='128,0,5',  # brick red
            autoScale='on',  # allow the track to autoscale
            tracktype='bigWig',  # required when making a track
            negateValues=negate_values)
        trackdb.add_tracks(track)
    trackhub.upload.upload_hub(hub=hub, host=host, remote_dir=remote_dir)


def get_compositetrack_text(track_name, parent):
    header = """\
    track {0}
    compositeTrack on
    group {0}
    parent {1}
    shortLabel {0}
    longLabel {0}
    visibility full
    \n""".format(track_name, parent)
    return dedent(header)


def get_multiwigtrack_text(track_name, parent):
    """Create a multiWig track.

    Example
    track myMultiWig
    container multiWig
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    type bigWig 0 1000
    viewLimits 0:10
    maxHeighPixels 100:32:8

        track myFirstOverlaySig
        parent myMultiWig
        color 255,128,128
        type bigWig 0 1139

        track myFirstBigWig
        parent myMultiWig
        color 120,235,204

    """
    header = """\
    track {0}
    container multiWig
    parent {1}
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    maxHeightPixels 500:100:8
    viewLimits 0:20
    shortLabel {0}
    longLabel {0}
    \n""".format(track_name, parent)
    return dedent(header)


def get_supertrack_text(track_name):
    header = """\
    track {0}
    superTrack on
    group {0}
    shortLabel {0}
    longLabel {0}
    visibility full
    \n""".format(track_name)
    return dedent(header)


def get_bigwigtrack_text(track_name, parent, big_data_url):
    text = """\
    track {0}
    parent {1}
    priority 300
    longLabel {0}
    shortLabel {0}
    type bigWig
    itemRgb On
    windowingFunction mean
    autoScale on
    gridDefault on
    color 24,90,197
    visibility full
    bigDataUrl {2}
    \n""".format(track_name, parent, big_data_url)
    return dedent(text)


def create_trackdb(bwdir, srp):
    """bwdir is the root directory
    where all fragment lengths sit

    bwdir: string


    SRPXXX
          |_____SRXYYYY
          |________29
          |         | __5prime_pos.bw
          |         | __5prime_
          |
          |________28


    """

    bw_types = ['5prime_pos', '5prime_neg', '3prime_pos', '3prime_neg']
    # Step 0. Create super track

    srp_header = get_supertrack_text(srp)
    master_text = '###########################################\n'+srp_header
    # Step 1. Create composite track for SRX
    for srx in os.listdir(bwdir):

        master_text += '\n\n###############SRXHeader begin####################\n\n'
        srx_header = indent(get_compositetrack_text(srx, srp), INDENT_SPACES)
        master_text += srx_header
        # Step 2. Inside each composite track for a SRX, create another
        # composite track for different orientations/strand inside which
        # we need another multiwig track comprising all fragments for this
        # particular orientation/strand
        for orientation_strand in bw_types:
            orientation_header = indent(get_compositetrack_text(
                srx + '_' + orientation_strand, srx), 2*INDENT_SPACES)
            multiwig_header = indent(get_multiwigtrack_text(
                srx+ '_' + orientation_strand + '_multiWig', srx + '_' + orientation_strand), 3*INDENT_SPACES)
            bigwig_text = ''
            for read_length in sorted(os.listdir(os.path.join(bwdir, srx))):
                bwpath = os.path.join(bwdir, srx, read_length,
                                      orientation_strand + '.bw')
                if not os.path.isfile(bwpath):
                    sys.stderr.write(bwpath+'\n')
                    continue
                if os.stat(bwpath) == 0:
                    continue
                track_name = '{}_{}_{}'.format(srx, read_length,
                                               orientation_strand)
                # Step 2a: Create multiwig track
                bwpath = bwpath.replace('/staging/as/skchoudh/re-ribo-analysis/', 'http://smithlab.usc.edu/lab/public/skchoudh/riboraptor_trackhub/hg38/').replace('mapped/', '')
                bigwig_text += indent(get_bigwigtrack_text(
                    track_name, srx + '_' + orientation_strand + '_multiWig', bwpath), 4*INDENT_SPACES)
            master_text += '\n\n' + orientation_header + '\n\n' + multiwig_header + '\n\n' + bigwig_text
        master_text += '\n\n###############SRXHeader end####################'
    return master_text


if __name__ == '__main__':
    print(create_trackdb(sys.argv[1], sys.argv[2]))
