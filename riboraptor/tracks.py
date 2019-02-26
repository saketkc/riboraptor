import os
from textwrap import indent, dedent
import sys

INDENT_SPACES = "\t"


def get_compositetrack_text(track_name, parent):
    """Create composite track text"""
    header = """\
    track {0}
    compositeTrack on
    group {0}
    parent {1}
    shortLabel {0}
    longLabel {0}
    visibility full
    \n""".format(
        track_name, parent
    )
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
    visibility full
    type bigWig
    \n""".format(
        track_name, parent
    )
    return dedent(header)


def get_supertrack_text(track_name):
    header = """\
    track {0}
    superTrack on
    group {0}
    shortLabel {0}
    longLabel {0}
    visibility full
    \n""".format(
        track_name
    )
    return dedent(header)


def get_bigwigtrack_text(track_name, parent, big_data_url, negate_values):
    """Create bigwig track text"""
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
    negateValues {3}
    bigDataUrl {2}
    \n""".format(
        track_name, parent, big_data_url, negate_values
    )
    return dedent(text)


def create_trackdb(bwdir, srp, orientation="5prime"):
    """Create track file
    """

    # Step 0. Create super track
    srp_header = get_supertrack_text(srp)
    master_text = "###########################################\n" + srp_header
    # Step 1. Create composite track for SRX
    for srx in sorted(os.listdir(bwdir)):
        master_text += "\n\n###############SRXHeader begin####################\n\n"
        srx_orientation_key = srx + "_" + orientation + "_multiWig"
        multiwig_header = indent(
            get_multiwigtrack_text(srx_orientation_key, srp), INDENT_SPACES
        )
        # Step 2. Inside each composite track for a SRX, create another
        # composite track for different orientations/strand inside which
        # we need another multiwig track comprising all fragments for this
        # particular orientation/strand
        add_header = 0
        for read_index, read_length in enumerate(
            sorted(os.listdir(os.path.join(bwdir, srx)))
        ):
            bigwig_text = ""
            for orientation_index, orientation_strand in enumerate(
                [orientation + "_pos", orientation + "_neg"]
            ):
                bwpath = os.path.join(
                    bwdir, srx, read_length, orientation_strand + ".bw"
                )
                track_name = "{}_{}_{}".format(srx, read_length, orientation_strand)
                # Step 2a: Create multiwig track
                if "neg" in bwpath:
                    negate_values = "on"
                else:
                    negate_values = "off"
                if os.path.isfile(bwpath) and os.stat(bwpath).st_size:
                    add_header += 1
                    if add_header == 1:
                        # Add header only if there is atleast one
                        # file with atleast one strand
                        master_text += multiwig_header + "\n\n"
                    bigwig_text += indent(
                        get_bigwigtrack_text(
                            track_name,
                            srx_orientation_key,
                            "http://ribopod.usc.edu" + bwpath,
                            negate_values,
                        ),
                        2 * INDENT_SPACES,
                    )
            master_text += "\n\n" + bigwig_text
        master_text += "\n\n###############SRXHeader end####################"
    return master_text


if __name__ == "__main__":
    print(create_trackdb(sys.argv[1], sys.argv[2]))
