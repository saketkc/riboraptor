import glob, os
import trackhub


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
