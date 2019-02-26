import json
import re
import logging as log
from .helpers import path_leaf

regexes = {
    "bp_processed": "Total basepairs processed:\s*([\d,]+) bp",
    "bp_written": "Total written \(filtered\):\s*([\d,]+) bp",
    "quality_trimmed": "Quality-trimmed:\s*([\d,]+) bp",
    "r_processed": "Total reads processed:\s*([\d,]+)",
    "r_with_adapters": "Reads with adapters:\s*([\d,]+)",
}


def cutadapt_to_json(filepath, savetofile=None):
    """Convert cutadapt/trim_galore output to json

    Parameters
    ----------
    filepath: string
              Path to trim_galore/cutadapt output.txt

    Returns
    -------
    json_data: dict
    """
    fh = open(filepath, "r")
    trim_info = {}
    length_counts = {}
    length_exp = {}
    length_obsexp = {}
    adapters = {}
    sample = None
    for l in fh:
        if "cutadapt" in l:
            sample = None
        if l.startswith("Used user"):
            # Used user provided input and hence no second pass
            adapters = "User provided"
            break
        if l.startswith("No adapter"):
            adapters = "None found (second pass)"
            break
        if l.startswith("Command line parameters"):
            sample = l.split()[-1]
            sample = path_leaf(sample).replace(".fq.gz", "").replace(".fastq.gz", "")
            if sample in trim_info:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample))
            trim_info[sample] = dict()
        if sample is not None:
            for k, r in list(regexes.items()):
                match = re.search(r, l)
                if match:
                    trim_info[sample][k] = int(match.group(1).replace(",", ""))

            if "===" in l:
                log_section = l.strip().strip("=").strip()
            if l.startswith("Sequence:"):
                plot_sname = "{} - {}".format(sample, log_section)
                adapters[plot_sname] = l.split(";")[0].strip("Sequence: ")

            if "length" in l and "count" in l and "expect" in l:
                plot_sname = sample
                if log_section is not None:
                    plot_sname = "{} - {}".format(sample, log_section)
                length_counts[plot_sname] = dict()
                length_exp[plot_sname] = dict()
                length_obsexp[plot_sname] = dict()
                for l in fh:
                    r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                    if r_seqs:
                        a_len = int(r_seqs.group(1))
                        length_counts[plot_sname][a_len] = int(r_seqs.group(2))
                        length_exp[plot_sname][a_len] = float(r_seqs.group(3))
                        if float(r_seqs.group(3)) > 0:
                            length_obsexp[plot_sname][a_len] = float(
                                r_seqs.group(2)
                            ) / float(r_seqs.group(3))
                        else:
                            length_obsexp[plot_sname][a_len] = float(r_seqs.group(2))
                    else:
                        break
    fh.close()
    json_data = {
        "adapters": adapters,
        "trim_info": trim_info,
        "length_exp": length_exp,
        "length_obsexp": length_obsexp,
        "length_counts": length_counts,
    }
    if savetofile:
        json.dump(json_data, savetofile)
    return json_data
