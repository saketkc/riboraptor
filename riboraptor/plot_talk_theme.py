import matplotlib
import seaborn as sns

PAPER_FONTSIZE = 22


def load_pub_theme(label_font_size=PAPER_FONTSIZE):
    PAPER_PRESET = {
        "style": "ticks",
        "font": "Helvetica",
        "context": "paper",
        "rc": {
            "font.size": label_font_size,
            "axes.titlesize": label_font_size,
            "axes.labelsize": label_font_size,
            "axes.linewidth": 2,
            "legend.fontsize": label_font_size,
            "xtick.labelsize": label_font_size,
            "ytick.labelsize": label_font_size,
            "xtick.major.size": 8.0,
            "ytick.major.size": 8.0,
            "axes.edgecolor": "black",
            "xtick.major.pad": 3.0,
            "ytick.major.pad": 3.0,
        },
    }
    sns.set(**PAPER_PRESET)
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
