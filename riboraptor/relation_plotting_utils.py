import numpy as np
from scipy import stats
import seaborn as sns


def r2_pearson(x, y):
    return stats.pearsonr(x, y)[0] ** 2


def r_pearson(x, y):
    return stats.pearsonr(x, y)[0]


def r2_spearman(x, y):
    return stats.spearmanr(x, y)[0] ** 2


def r_spearman(x, y):
    return stats.spearmanr(x, y)[0]


def plot_jointplot(
    x,
    y,
    df,
    kind="hex",
    correlation="spearman",
    logx=False,
    logy=False,
    reg_line_color="#e34a33",
    point_color="#4CB391",
    alpha=0.5,
):

    if logx:
        df["log({})".format(x)] = np.log10(df[x] + 1)
        x = "log({})".format(x)
    if logy:
        df["log({})".format(y)] = np.log10(df[y] + 1)
        y = "log({})".format(y)

    if correlation not in ["spearman", "pearson"]:
        raise Exception("Valid correlation: pearson, spearman")

    if correlation == "spearman":
        r = r_spearman
    else:
        r = r_pearson
    if kind == "hex":
        g = sns.jointplot(
            x=x, y=y, data=df, kind="hex", color=point_color, height=8, stat_func=r
        )
    else:
        g = sns.jointplot(
            x=x, y=y, data=df, color=point_color, height=8, stat_func=r, alpha=alpha
        )
    g = g.annotate(
        r,
        template="{stat}: {val:.2f}",
        stat="$R$ {}".format(correlation),
        loc="upper right",
        fontsize=18,
        frameon=False,
    )
    sns.regplot(x, y, data=df, ax=g.ax_joint, scatter=False, color=reg_line_color)
    return g
