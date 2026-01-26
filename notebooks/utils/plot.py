import seaborn as sns


def haar_scatter(ax, x, y):
    sns.set(style="whitegrid", font_scale=1.2)

    # scatter points colored by value with better aesthetics
    ax.scatter(
        x, y,
        s=40,
        alpha=0.8,
        edgecolors='black',
        linewidth=0.3,
        zorder=3
    )

    # draw vertical lines behind scatter points
    ax.vlines(
        x=x,
        ymin=0,
        ymax=y,
        colors='gray',
        alpha=0.5,
        linewidth=1,
        zorder=1
    )

    ax.yaxis.grid(True, linestyle='--', linewidth=0.7, alpha=0.7)
    ax.xaxis.grid(False)

    sns.despine(ax=ax)