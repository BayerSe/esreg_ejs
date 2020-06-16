import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
from arch.bootstrap import MCS

sns.set(style="white",
        rc={'font.family': 'serif', 'text.usetex': True,
            'axes.labelsize': 7, 'axes.titlesize': 7, 'font.size': 7, 'legend.fontsize': 7,
            'xtick.labelsize': 7, 'ytick.labelsize': 7})

markers = ["s", "^", "v", "o", "D", "."]

map_g1 = {
    1: "$G_1(z) = z$",
    2: "$G_1(z) = 0$"
}

map_g2 = {
    1: "$-\log(-z)$",
    2: "$-\sqrt{-z}$",
    3: "$-1/z$",
    4: "$\log(1 + \exp(z))$",
    5: "$\exp(z)$"
}

map_design = {
    i: 'DGP-(%s)' % i for i in range(1, 5)
}

map_cov = {
    "cov_iid_ind": "iid / ind",
    "cov_nid_scl_N": "nid / scl-N",
    "cov_nid_scl_t": "nid / scl-t",
    "cov_nid_scl_sp": "nid / scl-sp",
    "cov_true": "true",
    "cov_boot": "bootstrap"
}

g2_order = [
    "$-\log(-z)$",
    "$-\sqrt{-z}$",
    "$-1/z$",
    "$\log(1 + \exp(z))$",
    "$\exp(z)$"
]


def compute_mcs(loss, reps, block_size, bootstrap, statistic):
    mcs = MCS(losses=loss, size=0.05, bootstrap=bootstrap, method=statistic, reps=reps, block_size=block_size)
    mcs.seed(1)
    mcs.compute()
    return mcs.pvalues


def figure_size(width, height, horizontal_panels, vertical_panels):
    height_each_panel = height / horizontal_panels
    width_each_panel = width / vertical_panels
    aspect = width_each_panel / height_each_panel
    return height_each_panel, aspect


def plot_murphy_diagram(loss_mean, loss_sd, cl=0.95, file="output/murphy_diagram.pdf"):

    loss_mean.set_index('x', inplace=True)
    loss_sd.set_index('x', inplace=True)

    n_models = loss_mean.shape[1]
    fig, axes = plt.subplots(nrows=1, ncols=n_models, figsize=(1.5 * n_models, 1.3))

    for i, ax in enumerate(axes):
        m = loss_mean.iloc[:, i]
        se = loss_sd.iloc[:, i]
        q = scipy.stats.norm.ppf((1+cl)/2)
        low, up = m - q * se, m + q * se
        m.plot(ax=ax)
        ax.fill_between(loss_mean.index, low, up, alpha=0.3)
        ax.axhline(0, color=".2", ls="--", lw=1)

    lim = pd.DataFrame([ax.get_ylim() for ax in axes])
    for i, ax in enumerate(axes):
        ax.set_ylim((lim.iloc[:, 0].min(), lim.iloc[:, 1].max()))
        ax.locator_params(tight=True, nbins=5, axis="y")
        ax.locator_params(tight=True, nbins=4, axis="x")
        ax.set_xlabel("Threshold")
        if i == 0:
            ax.set_ylabel("Score Difference")
        else:
            ax.set_yticklabels([])
        ax.set_title("Difference to " + loss_mean.columns[i], y=0.9)

    sns.despine()
    plt.tight_layout(pad=0.1)

    plt.savefig(file)
    plt.close("all")


def plot_mse(df, file):
    df = df.astype({'n': int})
    df['g1'].replace(map_g1, inplace=True)
    df['g2'].replace(map_g2, inplace=True)
    df['design'].replace(map_design, inplace=True)

    height_each_panel, aspect = figure_size(width=4, height=5, horizontal_panels=4, vertical_panels=2)

    g = sns.catplot(x="n", y="value", row="design", col='g1', hue="g2", data=df, kind='point', sharey='row',
                    sharex=False,
                    legend=False, markers=markers, height=height_each_panel, aspect=aspect, scale=0.7)
    g.set_axis_labels("Sample Size", "MSE")
    g.set_titles("{row_name} $\mid$ {col_name}")

    for ax in g.axes.flatten():
        plt.setp(ax.collections, sizes=[30], zorder=100, edgecolor=["black"], lw=[0.3])
        ax.get_yaxis().set_label_coords(-0.2, 0.5)
        ax.yaxis.grid()
        ax.locator_params(nbins=5, axis="y")
        ax.set_ylim(0, ax.get_ylim()[1])

    g.fig.legend(handles=g._legend_data.values(), labels=g._legend_data.keys(),
                 loc='upper center', ncol=5, framealpha=0.5, title="$\mathcal{G}_2(z) = $",
                 columnspacing=0.1, labelspacing=0.1, handletextpad=0.1)

    g.fig.subplots_adjust(top=0.88, bottom=0.08, left=0.11, right=0.96, wspace=0.2, hspace=1)

    plt.savefig(file)
    plt.close("all")


def plot_norm(df, g1, file):
    df = df.astype({'n': int})
    df = df[df['g1'] == g1]
    df = df[df.g2.isin([1, 2, 3])]

    df['g1'].replace(map_g1, inplace=True)
    df['g2'].replace(map_g2, inplace=True)
    df['design'].replace(map_design, inplace=True)
    df['cov'].replace(map_cov, inplace=True)
    df['g2'] = '$\mathcal{G}_2(z)=' + df.g2.str[1:]

    height_each_panel, aspect = figure_size(width=4.75, height=6, horizontal_panels=4, vertical_panels=3)

    g = sns.catplot(x='n', y='value', row='design', col='g2', hue='cov', data=df, kind='point',
                    height=height_each_panel, aspect=aspect, sharey='row', sharex=True, scale=0.7,
                    legend=False, markers=markers)

    g.set_axis_labels('Sample Size', 'Norm')
    g.set_titles('{row_name} \n {col_name}')

    for ax in g.axes.flatten():
        plt.setp(ax.collections, sizes=[30], zorder=100, edgecolor=["black"], lw=[0.3])
        ax.yaxis.grid()
        ax.locator_params(nbins=5, axis='y')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

    g.fig.legend(handles=g._legend_data.values(), labels=g._legend_data.keys(),
                 loc='upper center', ncol=4, framealpha=0.5, title='Covariance Estimator',
                 columnspacing=1, labelspacing=0.2, handletextpad=0.1)

    g.fig.subplots_adjust(top=0.88, left=0.11, right=0.99, wspace=0.2, hspace=0.6)

    plt.savefig(file)
    plt.close('all')


def plot_estimation_time(df, file):
    df = df.astype({'n': int})
    df['g1'].replace(map_g1, inplace=True)
    df['g2'].replace(map_g2, inplace=True)
    df['design'].replace(map_design, inplace=True)
    df = df.groupby(['n', 'g1', 'g2']).mean().reset_index()

    height_each_panel, aspect = figure_size(width=4.75, height=1.5, horizontal_panels=1, vertical_panels=2)

    g = sns.catplot(x='n', y='time', hue='g2', col='g1', kind='point', data=df, dodge=True,
                    col_order=[map_g1[i] for i in [1, 2]], hue_order=g2_order,
                    height=height_each_panel, aspect=aspect, scale=0.7, legend=False, markers=markers)

    g.set_axis_labels('Sample Size', 'Avg. Est. Time (s)')
    g.set_titles('{col_name}')

    for ax in g.axes.flatten():
        plt.setp(ax.collections, sizes=[30], zorder=100, edgecolor=["black"], lw=[0.3])
        ax.yaxis.grid()
        ax.locator_params(nbins=5, axis="y")

    g.fig.legend(handles=g._legend_data.values(), labels=g._legend_data.keys(),
                 loc='right', ncol=1, framealpha=0, title="$\mathcal{G}_2(z) = $")
    g.fig.subplots_adjust(top=0.88, bottom=0.25, left=0.08, right=0.75, wspace=0.1, hspace=0.6)

    plt.savefig(file)
    plt.close('all')


def plot_mse_decomposed(df, design, sample_size, file):
    df1 = df.query('design == "%s" & n == %s & g1 == %s' % (design, sample_size, 1)).pivot('g2', 'par', 'value')
    df2 = df.query('design == "%s" & n == %s & g1 == %s' % (design, sample_size, 2)).pivot('g2', 'par', 'value')

    n_par = df1.shape[1]
    df1.index = [map_g2[i] for i in df1.index]
    df2.index = df1.index

    df1.columns = ['$\\beta^q_%s$' % i for i in range(int(n_par / 2))] + \
                  ['$\\beta^e_%s$' % i for i in range(int(n_par / 2))]
    df2.columns = df1.columns

    dfall = [df1, df2]
    labels = [map_g1[i] for i in [1, 2]]

    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)
    color = sns.color_palette("Blues", int(n_col / 2)) + sns.color_palette("Reds", int(n_col / 2))

    f = plt.figure(figsize=(4, 2))
    axe = f.add_subplot(111)

    for i, df in enumerate(dfall):
        axe = df.plot(ax=axe, kind="bar", linewidth=0.5, edgecolor='k', stacked=True,
                      legend=False, grid=False, color=color)
    axe.set_ylabel('Mean Squared Error')
    axe.set_xlabel('$\mathcal{G}_2(z)$')
    axe.yaxis.grid()

    hatch = "//"

    h,l = axe.get_legend_handles_labels()
    for i in range(0, n_df * n_col, n_col):
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches:
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(hatch * int(i / n_col))
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation=0)

    n = []
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=hatch * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc='upper left', ncol=2, framealpha=1)
    if labels is not None:
        if n_col == 4:
            l2 = plt.legend(n, labels, loc='upper left', bbox_to_anchor=(0, 0.7), framealpha=1)
        elif n_col == 6:
            l2 = plt.legend(n, labels, loc='upper left', bbox_to_anchor=(0, 0.6), framealpha=1)
    axe.add_artist(l1)

    sns.despine()
    plt.tight_layout(pad=0.1)

    plt.savefig(file)
    plt.close('all')


def plot_relse(data, g1, g2, sample_size, file):
    data = data.astype({'design': int})
    data['cov'].replace(map_cov, inplace=True)
    df = data.query('n == %s & g1 == %s & g2 == %s' % (sample_size, g1, g2)).copy()

    for design in df['design'].unique():
        n_par = df.loc[df['design'] == design, 'par'].max()
        d = {}
        for i in range(int(n_par/2)):
            d[i+1] = '$\\beta^q_%s$' % i
        for i in range(int(n_par / 2), n_par):
            d[i+1] = '$\\beta^e_%s$' % (i - int(n_par / 2))

        df.loc[df['design'] == design, 'par'] = df.loc[df['design'] == design, 'par'].replace(d)

    height_each_panel, aspect = figure_size(width=4.75, height=2, horizontal_panels=1, vertical_panels=4)

    g = sns.catplot(x='par', y='value', hue='cov', col='design', data=df,
                    kind='point', height=height_each_panel, aspect=aspect, sharex=False,
                    legend=False, markers=markers, scale=0.7)

    g.set_axis_labels('Parameter', 'RELSE')
    g.set_titles('DGP-({col_name})')

    for i, ax in enumerate(g.axes.flatten()):
        ax.axhline(1, color='0.5')
        plt.setp(ax.collections, sizes=[30], zorder=100, edgecolor=["black"], lw=[0.3])
        ax.yaxis.grid()
        ax.locator_params(nbins=7, axis='y')
        if i in [0, 1, 3]:
            ax.set_xticks([0, 1, 2, 3])
            ax.set_xlim(-0.5, 3.5)

    g.fig.legend(handles=g._legend_data.values(), labels=g._legend_data.keys(),
                 loc='upper center', ncol=4, framealpha=0.5, title='Covariance Estimator')
    g.fig.subplots_adjust(top=0.7, bottom=0.2, left=0.1, right=0.98, wspace=0.1, hspace=0.6)

    plt.savefig(file)
    plt.close('all')


def plot_relse_all(data, g1, g2, file):
    data = data.astype({'design': int, 'n': int})
    data['cov'].replace(map_cov, inplace=True)
    df = data.query('g1 == %s & g2 == %s' % (g1, g2)).copy()

    for design in df['design'].unique():
        n_par = df.loc[df['design'] == design, 'par'].max()
        d = {}
        for col_idx in range(int(n_par/2)):
            d[col_idx+1] = '$\\beta^q_%s$' % col_idx
        for col_idx in range(int(n_par / 2), n_par):
            d[col_idx+1] = '$\\beta^e_%s$' % (col_idx - int(n_par / 2))

        df.loc[df['design'] == design, 'par'] = df.loc[df['design'] == design, 'par'].replace(d)

    height_each_panel, aspect = figure_size(width=4.75, height=6, horizontal_panels=5, vertical_panels=4)

    g = sns.catplot(x='par', y='value', hue='cov', col='design', row='n', data=df,
                    kind='point', height=height_each_panel, aspect=aspect, sharex=False, sharey='row',
                    legend=False, markers=markers, scale=0.7)

    g.set_axis_labels('Parameter', 'RELSE')
    g.set_titles('DGP-({col_name}) \n $n={row_name}$')

    for row_idx, row_ax in enumerate(g.axes):
        for col_idx, ax in enumerate(row_ax):
            ax.axhline(1, color='0.5')
            ax.set_title(ax.get_title(), pad=-0.4)
            plt.setp(ax.collections, sizes=[30], zorder=100, edgecolor=["black"], lw=[0.3])
            ax.yaxis.grid()
            ax.locator_params(nbins=5, axis='y')
            if col_idx in [0, 1, 3]:
                ax.set_xticks([0, 1, 2, 3])
                ax.set_xlim(-0.5, 3.5)
            if row_idx in [1, 2, 3, 4]:
                ax.set_yticks([0.5, 0.75, 1, 1.25, 1.5])

    g.fig.legend(handles=g._legend_data.values(), labels=g._legend_data.keys(),
                 loc='upper center', ncol=4, framealpha=0.5, title='Covariance Estimator')
    g.fig.subplots_adjust(top=0.9, bottom=0.08, left=0.11, right=0.98, wspace=0.2, hspace=0.8)

    plt.savefig(file)
    plt.close('all')
