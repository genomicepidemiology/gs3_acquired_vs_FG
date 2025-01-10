from matplotlib.patches import Patch
from matplotlib import gridspec
import matplotlib.pyplot as plt

min_map_kwds = {
        "color": "white",
        "hatch": "////",
        "label": "No samples",
        "alpha": 0.5,
        "edgecolor": "lightgrey"
    }

def make_palette(labels, cmap=None, colors=None):

    if colors is not None:
        colorDict = dict(zip(labels, colors))
    elif cmap is not None:
        colorDict = dict(zip(labels, cmap))
    
    colorHandles = [
        Patch(facecolor=v, label=l) for l, v in colorDict.items()
    ]
    return colorDict, colorHandles

db_palette, db_handles = make_palette(
    labels=['resfinder', 'functional_amr', 'csabapal'], 
    colors=['#FC766AFF', '#B83B5E', '#F08A5D']
)
region_palette, region_handles = make_palette(
    labels = [
        'East Asia & Pacific', 'Europe & Central Asia', 
        'Latin America & Caribbean', 'Middle East & North Africa', 
        'North America', 'South Asia', 'Sub-Saharan Africa'
    ],
    colors = [
        '#FEFF32', 
        '#4BAD49',
        '#A75529',
        '#974EA2', 
        '#FF8000', 
        '#3A7EB5',
        '#E4191C'
    ]
)
region_order = list(region_palette.keys())


region_palette2, region_handles2 = make_palette(
    labels = [
        'East Asia & Pacific (EAP)', 'Europe & Central Asia (ECA)', 
        'Latin America & Caribbean (LAC)', 'Middle East & North Africa (MENA)', 
        'North America (NA)', 'South Asia (SA)', 'Sub-Saharan Africa (SSA)'
    ],
    colors = [
        '#FEFF32', 
        '#4BAD49',
        '#A75529',
        '#974EA2', 
        '#FF8000', 
        '#3A7EB5',
        '#E4191C'
    ]
)


dbgroup_palette, dbgroup_handles = make_palette(
    labels = ['ResFinder', 'Functional', 'Genera'],
    colors = ['#FC766AFF', '#5B84B1FF', '#1f77b4']
)

dbgroup_palette2, dbgroup_handles2 = make_palette(
    labels = ['Acquired', 'FG', 'Genera'],
    colors = ['#FC766AFF', '#5B84B1FF', '#1f77b4']
)



alphabet = list('abcdefghijklmnopqrstuvwxyz')

def plot_map(fig, subplot_spec, df, col, all_countries, subplot_label=None, **kwargs):

    height_ratios=[1]
    if kwargs.get('add_cbar', True):
        height_ratios.append(.025)

    gs_map_meta = gridspec.GridSpecFromSubplotSpec(
        nrows=1 + int(kwargs.get('add_cbar', True)), ncols=1, 
        subplot_spec=subplot_spec,
        height_ratios=height_ratios,
        hspace=-.7
    )

    ax = fig.add_subplot(gs_map_meta[0])

    vmin = kwargs.get('vmin', df[col].min())
    vmax = kwargs.get('vmax', df[col].max())

    
    if df.shape[0] > 1:
        df.plot(
            ax = ax,
            column = col,
            cmap = kwargs.get('cmap', 'Greens'),
            alpha = kwargs.get('alpha', .9),
            missing_kwds = kwargs.get('min_map_kwds', None),
        )

        all_countries = all_countries.loc[~((
                all_countries.name.isin(df.country)) |
                (all_countries.name.isin(df.country_alt))),
        ]

    all_countries.plot(
        color='white',
        edgecolor='lightgray',
        hatch='////',
        alpha = .5,
        ax = ax
    )


    if kwargs.get('add_cbar', True):
        cax = fig.add_subplot(gs_map_meta[1])

        sm = plt.cm.ScalarMappable(
            cmap=kwargs.get('cmap', 'Greens'), 
            norm=plt.Normalize(vmin=vmin, vmax=vmax)
        )
        plt.colorbar(
            sm, 
            shrink=.01,
            cax=cax,
            orientation='horizontal'
        ).set_label(label=kwargs.get('clab', ''), size=9)
        
    ax.axis('off')

    ax.set_title(kwargs.get('title', ''),fontdict={'size': 10}) 

    if subplot_label is not None:
        ax.text(-0.05, 1.05, subplot_label + '.', transform=ax.transAxes, fontsize=12, va='top', ha='right', fontweight='bold')
