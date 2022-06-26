
from matplotlib import rcParams

# rcParams.keys()

rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{bm}'

rcParams['lines.linewidth'] = 1.2

rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 15

rcParams['axes.linewidth'] = 0.8

rcParams['xtick.top'] = True
rcParams['xtick.direction'] = 'in'
rcParams['xtick.minor.visible'] = True
rcParams['xtick.major.size'] = 6
rcParams['xtick.minor.size'] = 3
rcParams['xtick.major.width'] = 0.75
rcParams['xtick.minor.width'] = 0.75
rcParams['xtick.major.pad'] = 5
rcParams['xtick.minor.pad'] = 5

rcParams['ytick.right'] = True
rcParams['ytick.direction'] = 'in'
rcParams['ytick.minor.visible'] = True
rcParams['ytick.major.size'] = 6
rcParams['ytick.minor.size'] = 3
rcParams['ytick.major.width'] = 0.75
rcParams['ytick.minor.width'] = 0.75
rcParams['ytick.major.pad'] = 5
rcParams['ytick.minor.pad'] = 5

rcParams['legend.fontsize'] = rcParams['font.size']
rcParams['legend.labelspacing'] = 0.2
rcParams['legend.loc'] = 'upper left'
rcParams['legend.frameon'] = False

rcParams['figure.figsize'] = (8.0, 5.0)
rcParams['figure.dpi'] = 150

rcParams['savefig.dpi'] = 200
rcParams['savefig.bbox'] = 'tight'
