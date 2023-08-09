import matplotlib.colors as mcolors


def make_colormap(seq):
    """
    Parameters
    ----------
    seq : 
        a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
    
    Returns
    -------
    LinearSegmentedColormap
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {"red": [], "green": [], "blue": []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict["red"].append([item, r1, r2])
            cdict["green"].append([item, g1, g2])
            cdict["blue"].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap("CustomMap", cdict)


def get_colormap(objname):
    """
    Custom colormaps for making pretty pictures for Twilight Zone website
    
    Parameters
    ----------
    objname : str, required. name of target. currently supported:
        io, neptune, titan, uranus
    """
    if objname.lower() == "io":
        objmap = "gist_heat"
    elif objname.lower() == "neptune":
        c = mcolors.ColorConverter().to_rgb
        rvb = make_colormap(
            [c("black"), c("blue"), 0.33, c("blue"), c("white"), 0.66, c("white")]
        )
        objmap = rvb
    elif objname.lower() == "titan":
        c = mcolors.ColorConverter().to_rgb
        rvb = make_colormap(
            [
                c("black"),
                c("goldenrod"),
                0.6,
                c("goldenrod"),
                c("white"),
                0.95,
                c("white"),
            ]
        )
        objmap = rvb
    elif objname.lower() == "uranus":
        cut1 = 0.33
        cut2 = 1.0
        c = mcolors.ColorConverter().to_rgb
        rvb = make_colormap(
            [c("black"), c("purple"), cut1, c("purple"), c("white"), cut2, c("white")]
        )
        objmap = rvb
    return objmap
