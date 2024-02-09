#Code to make the png images of fundamental domains. Modified based on David Lowry-Duda's LMFDB code.
gp.read("fdom.gp")
from sage.misc.decorators import options, rename_keyword
from sage.plot.hyperbolic_polygon import HyperbolicPolygon
from sage.libs.pari.convert_sage import gen_to_sage


#SECTION 1: GRAPHICS


#To encode the png file as a string.
def encode_mcurve_plot(P, transparent = True):
    from io import BytesIO as IO
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from base64 import b64encode
    from urllib.parse import quote

    virtual_file = IO()
    fig = P.matplotlib(axes_pad=None)
    ax = fig.axes[0]
    ax.set_xlim(xmin = -1, xmax = 1)
    ax.set_ylim(ymin = -1, ymax = 1)
    fig.set_canvas(FigureCanvasAgg(fig))
    fig.set_size_inches(2.5, 2.5) # images are 200 x 200 on the website
    fig.savefig(virtual_file, format = 'png', bbox_inches = 'tight', transparent = transparent, dpi = 120)
    virtual_file.seek(0)
    buf = virtual_file.getbuffer()
    return "data:image/png;base64," + quote(b64encode(buf))

#To make a hyperbolic polygon given the vertices.
@rename_keyword(color='rgbcolor')
@options(alpha = 1, fill = True, thickness = 0, rgbcolor = "blue", zorder = 2, linestyle = 'solid')
def my_hyperbolic_polygon(pts, model = "PD", resolution = 200, circlecolor = 'black', add_circle = False, **options):
    r"""
    Return a hyperbolic polygon in the hyperbolic plane with vertices ``pts``.

    This is on the poincare disk, but doesn't add random black circles
    repeatedly like the default version.

    INPUT:

    - ``pts`` -- a list or tuple of complex numbers

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) the style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or
      ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively
    """
    if not model == "PD":
        raise ValueError("Only the disk model is supported")
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))

    g.add_primitive(HyperbolicPolygon(pts, model, options))
    g.set_aspect_ratio(1)
    if add_circle:
        g = g + circle((0, 0), 1, rgbcolor=circlecolor)
    return g

#Given the vertices for one fundamental domain, makes the domain as a Graphics object.
def onefdomgraphicsobject(verts):
    g = Graphics()
    g += my_hyperbolic_polygon(verts, add_circle = True);
    g.axes(False)
    g.set_axes_range(-1, 1, -1, 1)
    return g


#SECTION 2: FUNDAMENTAL DOMAINS


#Return the fundamental domain for D
def level1fdom(D):
    gp.setrand(1)
    A = gp.alginit_Qdisc(D)
    gp.setrand(1)
    X = gp.afuchinit(A);
    return X
    
#Make the level 1 pictures for all D in Dlist. We need to also input the size of Aut_mu(O), so the pairs are [D, #Aut_mu] where mu is a degree 1 polarization.
def make_pictures_level1(Dlistwithmu):
    fil = open("../data/level1_pictures.fdom", "a")
    for pair in Dlistwithmu:
        D = pair[0]
        sizeautmu = pair[1]
        X = level1fdom(D)
        verts = gen_to_sage(pari(gp.afuchvertices(X, 1)));
        g = onefdomgraphicsobject(verts)
        genus = gp.afuchsignature(X)[1]
        pngstr = encode_mcurve_plot(g)
        fil.write(f"{D}.1.1.{sizeautmu}.{genus}.a.1?{pngstr}\n")
    fil.close()


#SECTION 3: PREPROCESSING TABLES


"""
Converts the quaternion orders saved in ../data/quaternion-orders.m into quaternion-orders-gp.dat, so that the data can be gp-read.
We save the order as: [discB, Olevel, [a, b], [bas1, bas2, bas3, bas4]]
discB = discriminant of the quaternion algebra
Olevel = level of the order O (often 1)
[a, b] = the quaternion algebra is B = (a, b / Q)
basi = A basis for O is [bas1, ..., bas4], where basi gives the basis element in terms of the basis (1, i, j, k).

The data is stored in "./quaternion-orders-gp.dat", and is sorted first by discB, then by Olevel.
"""
def make_gp_orders():
    f1 = open("../data/quaternion-orders/quaternion-orders.m", "r")
    lineno = 0
    allords = []
    for line in f1:
        lineno += 1
        if lineno <= 3:
            continue
        x = line.split("?")
        a = int(x[1])
        b = int(x[2])
        discO = int(x[3])
        discB = int(x[4])
        Olevel = int(discO / discB)
        bn = list(map(ZZ, x[5].replace("{", "").replace("}", "").split(",")))
        bd = list(map(ZZ, x[6].replace("{", "").replace("}", "").split(",")))
        bas1 = [bn[0]/bd[0], bn[1]/bd[0], bn[2]/bd[0], bn[3]/bd[0]]
        bas2 = [bn[4]/bd[1], bn[5]/bd[1], bn[6]/bd[1], bn[7]/bd[1]]
        bas3 = [bn[8]/bd[2], bn[9]/bd[2], bn[10]/bd[2], bn[11]/bd[2]]
        bas4 = [bn[12]/bd[3], bn[13]/bd[3], bn[14]/bd[3], bn[15]/bd[3]]
        entry = [discB, Olevel, [a, b], [bas1, bas2, bas3, bas4]]
        allords.append(entry)
    f1.close()
    allords = sorted(allords, key = lambda element: (element[0], element[1]))
    f2 = open("quaternion-orders-gp.dat", "w")
    for x in allords:
        f2.write(f"{x}\n")
    f2.close()

"""
Converts the polarized quaternion orders saved in ../data/quaternion-orders-polarized.m into quaternion-orders-polarized-gp.dat, so that the data can be gp-read.
We save the polarized order as: [discB, Olevel, degmu, cardautmuO, mu, gens]
discB = discriminant of the quaternion algebra
Olevel = level of the order O (often 1)
mu = the polarizing element, written in terms of the basis (1, i, j, k)
degmu = degree of mu, i.e. mu^2=-discB*Olevel*degmu
cardautmuO = Cardinality of Aut_{+/-mu}(O), should be even and at most 12?
gens = generators of Aut_{+/-mu}(O), taken in the basis (1, i, j, k), with norms dividing discB (at least for maximal).

The data is stored in "./quaternion-orders-polarized-gp.dat", and is sorted first by discB, then Olevel, and then degmu.
"""
def make_gp_polarized_orders():
    f1 = open("../data/quaternion-orders/quaternion-orders-polarized.m", "r")
    lineno = 0
    allpolords = []
    for line in f1:
        lineno += 1
        if lineno <= 3:
            continue
        x = line.split("?")
        #The next line will break with non-maximal orders once added. At this point, fix this method as appropriate
        discO = int(x[1])
        mu = list(map(ZZ, x[2].replace("{", "").replace("}", "").split(",")))
        degmu = int(x[3])
        cardautmuO = int(x[5])
        allgens = list(map(ZZ, x[8].replace("{", "").replace("}", "").split(",")))
        gens = []
        ind = -1
        for ent in allgens:
            ind += 1
            ind %= 4
            if ind == 0:
                elt = [ent, 0, 0, 0]
                continue
            elt[ind] = ent
            if ind == 3:
                gens.append(elt)
        entry = [discO, 1, degmu, cardautmuO, mu, gens]
        allpolords.append(entry)
    f1.close()
    allpolords = sorted(allpolords, key = lambda element: (element[0], element[1], element[2]))
    f2 = open("quaternion-orders-polarized-gp.dat", "w")
    for x in allpolords:
        f2.write(f"{x}\n")
    f2.close()

