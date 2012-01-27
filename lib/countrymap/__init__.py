
from os.path import dirname, exists, join, abspath

from matplotlib.collections import LineCollection

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap.pyproj import Geod
from mpl_toolkits.basemap.shapefile import Reader as ShapeReader


__version__ = '0.0.1'


__all__ = [
    'colors',
    'alphas',
    'linewidths',
    'shapefile2data',
    'data2isodict',
    'shape2collection',
    'draw_gc',
    'CountryMap'
]


# based on google colors, for consistency
colors = {
    'park_color': '#c9dfaf',
    'land_color': '#f4f3f0',
    'lake_color': '#a5bfdd',
    'edge_color': 'k',
    'border_color': '#9d9d9d',
    'fade_color': '#d3d2d1' # '#c6c5c2'
}


linewidths = {
    'border_width': 0.09375,
    'edge_width': 0.09375,
}


alphas = {
    'edge_alpha': 0.75,
    'border_alpha': 1.
}


_dead_countries = {
    'CS': ['ME', 'RS'],
    'YU': ['BA', 'HR', 'ME', 'MK', 'RS', 'SI']
}


def shapefile2data(shapefile):
    shf = ShapeReader(shapefile)
    fields = [f[0] for f in shf.fields if isinstance(f, list)]
    records = []
    for record in shf.shapeRecords():
        datadict = dict(zip(fields, record.record))
        datadict['SHAPE'] = [record.shape]
        records.append(datadict)
    return records


def data2isodict(data):
    new2dead = {}
    for dead, news in _dead_countries.items():
        for new in news:
            if new not in new2dead:
                new2dead[new] = []
            new2dead[new].append(dead)
    isodict = {}
    for row in data:
        if 'ISO2' not in row:
            raise ValueError('no ISO2 field found, shapedata cannot be organized by ISO code!')
        isocode = row['ISO2']
        isodict[isocode] = [row]
        if isocode in new2dead:
            deadcodes = new2dead[isocode]
            for deadcode in deadcodes:
                if deadcode not in isodict:
                    isodict[deadcode] = []
                isodict[deadcode].append(row)
    return isodict


def shape2collection(m, shapes):
    coords = []
    for shape in shapes:
        shapetype = shape.shapeType
        verts = shape.points
        if shapetype in (3, 5):
            parts = shape.parts.tolist()
            for idx1, idx2 in zip(parts, parts[1:] + [len(verts)]):
                lons, lats = list(zip(*verts[idx1:idx2]))
                if (max(lons) >  721. or
                    min(lons) < -721. or
                    max(lats) >   91. or
                    min(lats) <  -91.):
                    raise ValueError('shape must have lat/lon vertices')
                x, y = m(lons, lats)
                coords.append(list(zip(x, y)))
        else:
            raise ValueError('shape2collection only supports 2D shapes')
    return LineCollection(coords, antialiaseds=(1,))


def draw_gc(m, lon1, lat1, lon2, lat2, del_s=20., **kwargs):
    gc = Geod(a=m.rmajor, b=m.rminor)
    _, _, dist = gc.inv(lon1, lat1, lon2, lat2)
    npoints = int((0.5 * 1000. * del_s + dist) / (1000. * del_s))
    lonlats = gc.npts(lon1, lat1, lon2, lat2, npoints)
    lons = [lon1]
    lats = [lat1]
    xys = []
    for lon, lat in lonlats:
        # if we jump almost a whole globe away,
        # break the coords, resume afterwards
        if abs(lon - lons[-1]) > 170.:
            xys.append(m(lons, lats))
            lons = [lon]
            lats = [lat]
        else:
            lons.append(lon)
            lats.append(lat)
    xys.append(m(lons, lats))
    return [m.plot(x, y, **kwargs) for x, y in xys]


class CountryMap(Basemap):

    def __init__(self, shapefile=None, **kwargs):
        if shapefile is None:
            shapefile = join(
                dirname(abspath(__file__)),
                'data', 'tmwb-0.3',
                'TM_WORLD_BORDERS-0.3.shp'
            )
        if not exists(shapefile):
            raise ValueError('you must supply a valid country shapefile')
        if 'projection' not in kwargs:
            kwargs['projection'] = 'merc'
        self._isodict = data2isodict(shapefile2data(shapefile))
        for _, vs in self._isodict.items():
            for v in vs:
                if 'LAT' not in v or 'LON' not in v:
                    raise ValueError('we need lat/lon centroids for each country in the shapefile!')
        self._deaddict = {}
        super(CountryMap, self).__init__(**kwargs)

    def _isodata(self, iso):
        if (iso in self._isodict and
            len(self._isodict[iso]) > 1 and
            iso not in self._deaddict):
            parts = self._isodict[iso]
            lons = [d['LON'] for d in parts]
            lats = [d['LAT'] for d in parts]
            shapes = [d['SHAPE'][0] for d in parts]
            lon = sum(lons) / len(lons)
            lat = sum(lats) / len(lats)
            self._deaddict[iso] = {
                'LAT': lat,
                'LON': lon,
                'SHAPE': shapes
            }
        if iso in self._deaddict:
            return self._deaddict[iso]
        elif iso in self._isodict:
            return self._isodict[iso][0]
        else:
            raise ValueError('%s is not a valid country code on file!' % iso)

    def draw_gc(self, iso1, iso2, **kwargs):
        d1 = CountryMap._isodata(self, iso1)
        d2 = CountryMap._isodata(self, iso2)
        lon1, lat1 = d1['LON'], d1['LAT']
        lon2, lat2 = d2['LON'], d2['LAT']
        draw_gc(self, lon1, lat1, lon2, lat2, **kwargs)

    def draw_country(self, iso, color=None, edgecolor=None, linewidth=None, alpha=None):
        lines = shape2collection(self, CountryMap._isodata(self, iso)['SHAPE'])
        if color is not None:
            lines.set_facecolor(color)
        if edgecolor is not None:
            lines.set_edgecolor(edgecolor)
        if linewidth is not None:
            lines.set_linewidth(linewidth)
        if alpha is not None:
            lines.set_alpha(alpha)
        ax = self._check_ax()
        ax.add_collection(lines)
        self.set_axes_limits(ax=ax)
        return lines
