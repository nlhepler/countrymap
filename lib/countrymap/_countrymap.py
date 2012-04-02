
import json

from itertools import chain
from os.path import dirname, exists, join, abspath
from warnings import warn

import numpy as np

from matplotlib.collections import LineCollection, PolyCollection

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap.pyproj import Geod
from mpl_toolkits.basemap.shapefile import Reader as ShapeReader


__all__ = [
    'colors',
    'alphas',
    'linewidths',
    'enclaves',
    'shapefile2data',
    'data2isodict',
    'shape2collection',
    'drawgc',
    'CountryMap'
]


# based on google colors, for consistency
colors = {
    'park_color': '#c9dfaf',
    'land_color': '#f4f3f0',
    'lake_color': '#a5bfdd',
    'edge_color': 'k',
    'border_color': '#9d9d9d',
    'fade_color': '#c0c0bf' # '#d3d2d1' # '#c6c5c2'
}


linewidths = {
    'border_width': 0.28125,
    'edge_width': 0.28125,
}


alphas = {
    'edge_alpha': 0.75,
    'border_alpha': 1.
}


_dead_countries = {
    'CS': ['ME', 'RS'],
    'YU': ['BA', 'HR', 'ME', 'MK', 'RS', 'SI']
}


# completely locked within a single country
enclaves = [
    'LS',
    'SM',
    'VA'
]


def shapefile2data(shapefile, countries=False):
    shf = ShapeReader(shapefile)
    fields = [f[0] for f in shf.fields if isinstance(f, list)]
    latlons = None
    if countries and 'LAT' not in fields or 'LON' not in fields:
        with open(join(dirname(abspath(__file__)), 'data', 'latlons.json')) as fh:
            latlons = json.load(fh)
    records = []
    iso2s = set()
    for record in shf.shapeRecords():
        datadict = dict(zip(fields, record.record))
        if countries:
            iso2 = datadict['ISO_A2']
            iso2s.add(iso2)
            if latlons is not None and iso2 in latlons:
                lat, lon = latlons[iso2]
                datadict['LAT'] = lat
                datadict['LON'] = lon
        datadict['SHAPE'] = [record.shape]
        records.append(datadict)
    if countries:
        for iso2, latlon in latlons.items():
            if iso2 in iso2s:
                continue
            lat, lon = latlon
            datadict = { 'ISO_A2': iso2, 'LAT': lat, 'LON': lon }
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
        if 'ISO_A2' not in row:
            raise ValueError('no ISO_A2 field found, shapedata cannot be organized by ISO code!')
        isocode = row['ISO_A2']
        isodict[isocode] = [row]
        if isocode in new2dead:
            deadcodes = new2dead[isocode]
            for deadcode in deadcodes:
                if deadcode not in isodict:
                    isodict[deadcode] = []
                isodict[deadcode].append(row)
    return isodict


def interpolate(xy0, xy1, xval=None, yval=None):
    if (xval is None and yval is None) or (xval is not None and yval is not None):
        raise ValueError('must provide only one of xval or yval to interpolate the other')
    x0, y0 = xy0
    x1, y1 = xy1
    if xval is not None:
        return y0 + (xval - x0) * (y1 - y0) / (x1 - x0)
    else:
        return x0 + (yval - y0) * (x1 - x0) / (y1 - y0)


def shape2collection(m, shapes, collection=PolyCollection):
    segmentlist = []
    for shape in shapes:
        shapetype = shape.shapeType
        verts = shape.points
        if shapetype in (3, 5):
            parts = shape.parts.tolist()
            for idx1, idx2 in zip(parts, parts[1:] + [len(verts)]):
                lonlats = verts[idx1:idx2]
                for lon, lat in lonlats:
                    if (lon >  721. or
                        lon < -721. or
                        lat >   91. or
                        lat <  -91.):
                        raise ValueError('shape must have lat/lon vertices')
                # avoid shapes that cross axes
                # first, compile a list of segments
                # then join those segments together according to their boundary conditions
                lastll = lonlats[0]
                segs = []
                seg = [lastll]
                for ll in lonlats[1:]:
                    ew0, ns0 = lastll
                    ew1, ns1 = ll
                    if abs(ew0 - ew1) > 340.:
                        mod = 1
                        if ew0 < 0:
                            mod = -1
                        ewp = mod * 180
                        nsp = interpolate(lastll, ll, xval=ewp)
                        seg.append((ewp, nsp))
                        segs.append(seg)
                        ewp = mod * -180
                        # yp is the same both times
                        lastll = (ewp, nsp)
                        seg = [lastll]
                    seg.append(ll)
                    lastll = ll
                # make sure to append the last segment
                segs.append(seg)
                # if we have but one segment, we don't need to do anything else
                if len(segs) > 1:
                    # evaluate the junction of the start-end segments of segs
                    ll0 = segs[0][0]
                    ll1 = segs[-1][-1]
                    ew0, ns0 = ll0
                    ew1, ns1 = ll1
                    if abs(ew0 - ew1) > 340:
                        mod = 1
                        if ew0 < 0:
                            mod = -1
                        ewp = mod * 180
                        nsp = interpolate(ll0, ll1, xval=ewp)
                        segs[0].insert(0, (ewp, nsp))
                        ewp = mod * -180
                        # yp is the same both times
                        segs[-1].append((ewp, nsp))
                    else:
                        # don't span the boundary, so join!
                        segs[0].extend(segs.pop())
                    # if we have multiple segments, close them
                    if len(segs) > 1:
                        for i in range(len(segs))[::-1]:
                            ll0 = segs[i][0]
                            ll1 = segs[i][-1]
                            ew0, _ = ll0
                            ew1, _ = ll1
                            if abs(ew0 - ew1) < 2.:
                                segmentlist.append(segs[i] + [ll0])
                                del segs[i]
                        # now we would take care of joining segments together
                        # that fit in the same space, however this should be unnecessary
                        if len(segs):
                            warn("this code block should never be reached, so something bad happened")
                            pass
                    # this should never, ever happend, but just in case
                    else:
                        # antartica, TODO
                        pass
                # otherwise we have but a single region
                else:
                    segmentlist.extend(segs)
        else:
            raise ValueError('shape2collection only supports 2D shapes')
    coords = []
    for seg in segmentlist:
        lons, lats = list(zip(*seg))
        x, y = m(lons, lats)
        coords.append(list(zip(x, y)))
    return collection(coords, antialiaseds=(1,))


def drawgc(m, lon1, lat1, lon2, lat2, del_s=20., **kwargs):
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
        if abs(lon - lons[-1]) > 340.:
            xys.append(m(lons, lats))
            lons = [lon]
            lats = [lat]
        else:
            lons.append(lon)
            lats.append(lat)
    xys.append(m(lons, lats))
    return [m.plot(x, y, **kwargs) for x, y in xys]


class CountryMap(Basemap):

    def __init__(
        self,
        bordershapes=None,
        coastshapes=None,
        countryshapes=None,
        lakeshapes=None,
        landshapes=None,
        **kwargs):

        def getshapefile(shapefile):
            return join(
                dirname(abspath(__file__)),
                'data', 'ne', shapefile
            )

        def getshapesfromfile(shapefile):
            return [s for s in chain.from_iterable(r['SHAPE'] for r in shapefile2data(shapefile) if 'SHAPE' in r)]

        if bordershapes is None:
            bordershapes = getshapefile('ne_110m_admin_0_boundary_lines_land.shp')
        if coastshapes is None:
            coastshapes = getshapefile('ne_110m_coastline.shp')
        if countryshapes is None:
            countryshapes = getshapefile('ne_110m_admin_0_countries.shp')
        if lakeshapes is None:
            lakeshapes = getshapefile('ne_110m_lakes.shp')
        if landshapes is None:
            landshapes = getshapefile('ne_110m_land.shp')
        if not exists(countryshapes):
            raise ValueError('you must supply a valid country countryshapes')
        if 'projection' not in kwargs:
            kwargs['projection'] = 'merc'

        self._bordershapes = getshapesfromfile(bordershapes)
        self._coastshapes = getshapesfromfile(coastshapes)
        self._lakeshapes = getshapesfromfile(lakeshapes)
        self._landshapes = getshapesfromfile(landshapes)
        self._isodict = data2isodict(shapefile2data(countryshapes, countries=True))

        # fix sudan + south sudan
        self._isodict['SD'][0]['SHAPE'] += self._isodict['SS'][0]['SHAPE']

        # verify this data exists, we'll probably need it
        for iso2, vs in self._isodict.items():
            # silence these two countries in Natural Earth data,
            # we know they don't have centroids
            if iso2 in ('-99', 'SS'):
                continue
            for v in vs:
                if 'LAT' not in v or 'LON' not in v:
                    warn("missing LAT/LON centroid for country '%s'" % iso2)

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

    def drawgc(self, iso1, iso2, **kwargs):
        d1 = CountryMap._isodata(self, iso1)
        d2 = CountryMap._isodata(self, iso2)
        lon1, lat1 = d1['LON'], d1['LAT']
        lon2, lat2 = d2['LON'], d2['LAT']
        drawgc(self, lon1, lat1, lon2, lat2, **kwargs)

    def drawcoastlines(self, linewidth=None, color=None, antialiased=1, zorder=None):
        lines = shape2collection(self, self._coastshapes, collection=LineCollection)
        if color is not None:
            lines.set_edgecolor(color)
        if linewidth is not None:
            lines.set_linewidth(linewidth)
        if antialiased is not None:
            lines.set_antialiased(antialiased)
        if zorder is not None:
            lines.set_zorder(zorder)
        ax = self._check_ax()
        ax.add_collection(lines)
        self.set_axes_limits(ax=ax)
        return lines

    def fillcontinents(self, color=None, lake_color=None, zorder=None):
        land = shape2collection(self, self._landshapes)
        lakes = shape2collection(self, self._lakeshapes)
        land.set_linewidth(0.)
        lakes.set_linewidth(0.)
        if color is not None:
            land.set_color(color)
        if zorder is not None:
            land.set_zorder(zorder)
        ax = self._check_ax()
        if lake_color is None:
            lake_color = ax.get_axis_bgcolor()
        lakes.set_color(lake_color)
        lakes.set_zorder(land.get_zorder() + 1)
        ax.add_collection(land)
        ax.add_collection(lakes)
        self.set_axes_limits(ax=ax)
        return land

    def drawcountry(self, iso, color=None, edgecolor=None, linewidth=None, alpha=None):
        try:
            polys = shape2collection(self, CountryMap._isodata(self, iso)['SHAPE'])
        except KeyError:
            warn("country '%s' has no shape describing its borders" % iso)
            return None
        if color is not None:
            polys.set_facecolor(color)
        if edgecolor is not None:
            polys.set_edgecolor(edgecolor)
        if linewidth is not None:
            polys.set_linewidth(linewidth)
        if alpha is not None:
            polys.set_alpha(alpha)
        ax = self._check_ax()
        ax.add_collection(polys)
        self.set_axes_limits(ax=ax)
        return polys

    def drawmarkers(self, isos, s, c='b', marker='o', alpha=None, **kwargs):
        assert len(isos) == len(s), "isos and s must be of the same length"
        x, y = [np.zeros((len(isos),)) for _ in range(2)]
        for i, iso in enumerate(isos):
            d = CountryMap._isodata(self, iso)
            lon, lat = d['LON'], d['LAT']
            x[i], y[i] = self(lon, lat)
        return self.scatter(x, y, s, c, marker, alpha, **kwargs)
