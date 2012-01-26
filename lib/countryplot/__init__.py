

from mpl_toolkits.basemap import Basemap, pyproj
from mpl_toolkits.basemap.reader import Reader as ShapeReader

__version__ = '0.0.1'


# based on google colors, for consistency
_colors = {
    'continents': '#f4f3f0',
    'lake_color': '#a5bfdd'
    'edge_color': 'k',
    'border_color': '#c6c5c2'
}


def shapefile2data(shapefile):
    shf = ShapeReader(shapefile)
    fields = [f[0] for f in shf.fields if isinstance(f, list)]
    records = []
    for record in shf.shapeRecords():
        datadict = dict(zip(fields, record.record))
        datadict['SHAPE'] = record.shape
        records.append(datadict)
    return records


def shape2collection(m, shape):
    shapetype = shape.shapeType
    verts = shape.points
    if shapetype in (3, 5):
        parts = shp.parts.tolist()
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
    npoints = int((0.5 * 1000 * del_s + dist) / (1000. * del_s))
    lonlats = gc.npts(lon1, lat1, lon2, lat2, npoints)
    lons = [lon1]
    lats = [lat1]
    prevlon = lon1
    x1, y1 = None, None
    for lon, lat in lonlats:
        # if we jump almost a whole globe away, 
        if abs(lon - prevlon) > 170.:
            x1, y1 = m(lons, lats)
            lons = [lon]
            lats = [lat]
        else:
            lons.append(lon)
            lats.append(lat)
    x2, y2 = m(lons, lats)
    ret = [m.plot(x2, y2, **kwargs)]
    if x1 is not None and y1 is not None:
        ret.insert(0, m.plot(x1, y1, **kwargs))
    return ret
