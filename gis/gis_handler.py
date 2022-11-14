import geopandas as gpd
import rasterio
from shapely.geometry import Point

# precision of input coordinates, 5 decimals is about 1.1m of a commercial 'correcting' GPS unit
# this post is helpful: https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
precision = 5

class GISHandler:
    """A class to handle GIS raster data and optimizer points."""
        
    def __init__(self, files):
        self.rasters = {}
        self.points = gpd.GeoDataFrame(columns=['x', 'y', 'result', 'geometry'], geometry='geometry')
        
        for key, src in files.items():
            if key in self.rasters:
                print('raster {} already loaded!'.format(key))
            else:
                self.rasters[key] = rasterio.open(src)
    
    def query(self, x, y):
        x, y = self.clean(x, y)
                                    
        if not self.points.loc[(self.points.x==x) & (self.points.y==y)].empty:
            print('point exists, returning original data')
            return self.points.loc[(self.points.x==x) & (self.points.y==y)]
        
        conditions = {'x': x, 'y': y, 'geometry': Point(x, y)}
        
        for key, raster in self.rasters.items():
            index = raster.index(x, y)
            try:
                 conditions[key] = raster.read(1)[index] # yes, by default this only reads the first band, but this is probably okay
            except IndexError as error:
                print('failed to read {} raster: {}'.format(key, error))    
        
        self.points = self.points.append(conditions, ignore_index=True)
        return self.points.iloc[-1:]
    
    def record(self, x, y, value):
        x, y = self.clean(x, y)
        
        if not self.points.loc[(self.points.x==x) & (self.points.y==y)].empty:
            self.points.loc[(self.points.x==x) & (self.points.y==y), 'result'] = value
            return self.points.loc[(self.points.x==x) & (self.points.y==y)]
        
        conditions = {'x': x, 'y': y, 'geometry': Point(x, y), 'result': value}
        self.points = self.points.append(conditions, ignore_index=True)
        return self.points.iloc[-1:]       
        
    def clean(self, x, y):
        return round(x, precision), round(y, precision)