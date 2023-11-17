# Python GIS Handler
# SEA Lab at Cornell University, last updated: 11/14/23

# import necessary packages, namely geopandas and rasterio
import pandas as pd
import geopandas as gpd
import rasterio
from shapely.geometry import Point

precision = 5 # precision of input coordinates to prevent duplicate calls, 5dec~1.1m

class GISHandler:
    """A class to handle GIS data and optimizer points."""
        
    def __init__(self, conditions, conflicts, scope):
        """Initializes handler by creating a GeoDataFrame to store point measurement data and a dictionary of loaded raster files."""
        self.conditions = {}
        self.conflicts = {}
        self.scope = gpd.read_file(scope)
        # create initial geodataframe with anticipated fields
        self.points = gpd.GeoDataFrame(columns=['x', 'y', 'geometry', 'result', 'ok-conditions', 'ok-scope', 'ok-conflicts'], geometry='geometry')
    
        # check if raster already loaded as a condition, otherwise read source
        for key, src in conditions.items():
            if key in self.conditions:
                print('raster {} already loaded!'.format(key))
            else:
                self.conditions[key] = rasterio.open(src)
        
        # check if vector already loaded as a conflict, otherwise read source
        for key, src in conflicts.items():
            if key in self.conflicts:
                print('vector {} already loaded!'.format(key))
            else:
                self.conflicts[key] = gpd.read_file(src)
                
        # calculate boundaries of loaded files
        self.extent = self.extent()
    
    def query(self, x, y):
        """Gets condition data for a specified geography location (lon/lat), stores it in the GeoDataFrame, and returns the row."""
        x, y = self.coordinate(x, y)
                                
        # check if point has already been queried
        if not self.points.loc[(self.points.x==x) & (self.points.y==y)].empty:
            #print('point exists, returning original data')
            return self.points.loc[(self.points.x==x) & (self.points.y==y)]
        
        point = Point(x, y)
        conditions = {'x': x, 'y': y, 'geometry': point, 'ok-conditions': True, 'ok-scope': False, 'ok-conflicts': True}
        
        # check if point is offshore in desired scope
        for polygon in self.scope['geometry']:
            if point.intersects(polygon):
                conditions['ok-scope'] = True
        
        # check if point intersects with any of the conflicts loaded as vectors
        for key, vector in self.conflicts.items():
            for polygon in vector['geometry']:
                if point.intersects(polygon):
                    conditions['ok-conflicts'] = False
        
        # check if point exists for all condition datasets, and pull condition data for point
        for key, raster in self.conditions.items():
            index = raster.index(x, y)
            try:
                conditions[key] = raster.read(1)[index] # yes, by default this only reads the first band, but this is probably okay
            except IndexError:
                conditions[key] = 0
            if conditions[key] == 0: # assumes that zero is not a valid value, which isn't necessarily correct
                conditions['ok-conditions'] = False
        
        self.points = pd.concat([self.points, pd.DataFrame([conditions])], ignore_index=True)
        return self.points.iloc[-1:]
    
    def record(self, x, y, value):
        """Records a computed value from the optimizer to a geographic point, returns row recorded to."""
        x, y = self.coordinate(x, y)
        
        if not self.points.loc[(self.points.x==x) & (self.points.y==y)].empty:
            self.points.loc[(self.points.x==x) & (self.points.y==y), 'result'] = value
            return self.points.loc[(self.points.x==x) & (self.points.y==y)]
        
        conditions = {'x': x, 'y': y, 'geometry': Point(x, y), 'result': value}
        self.points.loc[len(self.points.index)] = conditions
        return self.points.iloc[-1:]
    
    def query_grid(self, x_min, x_max, y_min, y_max, xy_delta):
        """Iteratively builds a grid of points with desired bounds and resolution for brute-force analysis."""
        self.grid = gpd.GeoDataFrame(columns=['x', 'y', 'geometry', 'result', 'ok-conditions', 'ok-scope', 'ok-conflicts'], geometry='geometry')
        
        x = x_min
        y = y_min
        while x <= x_max:
            while y <= y_max:
                self.query(x, y)
                y += xy_delta
            y=y_min
            x += xy_delta
    
    def save(self, name):
        """Saves all loaded points and values into GIS format."""
        self.points.to_file(name, driver='GeoJSON')
        
    def load(self, name):
        """Loads previously saved points from GIS file."""
        self.points = gpd.read_file(name)
        
    def coordinate(self, x, y):
        """Rounds coordinates to given precision to prevent uneccessary duplication. In the future could handle projections."""
        return round(x, precision), round(y, precision)
    
    def extent(self):
        """Calculates largest rectangular extent that includes data from all loaded rasters."""
        
        extent = [-180, 180, -90, 90] # format: [W, E, S, N]
        for src in self.conditions.values():
            extent = [max(extent[0], src.bounds[0]), min(extent[1], src.bounds[2]), max(extent[2], src.bounds[1]), min(extent[3], src.bounds[3])]
            
        return extent