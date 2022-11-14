import geopandas
import rasterio
import rasterio.plot
import matplotlib.pyplot as plt
from shapely.geometry import Point

color_map = 'turbo'
crs = 4326

class GISHandler:
    
    def __init__(self):
        self.raster_data = {}
        self.sample_points = geopandas.GeoDataFrame()
    
    def load_raster_files(self, files):
        for key, src in files.items():
            if key in self.raster_data:
                raise Exception('raster {} already loaded!'.format(key))
            else:
                self.raster_data[key] = rasterio.open(
    
    def display_raster(self, key):
        raster = self.raster_data[key]
        fig, ax = plt.subplots()
        ax = rasterio.plot.show(raster.data, extent=raster.extent, ax=ax, cmap=color_map)
    
    def query_point(self, x, y, rasters=self.raster_data):
        if x in self.sample_points
        
        
        geometry = Point(x, y)
        
        point = geopandas.GeoDataFrame(id, geometry=geometry, crs=crs)
        value = [x for x in self.raster_data[key].data.sample([(x, y)])]
        point[key] = value
        self.sample_points = self.sample_points.append(point)
        return value[0]
    
    def record(x, y, value):
        return 0
        
class Raster:
    
    def __init__(self, name, src):
        self.name = name
        self.src = src
        
        self.data = rasterio.open(self.src)
        self.extent = [self.data.bounds[0], self.data.bounds[2], self.data.bounds[1], self.data.bounds[3]]